import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import (
    griddata, interp2d, RegularGridInterpolator, 
    RBFInterpolator, LinearNDInterpolator, 
    CloughTocher2DInterpolator
)
from scipy.spatial import Delaunay
import pandas as pd

def compare_interpolation_methods(x, y, z, test_points=None):
    """
    Compare different Python interpolation methods with focus on 
    matching Mathematica's Interpolation[] behavior
    
    Parameters:
    x, y: coordinate arrays or flattened coordinate points
    z: function values f(x,y) 
    test_points: optional test points for evaluation
    """
    z=z.T
    # Ensure we have the right data format
    if len(x.shape) == 1 and len(y.shape) == 1:
        # If x, y are 1D arrays, create meshgrid
        X, Y = np.meshgrid(x, y)
        points = np.column_stack([X.ravel(), Y.ravel()])
        values = z.ravel()
    else:
        # If x, y are already coordinate points
        points = np.column_stack([x.ravel(), y.ravel()])
        values = z.ravel()
    
    # Create test points if not provided
    if test_points is None:
        x_min, x_max = points[:, 0].min(), points[:, 0].max()
        y_min, y_max = points[:, 1].min(), points[:, 1].max()
        xi = np.linspace(x_min, x_max, 50)
        yi = np.linspace(y_min, y_max, 50)
        Xi, Yi = np.meshgrid(xi, yi)
        test_points = np.column_stack([Xi.ravel(), Yi.ravel()])
    
    interpolators = {}
    results = {}
    
    # Method 1: scipy.interpolate.griddata (closest to Mathematica default)
    print("Testing griddata methods...")
    for method in ['linear', 'cubic']:
        try:
            zi = griddata(points, values, test_points, method=method)
            interpolators[f'griddata_{method}'] = zi
            print(f"  ✓ griddata {method}: Success")
        except Exception as e:
            print(f"  ✗ griddata {method}: {e}")
    
    # Method 2: RegularGridInterpolator (if data is on regular grid)
    print("\nTesting RegularGridInterpolator...")
    try:
        if len(x.shape) == 1 and len(y.shape) == 1:
            # For regular grid data
            rgi = RegularGridInterpolator((x, y), z, method='linear')
            zi = rgi(test_points)
            interpolators['regular_grid_linear'] = zi
            
            rgi_cubic = RegularGridInterpolator((x, y), z, method='cubic')
            zi_cubic = rgi_cubic(test_points)
            interpolators['regular_grid_cubic'] = zi_cubic
            print("  ✓ RegularGridInterpolator: Success")
        else:
            print("  - RegularGridInterpolator: Skipped (irregular grid)")
    except Exception as e:
        print(f"  ✗ RegularGridInterpolator: {e}")
    
    # Method 3: RBF Interpolation (Radial Basis Functions)
    print("\nTesting RBF methods...")
    for kernel in ['linear', 'cubic', 'quintic', 'thin_plate_spline']:
        try:
            rbf = RBFInterpolator(points, values, kernel=kernel)
            zi = rbf(test_points)
            interpolators[f'rbf_{kernel}'] = zi
            print(f"  ✓ RBF {kernel}: Success")
        except Exception as e:
            print(f"  ✗ RBF {kernel}: {e}")
    
    # Method 4: Delaunay triangulation based methods
    print("\nTesting triangulation methods...")
    try:
        # Linear interpolation on Delaunay triangulation
        linear_interp = LinearNDInterpolator(points, values)
        zi_linear = linear_interp(test_points)
        interpolators['delaunay_linear'] = zi_linear
        
        # Clough-Tocher cubic interpolation
        cubic_interp = CloughTocher2DInterpolator(points, values)
        zi_cubic = cubic_interp(test_points)
        interpolators['clough_tocher'] = zi_cubic
        print("  ✓ Triangulation methods: Success")
    except Exception as e:
        print(f"  ✗ Triangulation methods: {e}")
    
    return interpolators, test_points

def export_for_mathematica_comparison(x, y, z, interpolators, test_points, filename_base="comparison"):
    """
    Export data in format suitable for Mathematica comparison
    """
    # Export original data
    if len(x.shape) == 1 and len(y.shape) == 1:
        X, Y = np.meshgrid(x, y)
        original_data = np.column_stack([X.ravel(), Y.ravel(), z.ravel()])
    else:
        original_data = np.column_stack([x.ravel(), y.ravel(), z.ravel()])
    
    np.savetxt(f"{filename_base}_original_data.csv", original_data, 
               delimiter=',', header='x,y,f', comments='')
    
    # Export test points
    np.savetxt(f"{filename_base}_test_points.csv", test_points, 
               delimiter=',', header='x,y', comments='')
    
    # Export interpolation results
    results_data = test_points.copy()
    method_names = []
    
    for method_name, values in interpolators.items():
        if values is not None:
            # Handle NaN values
            clean_values = np.nan_to_num(values, nan=0.0)
            results_data = np.column_stack([results_data, clean_values])
            method_names.append(method_name)
    
    header = 'x,y,' + ','.join(method_names)
    np.savetxt(f"{filename_base}_interpolation_results.csv", results_data, 
               delimiter=',', header=header, comments='')
    
    print(f"\nExported files:")
    print(f"  - {filename_base}_original_data.csv (for Mathematica Interpolation)")
    print(f"  - {filename_base}_test_points.csv (evaluation points)")
    print(f"  - {filename_base}_interpolation_results.csv (Python results)")

def analyze_interpolation_quality(interpolators, test_points):
    """
    Analyze and compare interpolation quality
    """
    print("\nInterpolation Quality Analysis:")
    print("=" * 50)
    
    for method_name, values in interpolators.items():
        if values is not None:
            valid_mask = ~np.isnan(values)
            valid_count = np.sum(valid_mask)
            total_count = len(values)
            
            if valid_count > 0:
                mean_val = np.mean(values[valid_mask])
                std_val = np.std(values[valid_mask])
                min_val = np.min(values[valid_mask])
                max_val = np.max(values[valid_mask])
                
                print(f"\n{method_name}:")
                print(f"  Valid points: {valid_count}/{total_count}")
                print(f"  Range: [{min_val:.6f}, {max_val:.6f}]")
                print(f"  Mean ± Std: {mean_val:.6f} ± {std_val:.6f}")
            else:
                print(f"\n{method_name}: No valid interpolated values")

def plot_comparison(x, y, z, interpolators, test_points, max_plots=6):
    """
    Create comparison plots of different interpolation methods
    """
    valid_methods = [(name, vals) for name, vals in interpolators.items() 
                    if vals is not None and not np.all(np.isnan(vals))]
    
    n_methods = min(len(valid_methods), max_plots)
    if n_methods == 0:
        print("No valid interpolation results to plot")
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.ravel()
    
    # Determine global color scale
    all_values = []
    for _, vals in valid_methods[:n_methods]:
        valid_vals = vals[~np.isnan(vals)]
        if len(valid_vals) > 0:
            all_values.extend(valid_vals)
    
    if len(all_values) > 0:
        vmin, vmax = np.min(all_values), np.max(all_values)
    else:
        vmin, vmax = 0, 1
    
    for i, (method_name, values) in enumerate(valid_methods[:n_methods]):
        ax = axes[i]
        
        # Reshape for plotting if needed
        if len(test_points) == len(values):
            # Create a regular grid for plotting
            xi = np.unique(test_points[:, 0])
            yi = np.unique(test_points[:, 1])
            if len(xi) * len(yi) == len(test_points):
                Xi, Yi = np.meshgrid(xi, yi)
                Zi = values.reshape(len(yi), len(xi))
                
                im = ax.contourf(Xi, Yi, Zi, levels=20, vmin=vmin, vmax=vmax)
                ax.set_title(f'{method_name}')
                ax.set_xlabel('x')
                ax.set_ylabel('y')
                plt.colorbar(im, ax=ax)
    
    # Hide unused subplots
    for i in range(n_methods, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.show()

# Example usage and Mathematica comparison workflow
def mathematica_comparison_workflow():
    """
    Complete workflow for comparing with Mathematica
    """
    print("Mathematica Interpolation Comparison Workflow")
    print("=" * 50)
    
    # Example: Create test data
    x = np.linspace(0, 5, 21)
    y = np.linspace(0, 3, 16)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(X) * np.cos(Y) + 0.1 * X * Y  # Example function
    
    print(f"Created test data: {len(x)} x {len(y)} grid")
    
    # Compare interpolation methods
    interpolators, test_points = compare_interpolation_methods(x, y, Z)
    
    # Analyze quality
    analyze_interpolation_quality(interpolators, test_points)
    
    # Export for Mathematica
    export_for_mathematica_comparison(x, y, Z, interpolators, test_points)
    
    # Create plots
    plot_comparison(x, y, Z, interpolators, test_points)
    
    print("\nNext steps for Mathematica comparison:")
    print("1. Import comparison_original_data.csv into Mathematica")
    print("2. Create interpolation: f = Interpolation[data]")
    print("3. Import test points and evaluate: f @@@ testPoints")
    print("4. Compare with Python results in comparison_interpolation_results.csv")
    
    return interpolators, test_points

'''
# Mathematica code template for comparison
mathematica_code = 
(* Mathematica comparison code *)
(* Import the data *)
data = Import["comparison_original_data.csv", "CSV", "HeaderLines" -> 1];
testPoints = Import["comparison_test_points.csv", "CSV", "HeaderLines" -> 1];

(* Create Mathematica interpolation *)
f = Interpolation[data];

(* Evaluate at test points *)
mmaResults = f @@@ testPoints;

(* Export results *)
Export["mathematica_results.csv", 
  Transpose[{testPoints[[All, 1]], testPoints[[All, 2]], mmaResults}],
  {"CSV", "HeaderLines" -> {"x", "y", "mathematica"}}];

(* Compare with Python results *)
pythonResults = Import["comparison_interpolation_results.csv", "CSV", "HeaderLines" -> 1];

(* Calculate differences for each Python method *)
nMethods = Length[pythonResults[[1]]] - 2; (* subtract x,y columns *)
differences = Table[
  Abs[mmaResults - pythonResults[[All, i + 2]]], 
  {i, nMethods}
];

(* Print comparison statistics *)
Print["Method comparison with Mathematica:"];
Do[
  Print[StringForm["Method ``: Mean error = ``, Max error = ``", 
    i, Mean[differences[[i]]], Max[differences[[i]]]]],
  {i, nMethods}
];
'''

def load_csv_data(x_file, y_file, z_file):
    """
    Load x, y, and f(x,y) data from separate CSV files
    
    Parameters:
    x_file: path to CSV file containing x coordinates
    y_file: path to CSV file containing y coordinates  
    z_file: path to CSV file containing f(x,y) values
    
    Returns:
    x, y, z arrays
    """
    print(f"Loading data from CSV files...")
    
    # Load the data
    x = np.loadtxt(x_file, delimiter=',').flatten()
    y = np.loadtxt(y_file, delimiter=',').flatten()
    z = np.loadtxt(z_file, delimiter=',')
    
    print(f"Loaded x: shape {x.shape}")
    print(f"Loaded y: shape {y.shape}")
    print(f"Loaded z: shape {z.shape}")
    
    # Handle different z formats
    if z.ndim == 1:
        # If z is flattened, reshape to grid
        if len(z) == len(x) * len(y):
            z = z.reshape(len(y), len(x))
            print(f"Reshaped z to: {z.shape}")
        else:
            print("Warning: z length doesn't match x*y. Using as-is.")
    
    return x, y, z

def csv_comparison_workflow(x_file, y_file, z_file, output_prefix="csv_comparison"):
    """
    Complete workflow for CSV data comparison with Mathematica
    
    Parameters:
    x_file: path to x coordinates CSV
    y_file: path to y coordinates CSV
    z_file: path to f(x,y) values CSV
    output_prefix: prefix for output files
    """
    print("CSV Data Interpolation Comparison Workflow")
    print("=" * 50)
    
    # Load data from CSV files
    try:
        x, y, z = load_csv_data(x_file, y_file, z_file)
    except Exception as e:
        print(f"Error loading CSV files: {e}")
        print("Make sure your CSV files exist and are properly formatted")
        return None, None
    
    # Compare interpolation methods
    interpolators, test_points = compare_interpolation_methods(x, y, z)
    
    # Analyze quality
    analyze_interpolation_quality(interpolators, test_points)
    
    # Export for Mathematica
    export_for_mathematica_comparison(x, y, z, interpolators, test_points, output_prefix)
    
    # Create plots
    plot_comparison(x, y, z, interpolators, test_points)
    
    print(f"\nFiles exported with prefix: {output_prefix}")
    print("Next steps for Mathematica comparison:")
    print(f"1. Import {output_prefix}_original_data.csv into Mathematica")
    print("2. Create interpolation: f = Interpolation[data]")
    print("3. Import test points and evaluate: f @@@ testPoints")
    print(f"4. Compare with Python results in {output_prefix}_interpolation_results.csv")
    
    return interpolators, test_points

# Alternative: Load from single CSV with x,y,z columns
def load_single_csv_data(filename, x_col=0, y_col=1, z_col=2):
    """
    Load x, y, z data from a single CSV file with columns
    
    Parameters:
    filename: path to CSV file
    x_col, y_col, z_col: column indices for x, y, z data
    """
    print(f"Loading data from single CSV file: {filename}")
    
    data = np.loadtxt(filename, delimiter=',', skiprows=1)  # Skip header if present
    
    x_points = data[:, x_col]
    y_points = data[:, y_col] 
    z_values = data[:, z_col]
    
    print(f"Loaded {len(x_points)} data points")
    print(f"X range: [{x_points.min():.3f}, {x_points.max():.3f}]")
    print(f"Y range: [{y_points.min():.3f}, {y_points.max():.3f}]")
    print(f"Z range: [{z_values.min():.3f}, {z_values.max():.3f}]")
    
    return x_points, y_points, z_values

def single_csv_comparison_workflow(filename, output_prefix="single_csv_comparison"):
    """
    Workflow for single CSV file with x,y,z columns
    """
    print("Single CSV Interpolation Comparison Workflow")
    print("=" * 50)
    
    try:
        x_points, y_points, z_values = load_single_csv_data(filename)
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        return None, None
    
    # For scattered data, we need to use the points directly
    interpolators, test_points = compare_interpolation_methods(x_points, y_points, z_values)
    
    # Analyze quality  
    analyze_interpolation_quality(interpolators, test_points)
    
    # Export for Mathematica
    export_for_mathematica_comparison(x_points, y_points, z_values, interpolators, test_points, output_prefix)
    
    print(f"\nFiles exported with prefix: {output_prefix}")
    
    return interpolators, test_points

if __name__ == "__main__":
    # Example usage for separate CSV files
    print("Choose your data format:")
    print("1. Separate CSV files for x, y, and f(x,y)")
    print("2. Single CSV file with x,y,z columns")
    print("3. Run example workflow")
    
    choice = input("Enter choice (1/2/3): ").strip()
    
    if choice == "1":
        x_file = input("Enter path to x coordinates CSV: ").strip()
        y_file = input("Enter path to y coordinates CSV: ").strip()
        z_file = input("Enter path to f(x,y) values CSV: ").strip()
        
        interpolators, test_points = csv_comparison_workflow(x_file, y_file, z_file)
        
    elif choice == "2":
        filename = input("Enter path to CSV file with x,y,z columns: ").strip()
        interpolators, test_points = single_csv_comparison_workflow(filename)
        
    else:
        # Run example workflow
        interpolators, test_points = mathematica_comparison_workflow()
    
    # Print Mathematica code
    print("\n" + "="*60)
    print("MATHEMATICA CODE FOR COMPARISON:")
    print("="*60)
    # print(mathematica_code)