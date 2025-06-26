import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, RegularGridInterpolator

def debug_interpolation_mismatch(x_file, y_file, z_file, mma_results_file):
    """
    Debug why Python and Mathematica interpolation results differ significantly
    """
    print("DEBUGGING INTERPOLATION MISMATCH")
    print("="*50)
    
    # Load your data
    x = np.loadtxt(x_file, delimiter=',').flatten()
    y = np.loadtxt(y_file, delimiter=',').flatten()
    z = np.loadtxt(z_file, delimiter=',')
    
    print(f"Original data shapes:")
    print(f"  X: {x.shape} - range [{x.min():.3f}, {x.max():.3f}]")
    print(f"  Y: {y.shape} - range [{y.min():.3f}, {y.max():.3f}]")
    print(f"  Z: {z.shape} - range [{z.min():.3f}, {z.max():.3f}]")
    
    # Check if Z needs reshaping
    if z.ndim == 1:
        if len(z) == len(x) * len(y):
            z = z.reshape(len(y), len(x))
            print(f"  Reshaped Z to: {z.shape}")
    
    # Check for dimension mismatch and fix
    X, Y = np.meshgrid(x, y)
    print(f"  Meshgrid X: {X.shape}, Y: {Y.shape}")
    
    if z.shape != X.shape:
        print(f"  ⚠️  Z shape {z.shape} doesn't match meshgrid {X.shape}")
        if z.shape == X.T.shape:
            print(f"  🔄 Transposing Z to match meshgrid")
            z = z.T
        elif z.shape == (len(x), len(y)):
            print(f"  🔄 Z appears to be (x,y) oriented, transposing to (y,x)")
            z = z.T
    
    print(f"  Final Z shape: {z.shape}")
    
    # Load Mathematica results
    try:
        mma_data = np.loadtxt(mma_results_file, delimiter=',', skiprows=1)
        test_x = mma_data[:, 0]
        test_y = mma_data[:, 1] 
        mma_values = mma_data[:, 2]
        print(f"\nMathematica results:")
        print(f"  Test points: {len(test_x)}")
        print(f"  Value range: [{mma_values.min():.3f}, {mma_values.max():.3f}]")
    except:
        # Create test points if no Mathematica file
        print("\nCreating test points (no Mathematica file found)")
        xi = np.linspace(x.min(), x.max(), 20)
        yi = np.linspace(y.min(), y.max(), 15)
        Xi, Yi = np.meshgrid(xi, yi)
        test_x = Xi.ravel()
        test_y = Yi.ravel()
        mma_values = None
    
    test_points = np.column_stack([test_x, test_y])
    
    # Test different data orientations
    print("\nTesting different data orientations:")
    
    orientations = [
        ("Original Z", z),
        ("Transposed Z", z.T),
        ("Flipped Y-axis", np.flipud(z)),
        ("Flipped X-axis", np.fliplr(z)),
        ("Both flipped", np.flipud(np.fliplr(z)))
    ]
    
    results = {}
    
    for name, z_test in orientations:
        try:
            # Method 1: RegularGridInterpolator (if regular grid)
            if z_test.shape == (len(y), len(x)):
                rgi = RegularGridInterpolator((y, x), z_test, 
                                            method='linear', bounds_error=False, fill_value=np.nan)
                # Note: RGI expects (y,x) order for coordinates
                values_rgi = rgi(np.column_stack([test_y, test_x]))
                results[f"{name}_RGI"] = values_rgi
            
            # Method 2: griddata
            X, Y = np.meshgrid(x, y)
            points = np.column_stack([X.ravel(), Y.ravel()])
            values = z_test.ravel()
            
            values_grid = griddata(points, values, test_points, method='linear')
            results[f"{name}_griddata"] = values_grid
            
            print(f"  ✓ {name}: Success")
            
        except Exception as e:
            print(f"  ✗ {name}: {e}")
    
    # Compare with Mathematica if available
    if mma_values is not None:
        print(f"\nComparison with Mathematica:")
        print(f"{'Method':<25} {'Mean Error':<15} {'Max Error':<15} {'Valid Points':<15}")
        print("-" * 70)
        
        for method_name, py_values in results.items():
            if py_values is not None:
                valid_mask = ~np.isnan(py_values)
                if np.sum(valid_mask) > 0:
                    valid_py = py_values[valid_mask]
                    valid_mma = mma_values[valid_mask]
                    
                    if len(valid_mma) > 0:
                        errors = np.abs(valid_py - valid_mma)
                        mean_error = np.mean(errors)
                        max_error = np.max(errors)
                        n_valid = len(valid_py)
                        
                        print(f"{method_name:<25} {mean_error:<15.2f} {max_error:<15.2f} {n_valid:<15}")
    
    # Visual comparison
    print(f"\nCreating visual comparison...")
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.ravel()
    
    # Plot original data - ensure X, Y, z are compatible
    try:
        im0 = axes[0].contourf(X, Y, z, levels=20)
        axes[0].set_title('Original Data')
        axes[0].set_xlabel('X')
        axes[0].set_ylabel('Y')
        plt.colorbar(im0, ax=axes[0])
        print("  ✓ Original data plot successful")
    except Exception as e:
        print(f"  ✗ Original data plot failed: {e}")
        axes[0].text(0.5, 0.5, f'Plot Error:\n{str(e)[:50]}...', 
                    transform=axes[0].transAxes, ha='center')
        axes[0].set_title('Original Data (Error)')
    
    # Plot best Python results
    plot_idx = 1
    for method_name, values in list(results.items())[:5]:
        if values is not None and plot_idx < 6:
            # Create grid for plotting
            xi = np.unique(test_x)
            yi = np.unique(test_y)
            
            if len(xi) * len(yi) == len(test_x):
                try:
                    Xi, Yi = np.meshgrid(xi, yi)
                    Zi = values.reshape(len(yi), len(xi))
                    valid_mask = ~np.isnan(Zi)
                    
                    if np.any(valid_mask):
                        # Replace NaN with interpolation for plotting
                        Zi_plot = np.where(np.isnan(Zi), np.nanmean(Zi), Zi)
                        im = axes[plot_idx].contourf(Xi, Yi, Zi_plot, levels=20)
                        axes[plot_idx].set_title(f'{method_name}')
                        axes[plot_idx].set_xlabel('X')
                        axes[plot_idx].set_ylabel('Y')
                        plt.colorbar(im, ax=axes[plot_idx])
                        print(f"  ✓ {method_name} plot successful")
                    else:
                        axes[plot_idx].text(0.5, 0.5, 'No valid data', 
                                          transform=axes[plot_idx].transAxes, ha='center')
                        axes[plot_idx].set_title(f'{method_name} (No data)')
                        print(f"  - {method_name}: No valid data")
                except Exception as e:
                    axes[plot_idx].text(0.5, 0.5, f'Plot error:\n{str(e)[:30]}...', 
                                      transform=axes[plot_idx].transAxes, ha='center')
                    axes[plot_idx].set_title(f'{method_name} (Error)')
                    print(f"  ✗ {method_name} plot failed: {e}")
            else:
                axes[plot_idx].text(0.5, 0.5, 'Irregular grid', 
                                  transform=axes[plot_idx].transAxes, ha='center')
                axes[plot_idx].set_title(f'{method_name} (Irregular)')
                print(f"  - {method_name}: Irregular grid, can't plot")
            
            plot_idx += 1
    
    # Hide unused subplots
    for i in range(plot_idx, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.show()
    
    return results

def check_coordinate_system(x, y, z):
    """
    Check if coordinate system matches expectations
    """
    print("COORDINATE SYSTEM CHECK")
    print("="*30)
    
    X, Y = np.meshgrid(x, y)
    
    print(f"Grid dimensions: {len(x)} x {len(y)}")
    print(f"Z matrix shape: {z.shape}")
    print(f"Expected shape: ({len(y)}, {len(x)}) [rows=y, cols=x]")
    
    if z.shape != (len(y), len(x)):
        print("⚠️  WARNING: Z shape doesn't match expected (y, x) convention")
        print("   Try transposing Z matrix: z = z.T")
    
    # Check corner values
    print(f"\nCorner value check:")
    print(f"  Bottom-left  [y[0], x[0]]  = z[0,0]   = {z[0,0]:.3f}")
    print(f"  Bottom-right [y[0], x[-1]] = z[0,-1]  = {z[0,-1]:.3f}")
    print(f"  Top-left     [y[-1], x[0]] = z[-1,0]  = {z[-1,0]:.3f}")
    print(f"  Top-right    [y[-1], x[-1]]= z[-1,-1] = {z[-1,-1]:.3f}")

# Mathematica debugging code
mathematica_debug_code = '''
(* MATHEMATICA DEBUGGING CODE *)

(* 1. Check your data orientation *)
data = Import["your_original_data.csv", "CSV", "HeaderLines" -> 1];
Print["Data dimensions: ", Dimensions[data]];
Print["X range: ", MinMax[data[[All, 1]]]];
Print["Y range: ", MinMax[data[[All, 2]]]];
Print["Z range: ", MinMax[data[[All, 3]]]];

(* 2. Try different interpolation methods *)
f1 = Interpolation[data]; (* Default *)
f2 = Interpolation[data, Method -> "Spline"];
f3 = Interpolation[data, InterpolationOrder -> 1]; (* Linear *)
f4 = Interpolation[data, InterpolationOrder -> 3]; (* Cubic *)

(* 3. Test at a few points *)
testPoint = {Mean[data[[All, 1]]], Mean[data[[All, 2]]]};
Print["Test point: ", testPoint];
Print["Default: ", f1 @@ testPoint];
Print["Spline: ", f2 @@ testPoint];
Print["Linear: ", f3 @@ testPoint];
Print["Cubic: ", f4 @@ testPoint];

(* 4. Check interpolation domain *)
Print["Interpolation domain: ", f1["Domain"]];

(* 5. Export results with different methods *)
testPoints = Import["your_test_points.csv", "CSV", "HeaderLines" -> 1];
results = Table[{
  testPoints[[i, 1]], 
  testPoints[[i, 2]], 
  f1 @@ testPoints[[i]], 
  f2 @@ testPoints[[i]], 
  f3 @@ testPoints[[i]], 
  f4 @@ testPoints[[i]]
}, {i, Length[testPoints]}];

Export["mathematica_methods_comparison.csv", results, 
  {"CSV", "HeaderLines" -> {"x", "y", "default", "spline", "linear", "cubic"}}];
'''

if __name__ == "__main__":
    print("Enter your CSV file paths:")
    x_file = input("X coordinates CSV: ").strip()
    y_file = input("Y coordinates CSV: ").strip() 
    z_file = input("Z values CSV: ").strip()
    mma_file = input("Mathematica results CSV (optional, press Enter to skip): ").strip()
    
    if not mma_file:
        mma_file = None
    
    results = debug_interpolation_mismatch(x_file, y_file, z_file, mma_file)
    
    print("\n" + "="*60)
    print("MATHEMATICA DEBUGGING CODE:")
    print("="*60)
    print(mathematica_debug_code)