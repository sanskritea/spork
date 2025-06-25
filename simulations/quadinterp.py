import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.spatial.distance import cdist

class QuadraticInterpolator2D:
    """
    2D Quadratic Interpolator matching Mathematica's behavior.
    
    This class implements quadratic interpolation on a 2D grid similar to 
    Mathematica's Interpolation function with InterpolationOrder -> 2.
    """
    
    def __init__(self, x_coords, y_coords, z_values, method='spline'):
        """
        Initialize the 2D quadratic interpolator.
        
        Parameters:
        -----------
        x_coords : array-like
            1D array of x coordinates (must be sorted)
        y_coords : array-like  
            1D array of y coordinates (must be sorted)
        z_values : array-like
            2D array of z values with shape (len(y_coords), len(x_coords))
        method : str
            Interpolation method: 'spline' (default) or 'polynomial'
        """
        self.x_coords = np.asarray(x_coords)
        self.y_coords = np.asarray(y_coords)
        self.z_values = np.asarray(z_values)
        self.method = method
        
        # Validate inputs
        if len(self.x_coords.shape) != 1 or len(self.y_coords.shape) != 1:
            raise ValueError("x_coords and y_coords must be 1D arrays")
            
        if self.z_values.shape != (len(self.y_coords), len(self.x_coords)):
            raise ValueError(f"z_values shape {self.z_values.shape} doesn't match "
                           f"expected shape ({len(self.y_coords)}, {len(self.x_coords)})")
        
        # Check if coordinates are sorted
        if not np.all(np.diff(self.x_coords) >= 0):
            raise ValueError("x_coords must be sorted in ascending order")
        if not np.all(np.diff(self.y_coords) >= 0):
            raise ValueError("y_coords must be sorted in ascending order")
            
        # Initialize interpolator
        if method == 'spline':
            self._init_spline_interpolator()
        elif method == 'polynomial':
            self._init_polynomial_interpolator()
        else:
            raise ValueError("Method must be 'spline' or 'polynomial'")
    
    def _init_spline_interpolator(self):
        """Initialize using scipy's RectBivariateSpline with quadratic order."""
        # RectBivariateSpline uses kx, ky for spline degrees (2 = quadratic)
        self.interpolator = RectBivariateSpline(
            self.y_coords, self.x_coords, self.z_values, 
            kx=min(2, len(self.x_coords)-1), 
            ky=min(2, len(self.y_coords)-1),
            s=0  # No smoothing, interpolate exactly through data points
        )
    
    def _init_polynomial_interpolator(self):
        """Initialize polynomial-based interpolator for small grids."""
        # For polynomial method, we'll use local quadratic fits
        # This is more similar to Mathematica's approach for small datasets
        pass
    
    def __call__(self, x, y):
        """
        Evaluate the interpolator at given points.
        
        Parameters:
        -----------
        x : float or array-like
            x coordinates where to evaluate
        y : float or array-like  
            y coordinates where to evaluate
            
        Returns:
        --------
        float or array
            Interpolated values
        """
        x = np.asarray(x)
        y = np.asarray(y)
        
        if self.method == 'spline':
            return self._evaluate_spline(x, y)
        else:
            return self._evaluate_polynomial(x, y)
    
    def _evaluate_spline(self, x, y):
        """Evaluate using spline interpolator."""
        # Handle scalar inputs
        if x.ndim == 0 and y.ndim == 0:
            return float(self.interpolator(y, x))
        
        # Handle array inputs
        x_flat = x.flatten()
        y_flat = y.flatten()
        
        # Check bounds
        x_min, x_max = self.x_coords[0], self.x_coords[-1]
        y_min, y_max = self.y_coords[0], self.y_coords[-1]
        
        # Clip to bounds (Mathematica-like behavior)
        x_clipped = np.clip(x_flat, x_min, x_max)
        y_clipped = np.clip(y_flat, y_min, y_max)
        
        # Evaluate
        result = np.array([self.interpolator(y_val, x_val) 
                          for x_val, y_val in zip(x_clipped, y_clipped)])
        
        # Reshape to match input shape
        if x.ndim > 0 or y.ndim > 0:
            target_shape = np.broadcast_shapes(x.shape, y.shape)
            result = result.reshape(target_shape)
        
        return result
    
    def _evaluate_polynomial(self, x, y):
        """Evaluate using local polynomial fits."""
        x = np.asarray(x)
        y = np.asarray(y)
        
        # For polynomial method, we find the nearest grid points
        # and fit a local quadratic surface
        
        if x.ndim == 0 and y.ndim == 0:
            return self._local_quadratic_fit(float(x), float(y))
        
        # Handle array inputs
        x_flat = x.flatten()
        y_flat = y.flatten()
        
        result = np.array([self._local_quadratic_fit(x_val, y_val) 
                          for x_val, y_val in zip(x_flat, y_flat)])
        
        if x.ndim > 0 or y.ndim > 0:
            target_shape = np.broadcast_shapes(x.shape, y.shape)
            result = result.reshape(target_shape)
            
        return result
    
    def _local_quadratic_fit(self, x, y):
        """Fit a local quadratic polynomial around the point (x, y)."""
        # Find the nearest grid points
        x_idx = np.searchsorted(self.x_coords, x)
        y_idx = np.searchsorted(self.y_coords, y)
        
        # Get a 3x3 neighborhood if possible (for quadratic fit)
        x_start = max(0, min(x_idx - 1, len(self.x_coords) - 3))
        y_start = max(0, min(y_idx - 1, len(self.y_coords) - 3))
        
        x_end = min(len(self.x_coords), x_start + 3)
        y_end = min(len(self.y_coords), y_start + 3)
        
        # Extract local grid
        local_x = self.x_coords[x_start:x_end]
        local_y = self.y_coords[y_start:y_end]
        local_z = self.z_values[y_start:y_end, x_start:x_end]
        
        # Create coordinate matrices
        X, Y = np.meshgrid(local_x, local_y, indexing='xy')
        
        # Flatten for polynomial fitting
        X_flat = X.flatten()
        Y_flat = Y.flatten()
        Z_flat = local_z.flatten()
        
        # Set up polynomial basis (quadratic: 1, x, y, x^2, xy, y^2)
        A = np.column_stack([
            np.ones(len(X_flat)),
            X_flat,
            Y_flat,
            X_flat**2,
            X_flat * Y_flat,
            Y_flat**2
        ])
        
        # Solve for coefficients
        try:
            coeffs = np.linalg.lstsq(A, Z_flat, rcond=None)[0]
        except np.linalg.LinAlgError:
            # Fallback to nearest neighbor if fitting fails
            distances = np.sqrt((X_flat - x)**2 + (Y_flat - y)**2)
            nearest_idx = np.argmin(distances)
            return Z_flat[nearest_idx]
        
        # Evaluate polynomial at (x, y)
        result = (coeffs[0] + coeffs[1]*x + coeffs[2]*y + 
                 coeffs[3]*x**2 + coeffs[4]*x*y + coeffs[5]*y**2)
        
        return result
    
    def gradient(self, x, y):
        """
        Compute the gradient (partial derivatives) at given points.
        
        Returns:
        --------
        tuple of arrays
            (dz/dx, dz/dy) at the given points
        """
        if self.method == 'spline':
            x = np.asarray(x)
            y = np.asarray(y)
            
            if x.ndim == 0 and y.ndim == 0:
                dx = float(self.interpolator(y, x, dx=1, dy=0))
                dy = float(self.interpolator(y, x, dx=0, dy=1))
                return dx, dy
            
            # Handle arrays
            x_flat = x.flatten()
            y_flat = y.flatten()
            
            dx_vals = np.array([self.interpolator(y_val, x_val, dx=1, dy=0) 
                               for x_val, y_val in zip(x_flat, y_flat)])
            dy_vals = np.array([self.interpolator(y_val, x_val, dx=0, dy=1) 
                               for x_val, y_val in zip(x_flat, y_flat)])
            
            if x.ndim > 0 or y.ndim > 0:
                target_shape = np.broadcast_shapes(x.shape, y.shape)
                dx_vals = dx_vals.reshape(target_shape)
                dy_vals = dy_vals.reshape(target_shape)
            
            return dx_vals, dy_vals
        else:
            raise NotImplementedError("Gradient not implemented for polynomial method")

