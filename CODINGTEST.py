import numpy as np
from scipy.interpolate import griddata

# Example data points
x = np.array([1, 2, 3])
y = np.array([4, 5, 6])
zi = np.array([10, 20, 30])

# Example target grid points
x_grid = np.linspace(0, 5, 10)
y_grid = np.linspace(0, 5, 10)

# Interpolation of data using griddata
zi_grid = griddata(x, y, zi, x_grid, y_grid, method='linear')

# Print the interpolated data
print(zi_grid)