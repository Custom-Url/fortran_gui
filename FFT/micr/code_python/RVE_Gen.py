

import sys
sys.path.append('/mnt/data/dg765/FFT/PEDS/src')

import matplotlib.pyplot as plt
import vtk
import numpy as np

from peds.distributions_fibres import FibreRadiusDistribution, FibreDistribution2d



# Create a fibre radius distribution instance
fibre_dist = FibreRadiusDistribution(
    r_avg=3.5e-3,
    r_min=3e-3,
    r_max=4e-3,
    sigma=0.5e-3,
    gaussian=True,
    seed=42,
)


# Draw 5 sample radii
samples = fibre_dist.draw(5)
print("Sample fibre radii:", samples)

# Create a 2D fibre distribution instance
fibres_2d = FibreDistribution2d(
    n=1023,
    domain_size=250e-3,
    r_fibre_dist=fibre_dist,
    volume_fraction=0.45,
    seed=42,
    fast_code=True    #Set to False if you don't have the C++ code or compiler
)

# Generate one fibre distribution alpha matrix from the iterator
alpha = next(iter(fibres_2d))
print("Alpha matrix shape:", alpha.shape)
alpha_bin = (alpha == 0).astype(np.uint8)

# Plot the alpha matrix
plt.figure(figsize=(8, 6))
plt.imshow(alpha_bin, cmap='viridis', origin='lower')
plt.colorbar(label='Alpha Values')
plt.title('2D Fibre Distribution Alpha Matrix')
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.show()


def write_vtk_structured_points(filename, alpha2d, spacing=0.0001):
    """
    Save a 2D numpy array as a legacy VTK STRUCTURED_POINTS file with CELL_DATA for AMITEX_FFTP.
    Automatically infers dimensions and pads to ensure consistent cell-point structure.
    """
    # Cell dimensions (ny, nx) from alpha shape
    ny, nx = alpha2d.shape

    # Point dimensions: one more than cell dimensions
    nx_points = nx + 1
    ny_points = ny + 1
    nz_points = 2  # just two slices to match expected input format
    nz_cells = nz_points - 1

    n_cells = nx * ny * nz_cells

    print("Unique alpha values: ", np.unique(alpha2d))
    print("Alpha dtype: ", alpha2d.dtype)

    # Add fake Z layer to make shape (nz_cells, ny, nx) = (1, ny, nx)
    volume = np.reshape(alpha2d, (1, ny, nx)).astype(np.uint8)

    with open(filename, 'wb') as f:
        f.write(b"# vtk DataFile Version 4.5\n")
        f.write(b"2D Fibre Alpha for AMITEX\n")
        f.write(b"BINARY\n")
        f.write(b"DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx_points} {ny_points} {nz_points}\n".encode())
        f.write(b"ORIGIN 0.000 0.000 0.000\n")
        f.write(f"SPACING {spacing:.6e} {spacing:.6e} {spacing:.6e}\n".encode())
        f.write(f"CELL_DATA {n_cells}\n".encode())
        f.write(b"SCALARS geom unsigned_char\n")
        f.write(b"LOOKUP_TABLE default\n")
        f.write(volume.tobytes(order='C'))


write_vtk_structured_points("2D_Fibre_Distribution_Alpha.vtk", alpha_bin)


