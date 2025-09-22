import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append('/mnt/data/dg765/FFT/PEDS/src')
from peds.distributions_fibres import FibreRadiusDistribution, FibreDistribution2d

# ----------------------------- CONFIGURATION ----------------------------- #
output_dir = "unseen_LR"
os.makedirs(output_dir, exist_ok=True)

fibre_dist = FibreRadiusDistribution(
    r_avg=3.5e-3,
    r_min=3e-3,
    r_max=4e-3,
    sigma=0.5e-3,
    gaussian=True,
    seed=42,
)
volume_fractions = [0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65]
n_geometries = 10 
porosity_targets = [0.03, 0.06, 0.09, 0.12, 0.15, 0.18]

n = 63 
domain_size = 50e-3
spacing = domain_size / n
# ------------------------------------------------------------------------ #


def initialize_pores(n, alpha_bin, existing_pores):
    pore_point = np.zeros((n + 1, n + 1), dtype=np.uint8)
    candidates = np.argwhere((alpha_bin == 1) & (existing_pores == 0))
    if len(candidates) == 0:
        raise RuntimeError("No available pore seed positions remaining.")
    x, y = candidates[np.random.randint(len(candidates))]
    pore_point[x, y] = 2
    return pore_point


def grow_pores(pore_point, alpha_bin, momentum=3):
    n_plus_1 = pore_point.shape[0]
    pore_mask = np.zeros_like(pore_point, dtype=np.uint8)

    center = np.argwhere(pore_point == 2)
    if len(center) != 1:
        raise ValueError("pore_point must contain exactly one value of 2.")
    cx, cy = center[0]

    radius = 1
    touched_fibre = False
    post_touch_steps = 0

    while True:
        rr, cc = np.ogrid[:n_plus_1, :n_plus_1]
        dx = (rr - cx + n_plus_1) % n_plus_1
        dx = np.minimum(dx, n_plus_1 - dx)
        dy = (cc - cy + n_plus_1) % n_plus_1
        dy = np.minimum(dy, n_plus_1 - dy)
        dist = np.sqrt(dx ** 2 + dy ** 2)

        mask = (dist <= radius)
        pore_mask[:] = 0
        pore_mask[mask] = 2

        if not touched_fibre and np.any((pore_mask == 2) & (alpha_bin == 0)):
            touched_fibre = True

        if touched_fibre:
            post_touch_steps += 1
            if post_touch_steps >= momentum:
                break

        radius += 1

    return pore_mask


def write_vtk_structured_points(filename, alpha2d, spacing, pore_mask=None):
    ny, nx = alpha2d.shape
    nx_points, ny_points, nz_points = nx + 1, ny + 1, 2
    nz_cells = nz_points - 1
    n_cells = nx * ny * nz_cells

    if pore_mask is None:
        pore_mask = np.zeros((n + 1, n + 1), dtype=np.uint8)

    volume2d = alpha2d.copy().astype(np.uint8)
    volume2d += pore_mask
    volume2d[volume2d == 2] = 0 #Map Fibre/Pore Overalp to Fibre
    volume2d[volume2d == 1] = 2 # Map Matrix to 2 as expected by AMITEX
    volume2d[volume2d == 0] = 1 # Map Fibre to 1 as expected by AMITEX

    volume = np.reshape(volume2d, (1, ny, nx)).astype(np.uint8)

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

def write_vtk_mask(filename, alpha2d, pore_mask):
    ny, nx = alpha2d.shape
    mask2d = (~((alpha2d == 1) & (pore_mask == 2))).astype(np.uint8)  # 1 = simulate; 0 = pore

    mask3d = np.reshape(mask2d, (1, ny, nx)).astype(np.uint8)

    with open(filename, 'wb') as f:
        f.write(b"# vtk DataFile Version 4.5\n")
        f.write(b"Pore Mask for AMITEX\n")
        f.write(b"BINARY\n")
        f.write(b"DATASET STRUCTURED_POINTS\n")
        f.write(f"DIMENSIONS {nx+1} {ny+1} 2\n".encode())
        f.write(b"ORIGIN 0.000 0.000 0.000\n")
        f.write(f"SPACING {spacing:.6e} {spacing:.6e} {spacing:.6e}\n".encode())
        f.write(f"CELL_DATA {(nx)*(ny)}\n".encode())
        f.write(b"SCALARS mask unsigned_char\n")
        f.write(b"LOOKUP_TABLE default\n")
        f.write(mask3d.tobytes(order='C'))

# ------------------------------ MAIN LOOP ------------------------------- #
for vf in volume_fractions:
    for i in range(n_geometries):
        fibres_2d = FibreDistribution2d(
            n=n,
            domain_size=domain_size,
            r_fibre_dist=fibre_dist,
            volume_fraction=vf,
            seed=42 + i,
            fast_code=True
        )

        alpha = next(iter(fibres_2d))
        alpha_bin = (alpha == 0).astype(np.uint8)
        print("alpha succesful")

        pore_mask = np.zeros((n + 1, n + 1), dtype=np.uint8)
        total_cells = (n + 1) * (n + 1)
        target_porosity = 0.10
        max_attempts = 1000
        print("pore mask succesful")

        for target_porosity in porosity_targets:
            pore_mask = np.zeros((n + 1, n + 1), dtype=np.uint8)
            total_cells = (n + 1) * (n + 1)
            max_attempts = 1000


            for attempt in range(max_attempts):
                current_porosity = np.sum((pore_mask == 2) & (alpha_bin == 1))/ total_cells
                print("current porosity = ", current_porosity)
                if current_porosity >= target_porosity:
                    break

                try:
                    pore_point = initialize_pores(n, alpha_bin, pore_mask)
                except RuntimeError:
                    print("Warning: Ran out of valid pore seed locations.")
                    break

                new_pore = grow_pores(pore_point, alpha_bin, momentum=4)
                pore_mask = np.maximum(pore_mask, new_pore)

            final_porosity = np.sum((pore_mask == 2) & (alpha_bin ==1))/ total_cells
            print(f"Final porosity: {final_porosity:.3%}")

            vf_str = f"{vf:.2f}".replace('.', '')
            domain_str = f"{domain_size:.2f}"
            spacing_str = f"{spacing:.6f}"

            subdir = f"unseen_LR/L{domain_str}_vf{vf:.2f}/h{spacing_str}"
            os.makedirs(subdir, exist_ok=True)
            
            
            filename = os.path.join(subdir, f"iUC{i + 1}_vpMIN{target_porosity}.vtk")
            write_vtk_structured_points(filename, alpha_bin, spacing, pore_mask)
            print(f"Saved: {filename}")
            

