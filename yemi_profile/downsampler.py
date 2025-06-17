#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def downsample_points(x, y, z, target_bins=(500, 500), reduce_method='mean'):
    """
    Downsample the (x, y, z) point cloud by binning into a 2D grid of size target_bins
    and computing one representative (x, y, z) per bin. This preserves the overall shape.

    Parameters:
    -----------
    x, y, z : 1D numpy arrays of length N (the original points)
    target_bins : tuple (nx_bins, ny_bins), the number of bins along x and y
    reduce_method : 'mean' or 'median', how to aggregate points in each bin

    Returns:
    --------
    x_ds, y_ds, z_ds : 1D numpy arrays of the downsampled points
    """
    nx, ny = target_bins
    x_edges = np.linspace(x.min(), x.max(), nx + 1)
    y_edges = np.linspace(y.min(), y.max(), ny + 1)

    # Digitize each point into a bin index (ix, iy)
    ix = np.searchsorted(x_edges, x, side='right') - 1
    iy = np.searchsorted(y_edges, y, side='right') - 1

    # Clip to valid indices
    ix = np.clip(ix, 0, nx - 1)
    iy = np.clip(iy, 0, ny - 1)

    sum_x = np.zeros((nx, ny), dtype=np.float64)
    sum_y = np.zeros((nx, ny), dtype=np.float64)
    sum_z = np.zeros((nx, ny), dtype=np.float64)
    count = np.zeros((nx, ny), dtype=np.int64)

    if reduce_method == 'median':
        bin_x_vals = [[[] for _ in range(ny)] for _ in range(nx)]
        bin_y_vals = [[[] for _ in range(ny)] for _ in range(nx)]
        bin_z_vals = [[[] for _ in range(ny)] for _ in range(nx)]
    else:
        bin_x_vals = bin_y_vals = bin_z_vals = None

    N = x.size
    for i in range(N):
        ixx = ix[i]
        iyy = iy[i]
        if reduce_method == 'median':
            bin_x_vals[ixx][iyy].append(x[i])
            bin_y_vals[ixx][iyy].append(y[i])
            bin_z_vals[ixx][iyy].append(z[i])
        else:
            sum_x[ixx, iyy] += x[i]
            sum_y[ixx, iyy] += y[i]
            sum_z[ixx, iyy] += z[i]
            count[ixx, iyy] += 1

    x_ds = []
    y_ds = []
    z_ds = []

    for ixx in range(nx):
        for iyy in range(ny):
            if reduce_method == 'median':
                if not bin_z_vals[ixx][iyy]:
                    continue
                x_rep = np.median(bin_x_vals[ixx][iyy])
                y_rep = np.median(bin_y_vals[ixx][iyy])
                z_rep = np.median(bin_z_vals[ixx][iyy])
            else:
                c = count[ixx, iyy]
                if c == 0:
                    continue
                x_rep = sum_x[ixx, iyy] / c
                y_rep = sum_y[ixx, iyy] / c
                z_rep = sum_z[ixx, iyy] / c

            x_ds.append(x_rep)
            y_ds.append(y_rep)
            z_ds.append(z_rep)

    return np.array(x_ds), np.array(y_ds), np.array(z_ds)


def draw_contour_map_downsampled(filename, reduced_file, num_contours=10,
                                 method='linear', grid_res=200, bin_x=500,
                                 bin_y=500, reduce_method='mean'):
    """
    Reads x, y, z from `filename`, downsamples to ~bin_x*bin_y points using a 2D grid,
    writes the downsampled points to `reduced_file`, then interpolates the downsampled
    points onto a grid and draws a contour map.

    - filename: path to whitespace-separated x y z file.
    - reduced_file: path where downsampled x y z will be written.
    - num_contours: number of contour levels.
    - method: interpolation method for griddata ('linear', 'nearest', 'cubic').
    - grid_res: resolution of the final interpolation grid (grid_res x grid_res).
    - bin_x, bin_y: number of bins to use when downsampling in x and y directions.
    - reduce_method: 'mean' or 'median' to aggregate points in each bin.
    """
    # Load all data
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    print(f"Original number of points: {x.size}")

    # Downsample
    x_ds, y_ds, z_ds = downsample_points(x, y, z,
                                         target_bins=(bin_x, bin_y),
                                         reduce_method=reduce_method)
    print(f"Reduced to {x_ds.size} points using {reduce_method} aggregation in a {bin_x}×{bin_y} grid.")

    # Write reduced points to file
    reduced_data = np.column_stack((x_ds, y_ds, z_ds))
    np.savetxt(reduced_file, reduced_data, fmt="%.6f %.6f %.6f")
    print(f"Downsampled points written to '{reduced_file}'")

    # Interpolate downsampled points onto a meshgrid for contouring
    xi = np.linspace(x_ds.min(), x_ds.max(), grid_res)
    yi = np.linspace(y_ds.min(), y_ds.max(), grid_res)
    XI, YI = np.meshgrid(xi, yi)

    ZI = griddata((x_ds, y_ds), z_ds, (XI, YI), method=method)

    # Plot
    plt.figure(figsize=(8, 6))
    contour_set = plt.contourf(XI, YI, ZI, levels=num_contours, cmap='viridis')
    plt.colorbar(contour_set, label='Z value')
    contour_lines = plt.contour(XI, YI, ZI, levels=num_contours, colors='k', linewidths=0.5)
    plt.clabel(contour_lines, inline=True, fontsize=8)

    plt.scatter(x_ds, y_ds, c='white', s=5, alpha=0.6, edgecolors='k', label='Downsampled points')
    plt.title("Contour Map (Downsampled)")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.legend(loc='upper right')
    plt.axis('equal')
    plt.grid(True)
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Read x y z from a large data file, downsample without losing shape, "
                    "write reduced points to file, and draw a contour map."
    )
    parser.add_argument(
        "input_file",
        help="Path to the whitespace-separated x y z data file."
    )
    parser.add_argument(
        "reduced_file",
        help="Path to output file for downsampled x y z points."
    )
    parser.add_argument(
        "-n", "--num_contours",
        type=int, default=10,
        help="Number of contour levels (default: 10)."
    )
    parser.add_argument(
        "-r", "--resolution",
        type=int, default=200,
        help="Grid resolution for interpolation (default: 200×200)."
    )
    parser.add_argument(
        "-m", "--method",
        choices=['linear', 'nearest', 'cubic'],
        default='linear',
        help="Interpolation method for griddata (default: linear)."
    )
    parser.add_argument(
        "--binx",
        type=int, default=500,
        help="Number of bins along X for downsampling (default: 500)."
    )
    parser.add_argument(
        "--biny",
        type=int, default=500,
        help="Number of bins along Y for downsampling (default: 500)."
    )
    parser.add_argument(
        "--reduce_method",
        choices=['mean', 'median'],
        default='mean',
        help="How to aggregate points in each bin: 'mean' or 'median' (default: mean)."
    )

    args = parser.parse_args()
    draw_contour_map_downsampled(
        args.input_file,
        args.reduced_file,
        num_contours=args.num_contours,
        method=args.method,
        grid_res=args.resolution,
        bin_x=args.binx,
        bin_y=args.biny,
        reduce_method=args.reduce_method
    )
