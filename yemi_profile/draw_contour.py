#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata


def draw_contour_from_xyz(filename, num_levels=10, interp_method='linear', grid_res=200, pdf_file="contour_plot.pdf"):
    """
    Reads x, y, z from `filename` (whitespace-separated columns),
    interpolates onto a grid, and draws a contour map with `num_levels` contours.
    """
    # 1) Load data
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    # 2) Build a regular grid over the x,y range
    xi = np.linspace(x.min(), x.max(), grid_res)
    yi = np.linspace(y.min(), y.max(), grid_res)
    XI, YI = np.meshgrid(xi, yi)

    # 3) Interpolate z values onto that grid
    ZI = griddata((x, y), z, (XI, YI), method=interp_method)

    # 4) Plot filled contours and contour lines
    plt.figure(figsize=(16, 12))
    cf = plt.contourf(XI, YI, ZI, levels=num_levels, cmap='coolwarm')
    cbar = plt.colorbar(cf, pad=0.02)
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label(label='Height [m]', fontsize=16)

    cl = plt.contour(XI, YI, ZI, levels=num_levels, colors='k', linewidths=0.5)
    plt.clabel(cl, inline=True, fontsize=8)

    # 5) Overlay the original points (optional)
    #plt.scatter(x, y, c='white', s=5, alpha=0.6,
    #            edgecolors='k', label='Data points')

    plt.xlabel("X [m]", fontsize=16)
    plt.ylabel("Y [m]", fontsize=16)
    #plt.title("Contour Map from (x, y, z) Data", fontsize=14)

    # 6) Force 1:1 aspect ratio
    plt.gca().set_aspect('equal', adjustable='box')

    #plt.legend(loc='upper right', fontsize=10)
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.6)
    plt.tight_layout()

    plt.savefig(pdf_file, format="pdf")
    
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Draw a contour map from a file of x y z points."
    )
    parser.add_argument(
        "input_file",
        help="Path to the whitespace-separated x y z data file."
    )
    parser.add_argument(
        "-n", "--num_levels",
        type=int, default=10,
        help="Number of contour levels (default: 10)."
    )
    parser.add_argument(
        "-m", "--method",
        choices=['linear', 'nearest', 'cubic'],
        default='linear',
        help="Interpolation method for griddata (default: linear)."
    )
    parser.add_argument(
        "-r", "--resolution",
        type=int, default=200,
        help="Grid resolution (grid will be resolution×resolution, default: 200)."
    )
    parser.add_argument(
        "-p", "--pdf_file",
        type=str, default="contour_plot.pdf",
        help="Print contour plot to pdf file."
    )    

    args = parser.parse_args()
    draw_contour_from_xyz(
        args.input_file,
        num_levels=args.num_levels,
        interp_method=args.method,
        grid_res=args.resolution,
        pdf_file= args.pdf_file
    )
