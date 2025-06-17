#!/usr/bin/env python3
import argparse
import numpy as np


def transform_xyz(x, y, z, origin, theta_rad):
    """
    Translate and rotate points (x, y, z).

    - origin: (x0, y0, z0) to subtract.
    - theta_rad: rotation angle (radians) around Z axis (counterclockwise).
    """
    x0, y0, z0 = origin

    # 1) Translate so that (x0, y0, z0) → (0, 0, 0)
    xt = x - x0
    yt = y - y0
    zt = z - z0

    # 2) Rotate about Z:
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)
    x_r = xt * cos_t - yt * sin_t
    y_r = xt * sin_t + yt * cos_t
    z_r = zt  # Z unchanged by rotation around Z

    return x_r, y_r, z_r


def main():
    parser = argparse.ArgumentParser(
        description="Transform (x, y, z) by translating origin and rotating around Z."
    )
    parser.add_argument(
        "input_file",
        help="Path to input contour‐map data file (whitespace-separated x y z)."
    )
    parser.add_argument(
        "output_file",
        help="Path to output file for transformed x y z."
    )

    args = parser.parse_args()

    # Load original data (expects three columns: x y z)
    data = np.loadtxt(args.input_file)
    if data.ndim == 1 and data.size == 3:
        # Single line case
        data = data[np.newaxis, :]
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]

    #
    # don't change
    origin = (169714.4, 510127.1, -118.7)
    theta_rad = np.deg2rad(-20.3)

    # Apply transform
    x_new, y_new, z_new = transform_xyz(x, y, z, origin, theta_rad)

    # Stack columns and save
    transformed = np.column_stack((x_new, y_new, z_new))
    np.savetxt(args.output_file, transformed, fmt="%.6f %.6f %.6f")
    print(f"Transformed data written to '{args.output_file}'.")


if __name__ == "__main__":
    main()
