'''
Script for generating object.grid.h5 files with NorSat-1 satellite model.
'''

import h5py
import numpy as np
import sys

def addSphere(space, radius, x, y, z):
    print(f"Added a sphere of radius {radius:.1f} at {x:1.f}, {y:1.f}, {z:1.f}")
    xs, ys, zs = np.indices(space.shape)
    return space | ((xs - x)**2 + (ys - y)**2 + (zs - z)**2 <= radius**2)

def addCylinder(space, radius, x, y, z, length, axis=0):
    if axis == 1:
        return addCylinder(space.swapaxes(0, 1), radius, y, x, z, length, axis=0).swapaxes(0, 1)
    elif axis == 2:
        return addCylinder(space.swapaxes(0, 2), radius, z, y, x, length, axis=0).swapaxes(0, 2)
    else:
        xs, ys, zs = np.indices(space.shape)
        return space | ((xs >= x) & (xs <= x + length) & ((ys - y) ** 2 + (zs - z) ** 2 <= radius ** 2))

def addCuboid(space, x, y, z, dx, dy, dz):
    print(f"Added a cuboid of dimensions {dx:.1f}, {dy:.1f}, {dz:.1f} at {x:.1f}, {y:.1f}, {z:.1f}")
    xs, ys, zs = np.indices(space.shape)
    return space | ((xs >= x) & (xs <= x + dx) & (ys >= y) & (ys <= y + dy) & (zs >= z) & (zs <= z + dz))

'''
NorSat-1 is defined in Emsis examples as:
 - Main body - cuboid stretching from (59,74,63) to (80,86,93),
 - Solar panel 1 - cuboid stretching from (59, 62, 63) to (60, 74, 93),
 - Solar panel 2 - cuboid stretching from (59, 86, 63) to (60, 98, 93),
 - Langmuir probe rod 1 - cylinder of radius 0.5, aligned with the y-axis, stretching from (80, 51, 63) to (80, 109, 63),
 - Langmuir probe rod 2 - cylinder of radius 0.5, aligned with the y-axis, stretching from (80, 51, 93) to (80, 109, 93),
 - Langmuir probe ball 1 - sphere of radius 0.5, centered at (80, 110, 63),
 - Langmuir probe ball 2 - sphere of radius 0.5, centered at (80, 110, 93),
 - Langmuir probe ball 3 - sphere of radius 0.5, centered at (80, 50, 63),
 - Langmuir probe ball 4 - sphere of radius 0.5, centered at (80, 50, 93)

Considering the length in x-axis as a normalising factor, we can shift everything by (-59, -80, -78) and divide all the dimensions by 21.
This results in:
 - Main body - cuboid stretching from (0, -2/7, -5/7) to (1, 2/7, 5/7),
 - Solar panel 1 - cuboid stretching from (0, -6/7, -5/7) to (1/21, -2/7, 5/7),
 - Solar panel 2 - cuboid stretching from (0, 2/7, -5/7) to (1/21, 6/7, 5/7),
 - Langmuir probe rod 1 - cylinder of radius 1/42, aligned with the y-axis, stretching from (1, -29/21, -5/7) to (1, 29/21, -5/7),
 - Langmuir probe rod 2 - cylinder of radius 1/42, aligned with the y-axis, stretching from (1, -29/21, 5/7) to (1, 29/21, 5/7),
 - Langmuir probe ball 1 - sphere of radius 1/42, centered at (1, 10/7, -5/7),
 - Langmuir probe ball 2 - sphere of radius 1/42, centered at (1, 10/7, 5/7),
 - Langmuir probe ball 3 - sphere of radius 1/42, centered at (1, -10/7, -5/7),
 - Langmuir probe ball 4 - sphere of radius 1/42, centered at (1, -10/7, 5/7)
'''


def addBody(space, length, front_x):
    x, y, z = space.shape
    space = addCuboid(space, front_x, -2/7 * length + y/2, -5/7 * length + z/2, length   , 4/7 * length, 10/7 * length)
    space = addCuboid(space, front_x, -6/7 * length + y/2, -5/7 * length + z/2, length/21, 4/7 * length, 10/7 * length)
    space = addCuboid(space, front_x,  2/7 * length + y/2, -5/7 * length + z/2, length/21, 4/7 * length, 10/7 * length)
    return space

def addLangmuirProbes(space, length, front_x):
    space = addCylinder(space, 1/42 * length, front_x + length, -29/21 * length + y/2, -5/7 * length + z/2, 58/21 * length)
    space = addCylinder(space, 1/42 * length, front_x + length, -29/21 * length + y/2,  5/7 * length + z/2, 58/21 * length)
    space = addSphere(space, 1/42 * length, front_x + length,  10/7 * length + y/2, -5/7 * length)
    space = addSphere(space, 1/42 * length, front_x + length,  10/7 * length + y/2,  5/7 * length)
    space = addSphere(space, 1/42 * length, front_x + length, -10/7 * length + y/2, -5/7 * length)
    space = addSphere(space, 1/42 * length, front_x + length, -10/7 * length + y/2,  5/7 * length)
    return space

if __name__ == "__main__":
    if len(sys.argv) < 8:
        print("Usage: python NorSat.py x y z length front_x mode filename")
        print("x, y, z - size of the simulation domain")
        print("length - length of the satellite in x-axis")
        print("front_x - x-coordinate of the front of the satellite")
        print("mode: 0 - only the main body, 1 - main body and Langmuir probes")
        print("filename - relative path of the output file")
        sys.exit(1)

    x, y, z = sys.argv[1:4]
    length = float(sys.argv[4])
    front_x = int(sys.argv[5])
    mode = int(sys.argv[6])
    filename = sys.argv[7]

    space = np.zeros((int(x), int(y), int(z)), dtype="int32")
    space = addBody(space, length, front_x)
    if mode == 1:
        space = addLangmuirProbes(space, length, front_x)


    file = h5py.File(filename, "w")
    file.create_dataset("Object", data=space, dtype="int32")