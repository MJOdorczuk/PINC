"""
** Uniform field creation for initial conditions **

@file                uniform.py
@author              Michał Jan Odorczuk <michaljo@uio.no>
@date                24.10.2024

"""


import h5py
import numpy as np
import sys

def uniform(x, y, z, value, filename):
    field = np.ones((x, y, z), dtype=np.float64) * value
    file = h5py.File(filename, "w")
    file.create_dataset("n=0.0", data=field, dtype="float64")
    file.flush()
    file.close()

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python uniform.py x y z value filename")
        print("x, y, z - size of the simulation domain")
        print("value - uniform value over the field")
        print("filename - relative path of the output file")
        sys.exit(1)
    x, y, z = [int(i) for i in sys.argv[1:4]]
    uniform(x, y, z, float(sys.argv[4]), sys.argv[5])