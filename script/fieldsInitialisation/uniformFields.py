"""
** Set of uniform fields creation for initial conditions **

@file                uniformFields.py
@author              Michał Jan Odorczuk <michaljo@uio.no>
@date                24.10.2024

"""

from uniform import uniform
import json
import sys

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python uniformFields.py x y z directory fields")
        print("x, y, z - size of the simulation domain")
        print("directory - relative path of the directory for the fields to be saved")
        print("fields - json string containing field names and values")
        exit(1)
    x, y, z = [int(i) for i in sys.argv[1:4]]
    directory = sys.argv[4]
    fields = json.loads(sys.argv[5])
    for field in fields:
        uniform(x, y, z, fields[field], f"{directory}/{field}.grid.h5")
