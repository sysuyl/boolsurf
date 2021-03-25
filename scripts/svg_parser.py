#! /usr/bin/env python3 -B

import sys
from xml.dom import minidom
import numpy as np


def parse(infile, outfile):
    doc = minidom.parse(infile) 
    path_strings = [path.getAttribute('d') for path
                    in doc.getElementsByTagName('path')]

    result = ""
    with open(outfile, "wb") as out:
        for path in path_strings:
            points = path[1:-2].lower().split("l")

            points1 = [point.split(",") for point in points]
            points1 = [(float(x), float(y)) for x, y in points1]

            np_points = np.array(points1, 'float32').flatten()
            np_points.tofile(out)

            result += str(len(points1)) + "\n"
            for x, y in points1:

                result += str(bin(x))+ " " + str(bin(y)) + " "
            result += "\n"
            
    print(result)

    # print(path_strings)
    doc.unlink()

if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    parse(infile, outfile)
