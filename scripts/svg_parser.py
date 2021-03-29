#! /usr/bin/env python3 -B

import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from xml.dom import minidom

def subdivide_bezier(points):
    points.append(points[0])
    result = []
    def subdivide_polygon(polygon):
        def midpoint(a, b): return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        
        Q0 = midpoint(polygon[0], polygon[1])
        Q1 = midpoint(polygon[1], polygon[2])
        Q2 = midpoint(polygon[2], polygon[3])
        R0 = midpoint(Q0, Q1)
        R1 = midpoint(Q1, Q2)
        S  = midpoint(R0, R1)
        return [polygon[0], Q0, R0, S, R1, Q2]

    for i in range(0, len(points)-1, 3):
        print(f'len {len(points)}: {i}, {i+4}')
        polygon = points[i:i+4]
        print('len', len(polygon))
        result += subdivide_polygon(polygon)

    return result

def bezier(points, num_subdivisions):
    result = []
    for _ in range(num_subdivisions):
        result = subdivide_bezier(points)
        points = result[:]
    return result

def draw_points(points):
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    plt.scatter(x, y)
    plt.show()

def parse(infile, outfile):
    doc = minidom.parse(infile) 
    path_strings = [path.getAttribute('d') for path
                    in doc.getElementsByTagName('path')]

    with open(outfile, "wb") as out:
        out.write(len(path_strings).to_bytes(8, "little"))
        for path in path_strings:
            points = path[1:-2].lower().split("l")

            points = [point.split(",") for point in points]
            points = [(float(x), float(y)) for x, y in points]
            points = points[:-1]

            out.write(len(points).to_bytes(8, "little"))
            np_points = np.array(points, 'float32').flatten()
            np_points.tofile(out)

    doc.unlink()

def create_json(infile, outfile):
    doc = minidom.parse(infile) 
    path_strings = [path.getAttribute('d') for path
                    in doc.getElementsByTagName('path')]
    doc.unlink()

    data = {'points_in_screenspace': True, 'points': [], 'polygons' : []}
    for path in path_strings:
        polygon = []

        points = path[1:-2].lower().split("l")
        points = [point.split(",") for point in points]
        points = [[float(x), float(y)] for x, y in points]
        points = points[:-1]

        for point in points:
            polygon.append(len(data['points']))
            data['points'].append(point)
        
        data['polygons'].append(polygon)

    with open(outfile, "w") as out:
        json.dump(data, out, indent=2)

    print(data)

if __name__ == "__main__":
    points = [(0,0), (0,1), (1, 0.5), (1, 0), (1, 0.25), (0,0.5)]
    p = bezier(points, 4)
    draw_points(p)

    # infile = sys.argv[1]
    # outfile = sys.argv[2]
    # create_json(infile, outfile)
