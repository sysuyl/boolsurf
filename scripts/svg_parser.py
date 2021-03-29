#! /usr/bin/env python3 -B

import sys
import json
import numpy as np
#import matplotlib.pyplot as plt
from xml.dom import minidom
import xml.etree.ElementTree as ET


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
        S = midpoint(R0, R1)
        return [polygon[0], Q0, R0, S, R1, Q2]

    for i in range(0, len(points)-3, 3):
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


# def draw_points(points):
#     x = [p[0] for p in points]
#     y = [p[1] for p in points]
#     plt.scatter(x, y)
#     plt.show()


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
    root = ET.parse(infile).getroot()
    data = {'points_in_screenspace': True, 'points': [], 'polygons': []}

    def parse_path(path):
        points = path[1:-2].lower().split("l")
        points = [point.split(",") for point in points]
        points = [[float(x), float(y)] for x, y in points]
        points = points[:-1]
        return points

    for element in root:
        if element.tag == "{http://www.w3.org/2000/svg}g":

            transform = element.attrib['transform'][7:-1].split(',')[-2:]
            translation = (float(transform[0]), float(transform[1]))
            # inserire anche altre trasformazioni

            for child in element:
                polygon = []
                path_points = parse_path(child.attrib['d'])
                path_points = bezier(path_points, 2)
                for point in path_points:
                    point = (point[0] + translation[0],
                             point[1] + translation[1])
                    polygon.append(len(data['points']))
                    data['points'].append(point)

                data['polygons'].append(polygon)

        elif element.tag == "{http://www.w3.org/2000/svg}path":
            polygon = []
            path_points = parse_path(element.attrib['d'])
            path_points = bezier(path_points, 2)

            for point in path_points:
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
