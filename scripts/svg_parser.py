#! /usr/bin/env python3 -B

import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET


def subdivide_bezier(points):
    result = [points[0]]

    def subdivide_polygon(polygon):
        def midpoint(a, b): return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)

        Q0 = midpoint(polygon[0], polygon[1])
        Q1 = midpoint(polygon[1], polygon[2])
        Q2 = midpoint(polygon[2], polygon[3])
        R0 = midpoint(Q0, Q1)
        R1 = midpoint(Q1, Q2)
        S = midpoint(R0, R1)
        return [Q0, R0, S, R1, Q2, polygon[3]]

    for i in range(0, len(points)-3, 3):
        # print(f'len {len(points)}: {i}, {i+4}')
        polygon = points[i:i+4]
        result += subdivide_polygon(polygon)

    return result


def bezier(points, num_subdivisions):
    result = []
    for _ in range(num_subdivisions):
        result = subdivide_bezier(points)
        points = result[:]
    return result[1:]


def draw_points(points):
    x = [p[0] for p in points]
    y = [p[1] for p in points]
    plt.scatter(x, y)
    plt.show()


def parse_points(points):
    points = [point.split(",") for point in points]
    points = [(float(x), float(y)) for x, y in points]
    return points


def parse_path_string(path_string, num_subdivisions):
    paths = path_string[1:-1].lower()
    paths = paths.split("zm")
    result_paths = []

    for path in paths:
        points = []
        tokens = path.split("l")
        for t in tokens:
            if "c" in t:
                chains = t.split("c")
                for c in chains:
                    if " " in c:
                        # Parsing bezier chain points
                        bezier_points = [points[-1]]
                        bezier_points += parse_points(c.split(" "))

                        assert(len(bezier_points) == 4)
                        bezier_points = bezier(bezier_points, num_subdivisions)
                        points += bezier_points
                    else:
                        # Parsing simple point
                        points += parse_points([c])
            else:
                # Parsing simple line point
                points += parse_points([t])

        result_paths.append(points[:-1])

    # for path in result_paths:
    #     draw_points(path)

    return result_paths


def create_json(infile, outfile, num_subdivisions):
    root = ET.parse(infile).getroot()
    data = {'screenspace': True, 'polygons': []}

    for element in root:
        if element.tag == "{http://www.w3.org/2000/svg}g":

            # Parsing translation matrix
            transform = element.attrib['transform'][7:-1].split(',')[-2:]
            translation = (float(transform[0]), float(transform[1]))

            for child in element:
                path_points = parse_path_string(
                    child.attrib['d'], num_subdivisions)
                for path in path_points:
                    for p in range(len(path)):
                        path[p] = (path[p][0] + translation[0],
                                   path[p][1] + translation[1])
                    data['polygons'].append(path)

        elif element.tag == "{http://www.w3.org/2000/svg}path":
            path_points = parse_path_string(
                element.attrib['d'], num_subdivisions)
            for path in path_points:
                data['polygons'].append(path)

    with open(outfile, "w") as out:
        json.dump(data, out, indent=2)


def main(infile, outfile, num_subdivisions):
    create_json(infile, outfile, num_subdivisions)


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    num_subdivisions = 2
    if len(sys.argv) > 3:
        num_subdivisions = int(sys.argv[3])

    create_json(infile, outfile, num_subdivisions)
