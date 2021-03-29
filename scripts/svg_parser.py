#! /usr/bin/env python3 -B

import sys
import json
import numpy as np
from xml.dom import minidom
import xml.etree.ElementTree as ET


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
                for point in path_points:
                    point[0] += translation[0]
                    point[1] += translation[1]
                    polygon.append(len(data['points']))
                    data['points'].append(point)

                data['polygons'].append(polygon)

        elif element.tag == "{http://www.w3.org/2000/svg}path":
            path_points = parse_path(element.attrib['d'])
            polygon = []
            for point in path_points:
                polygon.append(len(data['points']))
                data['points'].append(point)
            data['polygons'].append(polygon)

    with open(outfile, "w") as out:
        json.dump(data, out, indent=2)

    print(data)


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    create_json(infile, outfile)
