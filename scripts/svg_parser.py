from svg.path import parse_path
from svg.path.path import Line, CubicBezier, Close, Move
from xml.dom import minidom
import numpy as np
from typing import Iterable
import sys
import matplotlib.pyplot as plt
import json

def to_xy(point):
    return (point.real, point.imag)


def flatten(items):
    """Yield items from any nested iterable; see Reference."""
    for x in items:
        if isinstance(x, Iterable) and not isinstance(x, (str, bytes)):
            for sub_x in flatten(x):
                yield sub_x
        else:
            yield x

def subdivide_bezier(points):
    result = [points[0]]

    def subdivide_polygon(polygon):
        def midpoint(a, b):
            return ((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)

        Q0 = midpoint(polygon[0], polygon[1])
        Q1 = midpoint(polygon[1], polygon[2])
        Q2 = midpoint(polygon[2], polygon[3])
        R0 = midpoint(Q0, Q1)
        R1 = midpoint(Q1, Q2)
        S = midpoint(R0, R1)
        return [Q0, R0, S, R1, Q2, polygon[3]]

    for i in range(0, len(points) - 3, 3):
        # print(f'len {len(points)}: {i}, {i+4}')
        polygon = points[i: i + 4]
        result += subdivide_polygon(polygon)

    return result


def bezier(points, num_subdivisions):
    result = []
    for _ in range(num_subdivisions):
        result = subdivide_bezier(points)
        points = result[:]
    return result

def transform_points(points, transform):
    conv_points = []
    for point in points:
        nppoint = np.array(point)
        nppoint = np.append(nppoint, 1)
        nppoint = transform @ nppoint

        conv_points.append(tuple(nppoint.tolist()[:2]))
    return conv_points


def normalize(conv_group_paths, scale=1.0):
    max_value = max(flatten(conv_group_paths))

    for cgp_id, conv_group_path in enumerate(conv_group_paths):
        for cp_id, conv_paths in enumerate(conv_group_path):
            for p_id, path in enumerate(conv_paths):
                for id, point in enumerate(path):
                    conv_point = ((point[0] / max_value * 2 - 1) * scale, (point[1] / max_value * 2 - 1) * scale)
                    conv_group_paths[cgp_id][cp_id][p_id][id] = conv_point

    return conv_group_paths


def parse_svg(svg_file, num_subdivisions):
    doc = minidom.parse(svg_file)

    shapes = []
    for group in doc.getElementsByTagName('g'):

        transform = group.getAttribute('transform')
        if transform == "": continue
        matrix = transform[7:-1].split(",")

        rot_transform = np.array(matrix, dtype=float).reshape(3, 2).transpose()
        rot_transform = np.append(rot_transform, [[0, 0, 1]], axis=0)

        for path_group in group.getElementsByTagName('path'):
            shape = []

            current_polygon = []

            path = parse_path(path_group.getAttribute('d'))
            num_paths = len(path)

            for segment in path:
                if isinstance(segment, Line):
                    points = [to_xy(segment.start), to_xy(segment.end)]
                    points = transform_points(points, rot_transform)
                    current_polygon += points

                elif isinstance(segment, CubicBezier):
                    points = [to_xy(segment.start), to_xy(segment.control1), to_xy(segment.control2), to_xy(segment.end)]
                    points = bezier(points, num_subdivisions)
                    points = transform_points(points, rot_transform)

                    if len(current_polygon) == 0: current_polygon += points
                    else: current_polygon += points[1:]


                elif isinstance(segment, Close):
                    if len(current_polygon) > 0:
                        shape.append(current_polygon)
                        current_polygon = []

                    continue

            if len(current_polygon) > 0:
                shape.append(current_polygon)
            
            shapes.append(shape)

    #shapes = normalize(shapes)
    return shapes


def draw_points(points):
    x = []
    y = []

    for group_paths in points:
        for path in group_paths:
            for segment in path:
                x.append(segment[0])
                y.append(segment[1])

    plt.scatter(x, y)
    plt.show()


def create_json(outfile, shapes_paths):
    data = {"screenspace": True, "shapes": [], "polygons": [], "are_closed":[]}

    for shape in shapes_paths:
        shape_ids = []
        for polygon in shape:
            shape_ids.append(len(data["polygons"]))

            if (polygon[0]==polygon[-1]):
                data["polygons"].append(polygon[:-1])
                data["are_closed"].append(True)
            else:
                data["polygons"].append(polygon)
                data["are_closed"].append(False)

        data["shapes"].append(shape_ids)

    with open(outfile, "w") as out:
        json.dump(data, out, indent=2)


if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    num_subdivisions = 2

    if len(sys.argv) > 3:
        num_subdivisions = int(sys.argv[3])

    shape_paths = parse_svg(infile, num_subdivisions)

    create_json(outfile, shape_paths)
