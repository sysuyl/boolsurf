import sys
import json

def main():
    test_filename0 = sys.argv[1]
    test_filename1 = sys.argv[2]
    output_filename = sys.argv[3]

    with open(test_filename0) as f: test0 = json.load(f)
    with open(test_filename1) as f: test1 = json.load(f)

    num_points = len(test0["points"])
    test0["points"] += test1["points"]

    for polygon in test1["polygons"]:
        for i in range(len(polygon)):
            polygon[i] += num_points

    test0["polygons"] += test1["polygons"]

    with open(output_filename, 'w') as json_file:
        json.dump(test0, json_file, indent=2)


if __name__ == '__main__':
    main()