import sys
import json

def main():
    test_filename = sys.argv[1]
    model_filename = sys.argv[2]

    with open(test_filename) as f:
        js = json.load(f)

    splines = js['splines']

    result = {}
    result["points"] = []
    result["polygons"] = [[]]
    result["model"] = model_filename
    for spline in splines:
        polygon = []
        for point in spline['control_points']:
            point_id = len(result["points"])
            result["points"].append(point)
            polygon.append(point_id)
        
        if len(polygon) > 2:
            result['polygons'].append(polygon)
    
    with open('deleteme.json', 'w') as json_file:
        json.dump(result, json_file, indent=2)


if __name__ == '__main__':
    main()