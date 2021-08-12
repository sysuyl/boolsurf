import glob
import json
import sys


def clean(dirname):
    files = glob.glob(f"{dirname}/*.json")
    print(len(files))

    for file in files:
        with open(file) as infile:
            js = json.load(infile)

        model = js['model']
        model = model.replace("./", "")
        model = model.replace("\\", "/")
        model = model.replace("//", "/")
        js['model'] = model

        with open(file, 'w') as outfile:
            json.dump(js, outfile, indent=2)


def riclean(dirname):
    files = glob.glob(f"{dirname}/*.json")
    print(len(files))

    for file in files:
        with open(file) as infile:
            js = json.load(infile)

        model = js['model']
        model = model[:-3]
        model += ".ply"
        js['model'] = model

        with open(file, 'w') as outfile:
            json.dump(js, outfile, indent=2)


if __name__ == "__main__":
    dirname = sys.argv[1]
    riclean(dirname)
