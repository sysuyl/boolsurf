#! /usr/bin/env python3 -B

import click
import sys
import glob
import os
import subprocess
import json
import random


@click.group()
def cli():
    pass


@cli.command()
@click.argument('meshes_dir')
@click.argument('svgs_dir')
@click.argument('output_jsons_dir')
def gui(meshes_dir, svgs_dir, output_jsons_dir):
    mesh_names = glob.glob(f'{meshes_dir}/meshes/*.ply')
    mesh_num = len(mesh_names)

    svgs_names = glob.glob(f'{svgs_dir}/*.svg')
    svgs_num = len(svgs_names)

    try:
        os.mkdir(output_jsons_dir)
    except:
        pass

    output_jsons_names = glob.glob(f'{output_jsons_dir}/*.json')
    output_jsons_names = [os.path.basename(filename).split(
        '.')[0] for filename in output_jsons_names]

    print(output_jsons_names)
    for mesh_id, mesh_name in enumerate(mesh_names):
        name = os.path.basename(mesh_name).split('.')[0]
        if name in output_jsons_names:
            continue

        msg = f'[{mesh_id+1}/{mesh_num}] {mesh_name}'
        print(msg + ' ' * max(0, 78-len(msg)))

        first = random.randint(0, svgs_num)
        second = random.randint(0, svgs_num)

        while (first != second):
            second = random.randint(0, svgs_num)

        first_svg = svgs_names[first]
        second_svg = svgs_names[second]

        cmd = f'.\\bin\\gui.exe --svg {first_svg} --output-test {output_jsons_dir}/{name}.json {mesh_name}'
        print(cmd)

        subprocess.run(cmd, shell=True).returncode


cli()
