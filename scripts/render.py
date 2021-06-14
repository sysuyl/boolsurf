#! /usr/bin/env python3 -B

import click
import glob
import os
import subprocess


@click.group()
def cli():
    pass


@cli.command()
@click.argument('bin')
@click.argument('dirname')
@click.argument('spp', default=1, type=int)
def render(bin, dirname, spp):
    scene_dirs = glob.glob(f'{dirname}/scenes/*')
    num_scenes = len(scene_dirs)

    render_dir = f'{dirname}/renders'
    try:
        os.mkdir(render_dir)
    except:
        pass

    for id_scene, scene in enumerate(scene_dirs):
        msg = f'[{id_scene}/{num_scenes}] {scene}'
        print(msg + ' ' * max(0, 78-len(msg)))

        name = os.path.basename(scene)
        scene_file = f'{scene}/scene.json'

        cmd = f'{bin} render {scene_file} --output {render_dir}\\{name}.png --samples {spp}'
        options = ' --envhidden --envname .\\data\\hdr-images\\doge2.hdr'
        options += ' --denoise'
        options += ' --exposure 1.0'
        # options += ' --embreebvh'

        cmd += options
        print(cmd)

        subprocess.run(cmd, shell=True)


cli()
