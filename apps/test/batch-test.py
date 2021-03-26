#! /usr/bin/env python3 -B

import sys
import glob
import os
import subprocess
import json

def trace_meshes(dirname, output):
    mesh_names = glob.glob(f'{dirname}/meshes/*')
    mesh_num = len(mesh_names)

    images_dir = f'{output}/images'
    try:
        os.mkdir(images_dir)
    except:
        pass

    result = {}
    result['num_tests'] = 0
    result['num_errors'] = 0
    result['ok'] = []

    append = ''
    for mesh_id, mesh_name in enumerate(mesh_names):
        result['num_tests'] += 1
        name = os.path.basename(mesh_name).split('.')[0]
        # stats_name = f'{output}/stats/{name}.json'
        # scene_name = f'{output}/scenes/{name}/scene.json'
        # curve_name = name.replace('meshes/', 'curves/') + '.json'
        # scene_name = name.replace('meshes/', 'scenes/') + '.json'
        msg = f'[{mesh_id}/{mesh_num}] {mesh_name}'
        print(msg + ' ' * max(0, 78-len(msg)))

        cmd = f'./bin/Debug/test --model {mesh_name} --output {images_dir}/{name}.png data/svgs/abc.json'
        print(cmd)
        if append == '': append = '--append-timings'
        
        try:
            retcode = subprocess.run(cmd, timeout=60, shell=True).returncode
            if retcode < 0:
                result['num_errors'] += 1
            elif retcode > 0:
                result['num_errors'] += 1
            else:
                result['ok'] += [mesh_name]
        except OSError:
            result['num_errors'] += 1
        except subprocess.TimeoutExpired:
            result['num_errors'] += 1
    
    # with open(f'{output}/trace-result.json', 'wt') as f:
        # json.dump(result, f, indent=2)



def trace_jsons(dirname, output):
    jsons_names = glob.glob(f'{dirname}/tests/*')
    jsons_num = len(jsons_names)

    images_dir = f'{output}/images'
    try:
        os.mkdir(images_dir)
    except:
        pass

    result = {}
    result['num_tests'] = 0
    result['num_errors'] = 0
    result['ok'] = []

    append = ''
    for json_id, json_name in enumerate(jsons_names):
        result['num_tests'] += 1
        name = os.path.basename(json_name).split('.')[0]
        # stats_name = f'{output}/stats/{name}.json'
        # scene_name = f'{output}/scenes/{name}/scene.json'
        # curve_name = name.replace('jsones/', 'curves/') + '.json'
        # scene_name = name.replace('jsones/', 'scenes/') + '.json'
        msg = f'[{json_id}/{jsons_num}] {json_name}'
        print(msg + ' ' * max(0, 78-len(msg)))

        cmd = f'./bin/test {json_name} --output {images_dir}/{name}.png'
        print(cmd)
        if append == '': append = '--append-timings'
        
        try:
            retcode = subprocess.run(cmd, timeout=60, shell=True).returncode
            if retcode < 0:
                result['num_errors'] += 1
            elif retcode > 0:
                result['num_errors'] += 1
            else:
                result['ok'] += [json_name]
        except OSError:
            result['num_errors'] += 1
        except subprocess.TimeoutExpired:
            result['num_errors'] += 1



if __name__ == '__main__':
    dir = sys.argv[1]
    task = ''    
    if len(sys.argv) >= 3:
        task = sys.argv[2]

    output = f'{dir}/output'
    try:
        os.mkdir(output)
    except:
        pass

    if(task == 'jsons'):
        trace_jsons(dir, output)
    else:
        trace_meshes(dir, output)



