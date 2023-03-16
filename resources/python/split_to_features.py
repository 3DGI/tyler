"""Split a directory of cityjson files to cityjsonfeature files, one file per feature.
One main cityjson file (metadata.city.json) is written for all the features and the
'transform' property is computed from the extent of all files.

Copyright 2023 Balázs Dukai, Ravi Peters

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
from pathlib import Path
from sys import argv
from os import cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed
import json

from cjio import cityjson

# zwaartepunt in Nederland
TRANSLATE = [171800.0, 472700.0, 0.0]
IMPORTANT_DIGITS = 3

input_dir = Path(argv[1]).resolve()
output_dir = Path(argv[2]).resolve()
max_workers = int(argv[3]) if int(argv[3]) <= cpu_count() else cpu_count()

if input_dir.is_dir():
    cityjson_path_list = [cityjson_path for cityjson_path in input_dir.iterdir()]
else:
    cityjson_path_list = [input_dir, ]
if len(cityjson_path_list) <= max_workers:
    max_workers = len(cityjson_path_list)
print(f"Using maximum {max_workers} workers")


# --- Compute translation properties from the extent

# Compute the extent of all the files
def bbox_from_file(path: Path):
    with path.open(mode='r', encoding='utf-8-sig') as f:
        cm = cityjson.reader(file=f, ignore_duplicate_keys=False)
        return cm.get_bbox()


futures = []

# Init the extent
with cityjson_path_list[0].open(mode='r', encoding='utf-8-sig') as f:
    cm = cityjson.reader(file=f, ignore_duplicate_keys=False)
    extent = cm.calculate_bbox()

with ProcessPoolExecutor(max_workers=max_workers) as executor:
    for cj_path in cityjson_path_list:
        futures.append(executor.submit(bbox_from_file, cj_path))
    for i, future in enumerate(as_completed(futures)):
        minx, miny, minz, maxx, maxy, maxz = future.result()
        if minx < extent[0]:
            extent[0] = minx
        if miny < extent[1]:
            extent[1] = miny
        if minz < extent[2]:
            extent[2] = minz
        if maxx > extent[3]:
            extent[3] = maxx
        if maxy > extent[4]:
            extent[4] = maxx
        if maxz > extent[5]:
            extent[5] = maxz
del cm, i, future, futures, executor

# The features are centered around the center of the extent
dx = extent[3] - extent[0]
dy = extent[4] - extent[1]
dz = extent[5] - extent[2]
center = [extent[0] + dx * 0.5, extent[1] + dy * 0.5, extent[2] + dz * 0.5]
translate = TRANSLATE
print(f"Computed translation property: {translate}")


# --- Write the metadata file
with cityjson_path_list[0].open(mode='r', encoding='utf-8-sig') as f:
    cm = cityjson.reader(file=f, ignore_duplicate_keys=False)
    cm.upgrade_version("1.1", digit=IMPORTANT_DIGITS)
    cm.decompress()
    cm.compress(important_digits=IMPORTANT_DIGITS, translate=translate)
    outfile = output_dir / "metadata.city.json"
    with outfile.open("w") as fo:
        fo.write(cm.cityjson_for_features())
    print(f"Written {outfile}")
del cm, outfile

# --- Split to features
def file_to_feature_files(filepath: Path, out_dir: Path):
    with filepath.open(mode='r', encoding='utf-8-sig') as f:
        cm = cityjson.reader(file=f, ignore_duplicate_keys=False)
        cm.upgrade_version("1.1", digit=IMPORTANT_DIGITS)
        cm.decompress()
        cm.compress(important_digits=IMPORTANT_DIGITS, translate=translate)

        fail = []
        # e.g: 'gb2' in /home/cjio/gb2.city.json
        old_filename = filepath.name.replace("".join(filepath.suffixes), "")
        # e.g: '/home/cjio/gb2' in /home/cjio/gb2.city.json
        filedir = out_dir / old_filename
        filedir.mkdir(exist_ok=True)
        for feature in cm.generate_features():
            feature_id = feature.j['id']
            new_filename = f"{feature_id}.city.jsonl"
            filepath = filedir / new_filename
            try:
                with open(filepath, "w") as fout:
                    fout.write(json.dumps(feature.j, separators=(',', ':')))
            except IOError as e:
                print(f"Invalid output file: {filepath}\n{e}")
                fail.append(feature_id)
            except BaseException as e:
                print(e)
                fail.append(feature_id)
        return fail


failed = []
futures = []
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    for cj_path in cityjson_path_list:
        futures.append(executor.submit(file_to_feature_files, cj_path, output_dir))
    for i, future in enumerate(as_completed(futures)):
        failed.extend(future.result())
del i, future, futures, executor

print(f"Failed to export the CityObjects: {failed}")
