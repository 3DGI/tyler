"""Copyright 2023 Balázs Dukai, Ravi Peters"""
import json
import argparse
from sys import argv
from pathlib import Path
from copy import deepcopy

from cjio.cityjson import CityJSON


def merge(cityjson_path, path_features_input_file: Path):
    lcount = 1
    # -- read first line
    with cityjson_path.open("r") as fo:
        j1 = json.load(fo)
    cm = CityJSON(j=j1)
    if "CityObjects" not in cm.j:
        cm.j["CityObjects"] = {}
    if "vertices" not in cm.j:
        cm.j["vertices"] = []
    with path_features_input_file.open("r") as input_file:
        for p in input_file:
            path = Path(p.strip("\n")).resolve()
            if path.suffix == ".jsonl":
                with path.open("r") as fo:
                    j1 = json.load(fo)
                if not ("type" in j1 and j1["type"] == 'CityJSONFeature'):
                    raise IOError(
                        "Line {} is not of type 'CityJSONFeature'.".format(lcount))
                cm.add_cityjsonfeature(j1)
            else:
                print(f"Not a .jsonl file {path}, suffix: {path.suffix}")
    path_features_input_file.unlink()
    return cm


parser=argparse.ArgumentParser()
parser.add_argument("--output_format", help="Format to convert to")
parser.add_argument("--output_file", help="Where to save the output")
parser.add_argument("--path_metadata", help="The main .city.json file with the transformation properties")
parser.add_argument("--path_features_input_file", help="File with the list of feature paths")
parser.add_argument("--min_x", help="Bounding box minimum x coordinate")
parser.add_argument("--min_y", help="Bounding box minimum y coordinate")
parser.add_argument("--min_z", help="Bounding box minimum z coordinate")
parser.add_argument("--max_x", help="Bounding box maximum x coordinate")
parser.add_argument("--max_y", help="Bounding box maximum y coordinate")
parser.add_argument("--max_z", help="Bounding box maximum z coordinate")
parser.add_argument("--cotypes", help="Comma separated list of CityObject types to include in the tile")
parser.add_argument("--metadata_class", help="The name of the metadata class to create (for EXT_structural_metadata).")
parser.add_argument("--attribute_spec", help="The CityObject attribute to include in the output.")
parser.add_argument("--geometric_error")

if __name__ == "__main__":
    args = parser.parse_args()
    # format to convert to
    supported_formats = ["cityjson", "3dtiles"]
    output_format = args.output_format
    if output_format not in supported_formats:
        raise ValueError(f"Output format {output_format} is not supported. Supported formats: {supported_formats}")
    # where to save the output
    output_file = Path(args.output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    # the main .city.json file with the transformation properties
    cityjson_path = Path(args.path_metadata).resolve()
    # a file with the list of feature paths
    path_features_input_file = Path(args.path_features_input_file).resolve()
    # tile bbox (used in the terrain splitter)
    bbox = [args.min_x, args.min_y, args.min_z, args.max_x, args.max_y, args.max_z]
    # CityObject types to use
    cotypes = args.cotypes.split(",")

    cm = merge(cityjson_path, path_features_input_file)
    nr_co_after_merge = len(cm.j["CityObjects"])
    types_after_merge = cm.get_info()
    if output_format == "cityjson":
        with output_file.open("w") as fo:
            fo.write(json.dumps(cm.j, separators=(',', ':')))
    elif output_format == "3dtiles":
        cm.reproject(4978)
        cm = cm.get_subset_cotype(cotypes)
        nr_co_after_subset = len(cm.j["CityObjects"])
        if nr_co_after_subset == 0:
            if nr_co_after_merge > 0:
                print(f"CityModel {output_file} does not contain any {cotypes} objects, but it did contain {nr_co_after_merge} objects of type {types_after_merge}.")
            else:
                print(
                    f"CityModel {output_file} does not contain any {cotypes} objects and it did not contain any objects at all after merge either.")
        glb = cm.export2glb(do_triangulate=False)
        glb.seek(0)
        with output_file.open("wb") as bo:
            bo.write(glb.getvalue())

        # # TODO: we are cheating here, because we know that the data has 3 LoD-s and we
        # #  also hardcoded the same into Tileset.from() in tyler. Same for the tile/file
        # #  names.
        # # lod_file_names = [("1.2", ""), ("1.3", "-0"), ("2.2", "-0-0")] # grid with multi-lod
        # # lod_file_names = [("2.2", "-0-0"),] # for grid with single lod
        # lod_file_names = [("2.2", ""),] # for quadtree
        # for lod, suffix in lod_file_names:
        #     cm_copy = deepcopy(cm)
        #     cm_copy.filter_lod(lod)
        #     glb = cm_copy.export2glb(do_triangulate=False)
        #     glb.seek(0)
        #     output_file_tile = (output_file.parent / (output_file.stem + suffix)).with_suffix(".glb")
        #     with output_file_tile.open("wb") as bo:
        #         bo.write(glb.getvalue())
    else:
        raise ValueError("unsupported format and we should have reached this branch anyway")

    # /home/bdukai/software/cjio/venv_38/lib/python3.8/site-packages/pyproj/transformer.py:197: UserWarning: Best transformation is not available due to missing Grid(short_name=nl_nsgi_nlgeo2018.tif, full_name=, package_name=, url=https://cdn.proj.org/nl_nsgi_nlgeo2018.tif, direct_download=True, open_license=True, available=False)