import json
from sys import argv
from pathlib import Path
from copy import deepcopy

from cjio.cityjson import CityJSON


def merge(cityjson_path, features_paths):
    lcount = 1
    # -- read first line
    with cityjson_path.open("r") as fo:
        j1 = json.load(fo)
    cm = CityJSON(j=j1)
    if "CityObjects" not in cm.j:
        cm.j["CityObjects"] = {}
    if "vertices" not in cm.j:
        cm.j["vertices"] = []
    for path in features_paths:
        if path.suffix == ".jsonl":
            with path.open("r") as fo:
                j1 = json.load(fo)
            if not ("type" in j1 and j1["type"] == 'CityJSONFeature'):
                raise IOError(
                    "Line {} is not of type 'CityJSONFeature'.".format(lcount))
            cm.add_cityjsonfeature(j1)
        else:
            print(f"Not a .jsonl file {path}")
    return cm


if __name__ == "__main__":
    # format to convert to
    supported_formats = ["cityjson", "3dtiles"]
    output_format = argv[1]
    if output_format not in supported_formats:
        raise ValueError(f"Output format {output_format} is not supported. Supported formats: {supported_formats}")
    # where to save the output
    output_file = Path(argv[2])
    # the main .city.json file with the transformation properties
    cityjson_path = Path(argv[3]).resolve()
    # root directory for the .city.jsonl files
    cityjsonfeatures_paths_root = Path(argv[4]).resolve()
    # comma separated list of .city.jsonl files
    cityjsonfeatures_paths = [(cityjsonfeatures_paths_root / p).resolve() for p in argv[5].split(",")]

    cm = merge(cityjson_path, cityjsonfeatures_paths)
    if output_format == "cityjson":
        with output_file.open("w") as fo:
            fo.write(json.dumps(cm.j, separators=(',', ':')))
    elif output_format == "3dtiles":
        cm.reproject(4978)
        # TODO: we are cheating here, because we know that the data has 3 LoD-s and we
        #  also hardcoded the same into Tileset.from() in tyler. Same for the tile/file
        #  names.
        lod_file_names = [("1.2", ""), ("1.3", "-0"), ("2.2", "-0-0")]
        for lod, suffix in lod_file_names:
            cm_copy = deepcopy(cm)
            cm_copy.filter_lod(lod)
            glb = cm_copy.export2glb()
            glb.seek(0)
            output_file_tile = (output_file.parent / (output_file.stem + suffix)).with_suffix(".glb")
            with output_file_tile.open("wb") as bo:
                bo.write(glb.getvalue())
    else:
        raise ValueError("unsupported format and we should have reached this branch anyway")