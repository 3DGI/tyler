from pathlib import Path
import json

tileset_path = Path("/data/work/bdukai/3dbag_v2103_features/3dtiles/tileset_explicit.json")
with tileset_path.open("r") as fo:
    tileset = json.load(fo)


def printTree(root, markerStr="+- ", levelMarkers=[]):
    """https://simonhessner.de/python-3-recursively-print-structured-tree-including-hierarchy-markers-using-depth-first-search/"""
    emptyStr = " " * len(markerStr)
    connectionStr = "|" + emptyStr[:-1]
    level = len(levelMarkers)
    mapper = lambda draw: connectionStr if draw else emptyStr
    markers = "".join(map(mapper, levelMarkers[:-1]))
    markers += markerStr if level > 0 else ""
    if "content" in root:
        print(f"{markers}{level}: {root['geometricError']} {root['content']['uri']}")
    else:
        print(f"{markers}{level}: {root['geometricError']} {root['content']['uri']}")
    if "children" in root:
        for i, child in enumerate(root["children"]):
            isLast = i == len(root["children"]) - 1
            printTree(child, markerStr, [*levelMarkers, not isLast])


if __name__ == "__main__":
    printTree(tileset["root"])
