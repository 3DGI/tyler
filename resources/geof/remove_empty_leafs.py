import json
import os, sys
from pathlib import Path

TILESET_DIR = Path("/dev/shm/3dtiles-buildings")

# filename = TILESET_DIR / "tileset_og.json"
# output = TILESET_DIR / "tileset.json"

filename = sys.argv[1]
output = sys.argv[2]

TILESET_DIR = Path(filename).parent
print(TILESET_DIR)

def remove_empty_leafs(node):
  if node.get("children"):
    # print("node has {} children".format(len(node["children"])))
    # print(node["children"])
    new_children = []
    for child in node["children"]:
      if child.get("content"):
        # print("checking " + str(TILESET_DIR / child["content"]["uri"]))
        if(os.path.exists(TILESET_DIR / child["content"]["uri"])):
          new_children.append(child)
          # print("exists")
      else:
        new_children.append(child)
    if len(new_children) == 0:
      del node["children"]
    else:
      node["children"] = new_children
    # print(node["children"])
      for child in node["children"]:
        remove_empty_leafs(child)

with open(filename, "r") as in_file:
    tileset = json.load(in_file)

remove_empty_leafs(tileset["root"])

with open(output, "w") as file:
    json.dump(tileset, file)
