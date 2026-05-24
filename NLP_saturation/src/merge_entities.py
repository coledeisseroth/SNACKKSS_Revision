import sys
import os
from collections import defaultdict

for line in open(sys.argv[1]):
    line = line.strip().split("\t")
    positions = line[1].split(";")
    if not positions[0]: positions = []
    new_positions = []
    cur_start = None
    cur_end = None
    cur_type = None
    for i in range(len(positions)):
        item = positions[i].split(":")
        entity_type = item[1]
        pair = item[0].split("-")
        start = pair[0]
        end = pair[1]
        if start != cur_end:
            if cur_start is not None: new_positions.append(cur_start + "-" + cur_end + ":" + cur_type)
            cur_start = start
            cur_type = entity_type
        if entity_type == "GENE": cur_type = entity_type
        cur_end = end
    if cur_start is not None: new_positions.append(cur_start + "-" + cur_end + ":" + cur_type)
    print(line[0] + "\t" + ";".join(new_positions) + "\t" + line[2])

