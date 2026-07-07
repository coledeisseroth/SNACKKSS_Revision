import sys
import os
from collections import defaultdict

def is_numeric(string):
    try:
        int(string)
        return True
    except ValueError: return False

case_sensitive = False
if "-s" in sys.argv:
    sys.argv.remove("-s")
    case_sensitive = True

for line in open(sys.argv[1]):
    line = line.strip().split("\t")
    sampleid = line[0]
    if not line[1]: continue
    entities = line[1].split(";")
    text = line[2]
    if case_sensitive: text_lower = text
    else: text_lower = text.lower()
    entity_types = {}
    for entity in entities:
        entity = entity.split(":")
        name = ":".join(entity[:-1])
        if not case_sensitive: name = name.lower()
        if not name: continue
        if name in entity_types.keys() and entity_types[name] == "GENE": continue
        if name not in text_lower: continue
        entity_type = entity[-1]
        entity_types[name] = entity_type
    positions = []
    offset = 0
    while text_lower:
        champ_index = None
        champ_name = None
        for name in entity_types.keys():
            if name not in text_lower: continue
            index = text_lower.index(name)
            if champ_index is None or index < champ_index:
                champ_name = name
                champ_index = index
        if champ_name is None: break
        start = champ_index + offset
        end = champ_index + len(champ_name) + offset
        while end < len(text_lower) and is_numeric(text_lower[end - offset]): end += 1
        positions.append(str(start) + "-" + str(end) + ":" + entity_types[champ_name])
        text_lower = text_lower[end:]
        offset += end
    print(sampleid + "\t" + ";".join(positions) + "\t" + text)
    

