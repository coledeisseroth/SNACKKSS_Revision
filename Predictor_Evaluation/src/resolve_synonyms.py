#The table must have three columns: (1) The non-mapped data; (2) the name of the entity you're mapping; (3) the ID candidate this name has been mapped to. It must be sorted in order of column 1. 
import sys
from collections import defaultdict

id_to_names = defaultdict(list)
name_to_ids = defaultdict(list)
names = []
cur_subject = None

def longest_synonym_id(name, name_to_ids, id_to_names):
    champ_length = len(name)
    candidates = []
    for entity_id in name_to_ids[name]:
        for candidate_name in id_to_names[entity_id]:
            if len(candidate_name) < champ_length: continue
            if len(candidate_name) > champ_length:
                candidates = []
                champ_length = len(candidate_name)
            candidates.append(entity_id)
    candidates.sort()
    return candidates[0]

for line in open(sys.argv[1]):
    line = line.strip().split("\t")
    subject = line[0]
    name = line[1]
    if cur_subject is not None and subject != cur_subject:
        for curname in names:
            champ_id = longest_synonym_id(curname, name_to_ids, id_to_names)
            print(cur_subject + "\t" + curname + "\t" + str(champ_id))
        id_to_names = defaultdict(list)
        name_to_ids = defaultdict(list)
        names = []
    cur_subject = subject
    entity_id = int(line[2])
    names.append(name)
    name_to_ids[name].append(entity_id)
    id_to_names[entity_id].append(name)

if cur_subject is not None:
    for name in names:
        champ_id = longest_synonym_id(name, name_to_ids, id_to_names)
        print(cur_subject + "\t" + name + "\t" + str(champ_id))

