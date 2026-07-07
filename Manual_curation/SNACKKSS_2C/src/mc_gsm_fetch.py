import sys
import os
from collections import defaultdict

entries = []
for line in open(sys.argv[1]):
    line = line.strip("\n")
    line = line.split("\t")
    gse = line[1]
    pert_type = line[3]
    if pert_type == "N": continue
    term_type = line[4]
    control_samples = line[5]
    pert_samples = line[6]
    comparison_requirements = line[7].split(";")
    if not comparison_requirements[0]: comparison_requirements = []
    #We're not dealing with ChIP-Seq right now...
    if line[8]:
        print(gse + "\t#NOT_RNA")
        continue
    if "NOT_RNA" in line[11]:
        print(gse + "\t#NOT_RNA")
        continue
    if "NO_SOFT_FILE" in line[11]:
        print(gse + "\t#NO_SOFT_FILE")
        continue
    if "UNUSABLE" in line[11]:
        print(gse + "\t#UNUSABLE")
        continue
    required_features = line[9].split(";")
    if not required_features[0]: required_features = []
    if "&" in line[10] and "\\&" not in line[10]:
        print (gse + "\t#ERROR_MULTIPLE_PERTURBAGENS")
        continue
    perturbagen_terms = line[10].split(";")
    if not perturbagen_terms[0]:
        print(gse + "\t#ERROR_NO_PERTURBAGEN")
        continue
    control_terms = control_samples.split(";")
    if not control_terms[0]:
        print(gse + "\t#NO_CONTROL")
        continue
    pert_sample_terms = pert_samples.split(";")
    if not pert_sample_terms:
        print(gse + "\t#ERROR_NO_PERTURBED_SAMPLES")
        continue
    if term_type == "G":
        print(gse + "\t" + pert_type + "\t" + control_samples + "\t" + pert_samples + "\t" + ";".join(perturbagen_terms))
        continue
    entries.append([gse, pert_type, control_terms, pert_sample_terms, comparison_requirements, required_features, perturbagen_terms, defaultdict(list), defaultdict(list)])

for line in open(sys.argv[2]):
    line = line.strip().split("\t")
    gse = line[0]
    gsm = line[1]
    info = line[2]
    info = info.replace("; ", "\t")
    info_items = info.split("\t")
    for i in range(len(entries)):
        if gse != entries[i][0]: continue
        req_met = True
        for req in entries[i][5]:
            if req not in info:
                req_met = False
                break
        if not req_met: continue
        if entries[i][4]:
            group = []
            for comp_req in entries[i][4]:
                for item in info_items:
                    if comp_req in item: group.append(item)
            group.sort()
            group = "; ".join(group)
        else: group = ""
        control = False
        for cterm in entries[i][2]:
            if cterm in info:
                control = True
                continue
        if control:
            entries[i][7][group].append(gsm)
            continue
        for pterm in entries[i][3]:
            if pterm in info:
                entries[i][8][group].append(gsm)
                break

for entry in entries:
    if not entry[7].keys():
        print(entry[0] + "\t#NO_GROUPS_APPLICABLE")
        continue
    for group in entry[7].keys():
        if (not entry[7][group]) or (not entry[8][group]):
            print(entry[0] + "\t" + group + "\t#GROUP_NOT_FOUND")
            continue
        print(entry[0] + "\t" + entry[1] + "\t" + ";".join(entry[7][group]) + "\t" + ";".join(entry[8][group]) + "\t" + ";".join(entry[6]))



