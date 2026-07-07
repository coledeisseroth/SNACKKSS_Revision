import sys
import os
from collections import defaultdict

def parenthetical_nesting(termsfield):
    returnList = []
    nestings = 0
    curString = ""
    escape = False
    for i in range(len(termsfield)):
        if escape:
            if termsfield[i] in ["\\", "&", ";"]:
                curString += termsfield[i-1:i+1]
            elif termsfield[i] == "(":
                if nestings == 0:
                    if curString: returnList.append(curString)
                    curString = ""
                else: curString += termsfield[i-1:i+1]
                nestings += 1
            elif termsfield[i] == ")":
                nestings -= 1
                if nestings == 0:
                    returnList.append(parenthetical_nesting(curString))
                    curString = ""
                else: curString += termsfield[i-1:i+1]
            escape = False
            continue
        if termsfield[i] == "\\":
            escape = True
            continue
        if nestings > 0:
            curString += termsfield[i]
            continue
        elif termsfield[i] == "&":
            if curString: returnList.append(curString)
            returnList.append("&")
            curString = ""
        elif termsfield[i] == ";":
            if curString: returnList.append(curString)
            returnList.append(";")
            curString = ""
        else: curString += termsfield[i]
    if curString: returnList.append(curString)
    return returnList

def separate_terms(pnest):
    while len(pnest) == 1 and pnest[0].__class__ == list: pnest = pnest[0]
    if "&" in pnest: returnList = ["&"]
    elif ";" in pnest: returnList = [";"]
    else: return pnest
    curList = []
    for i in range(len(pnest)):
        if pnest[i] == returnList[0]:
            returnList.append(curList)
            curList = []
        else: curList.append(pnest[i])
    if curList: returnList.append(curList)
    for i in range(1, len(returnList)):
        returnList[i] = separate_terms(returnList[i])
    return returnList

def qualify_description(termreqs, evalstring):
    if not termreqs: return True
    if len(termreqs) == 1: return termreqs[0] in evalstring
    if termreqs[0] == "&":
        valid = True
        for i in range(1, len(termreqs)):
            if not qualify_description(termreqs[i], evalstring): 
                valid = False
                break
        return valid
    if termreqs[0] == ";":
        valid = False
        for i in range(1, len(termreqs)):
            if qualify_description(termreqs[i], evalstring):
                valid = True
                break
        return valid

#print(qualify_description(separate_terms(parenthetical_nesting(sys.argv[1])), sys.argv[2]))
#quit()

def compreq_split(terms, delimiter):
    terms = terms.replace("\\" + delimiter, "\t")
    terms = terms.split(delimiter)
    for i in range(len(terms)): terms[i] = terms[i].replace("\t", delimiter)
    return terms

#TODO parenthetical order of operations
entries = []
for line in open(sys.argv[1]):
    line = line.strip("\n")
    line = line.split("\t")
    gse = line[0]
    if "GSE" not in gse: continue
    pert_type = line[2]
    if pert_type == "N": continue
    if pert_type == "Not RNA-Seq": continue
    term_type = line[3].strip()
    control_samples = line[4]
    pert_samples = line[5]
    comparison_requirements = compreq_split(line[6], ";")
    if not comparison_requirements[0]: comparison_requirements = []
    required_features = separate_terms(parenthetical_nesting(line[7]))
    perturbagen_terms = line[8]
    if len(compreq_split(perturbagen_terms, "&")) > 1:
        print (gse + "\t#ERROR_MULTIPLE_PERTURBAGENS")
        continue
    perturbagen_terms = compreq_split(perturbagen_terms, ";")
    if not perturbagen_terms[0]:
        print(gse + "\t#ERROR_NO_PERTURBAGEN")
        continue
    control_terms = separate_terms(parenthetical_nesting(control_samples))
    #If the control field is blank in the new dataset, the testing conditions field must be populated.
    if not control_terms and not required_features:
        print(gse + "\t#ERROR_NO_CONTROL")
        continue
    pert_sample_terms = separate_terms(parenthetical_nesting(pert_samples))
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
        if not qualify_description(entries[i][5], info): continue
        if entries[i][4]:
            group = []
            for comp_req in entries[i][4]:
                for item in info_items:
                    if comp_req in item: group.append(item)
            group.sort()
            group = "; ".join(group)
        else: group = ""
        if qualify_description(entries[i][2], info):
            entries[i][7][group].append(gsm)
            continue
        if qualify_description(entries[i][3], info): entries[i][8][group].append(gsm)

for entry in entries:
    if not entry[7].keys():
        print(entry[0] + "\t#NO_GROUPS_APPLICABLE")
        continue
    for group in entry[7].keys():
        if (not entry[7][group]) or (not entry[8][group]):
            print(entry[0] + "\t" + group + "\t#GROUP_NOT_FOUND")
            continue
        print(entry[0] + "\t" + entry[1] + "\t" + ";".join(entry[7][group]) + "\t" + ";".join(entry[8][group]) + "\t" + ";".join(entry[6]))



