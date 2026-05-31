#This is just a script to remove confusing characters and lowercase-lock text in a given column of a tab-delimited table
import sys

col = int(sys.argv[2])
for line in open(sys.argv[1]):
    line = line.strip().split("\t")
    line[col] = line[col].lower()
    line[col] = line[col].replace(',', '')
    line[col] = line[col].replace(' ', '')
    line[col] = line[col].replace('-', '')
    print("\t".join(line))

