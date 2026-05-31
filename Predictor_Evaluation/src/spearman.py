import sys
import os
from collections import defaultdict
from scipy import stats

#x = [1.0, 1.0, -1.0, -1.0]
#y = [1.0, -1.0, 1.0, -1.0]
x = []
y = []

for line in open(sys.argv[1]):
    line = line.strip().split("\t")
    x.append(float(line[0]))
    y.append(float(line[1]))


results = stats.spearmanr(x, y)

print(str(results.correlation) + "\t" + str(results.pvalue))

