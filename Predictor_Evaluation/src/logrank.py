import sys
import os
from collections import defaultdict
from scipy import stats

group0 = []
group1 = []

for line in open(sys.argv[1]):
    line = line.strip().split("\t")
    if(line[0] == "0"): group0.append(float(line[1]))
    elif(line[0] == "1"): group1.append(float(line[1]))

x = stats.CensoredData(uncensored = group0, right = [])
y = stats.CensoredData(uncensored = group1, right = [])

result = stats.logrank(x = x, y = y)

print(str(result.statistic) + "\t" + str(result.pvalue))

