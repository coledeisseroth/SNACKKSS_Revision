for pert in gene drug; do
for machine in 64core 40core; do
for split in $(ls ../import/split/); do for study in $(cat ../import/split/$split | sort -u); do cat $(for preds in $(ls ../import/control/$machine/${pert}_*_*.1_${split}.txt); do 
cat $preds | awk 'BEGIN {FS = "_"} $1 != "'$study'"' | cut -f-2 | awk 'BEGIN {FS = "_"} {print $2 "\t" $0}' | sort -k1,1 | join -t$'\t' - <(cat $preds | awk 'BEGIN {FS = "_"} $1 != "'$study'"' | cut -d_ -f2 | sort | uniq -c | awk '{print $2 "\t" 1 / $1}' | sort -u | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "_" $4 "\t" $7 "\t" $8}' | sort -k1,1 | join -t$'\t' <(cat ../import/labeled_pairs/$pert.txt | cut -f1 | awk 'BEGIN {FS = "_"} {print $1 "_" $2 "_" $3 "\t" $4}' | sort -k1,1) - | sed 's/POSITIVE/1/g' | sed 's/NEGATIVE/0/g' | awk 'BEGIN {FS = "\t"} {if($3 == 1){$3 = 1} else{$3 = 0} print $2 "\t" $3 "\t" $4}' | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($1 == 1 && $2 == 1){tp += $3} else if($1 == 1){fn += $3} else if($2 == 1){fp += $3}} END {print precision = tp / (tp+fp+1); recall = tp / (tp+fn+1); print "'$preds'\t" 3 * precision * recall / (precision + recall + 1)}'; done | sort -k2,2gr -k1,1 | head -1 | cut -f1) | awk 'BEGIN {FS = "_"} $1 == "'$study'"'
done; done > ${machine}_${pert}_loo.txt &
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

for pert in gene drug; do
for machine in 64core 40core; do
cat ${machine}_${pert}_loo.txt | cut -f-2 | awk 'BEGIN {FS = "_"} {print $2 "\t" $0}' | sort -k1,1 | join -t$'\t' - <(cat ${machine}_${pert}_loo.txt | cut -d_ -f2 | sort | uniq -c | awk '{print $2 "\t" 1 / $1}' | sort -u | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "_" $4 "\t" $7 "\t" $8}' | sort -k1,1 | join -t$'\t' <(cat ../import/labeled_pairs/$pert.txt | cut -f1 | awk 'BEGIN {FS = "_"} {print $1 "_" $2 "_" $3 "\t" $4}' | sort -k1,1) - | sed 's/POSITIVE/1/g' | sed 's/NEGATIVE/0/g' | awk 'BEGIN {FS = "\t"} {if($3 == 1){$3 = 1} else{$3 = 0} print $2 "\t" $3 "\t" $4}' | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($1 == 1 && $2 == 1){tp += $3} else if($1 == 1){fn += $3} else if($2 == 1){fp += $3}} END {print "'$pert'\t'$machine'\t" tp "\t" fp "\t" fn}'
done
done > loo_stats.txt

#Unweighted stats
for pert in gene drug; do
for machine in 64core 40core; do
cat ${machine}_${pert}_loo.txt | cut -f-2 | awk 'BEGIN {FS = "_"} {print $2 "\t" $0}' | sort -k1,1 | join -t$'\t' - <(cat ${machine}_${pert}_loo.txt | cut -d_ -f2 | sort | uniq -c | awk '{print $2 "\t" 1 / $1}' | sort -u | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "_" $4 "\t" $7 "\t" $8}' | sort -k1,1 | join -t$'\t' <(cat ../import/labeled_pairs/$pert.txt | cut -f1 | awk 'BEGIN {FS = "_"} {print $1 "_" $2 "_" $3 "\t" $4}' | sort -k1,1) - | sed 's/POSITIVE/1/g' | sed 's/NEGATIVE/0/g' | awk 'BEGIN {FS = "\t"} {if($3 == 1){$3 = 1} else{$3 = 0} print $2 "\t" $3 "\t" $4}' | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($1 == 1 && $2 == 1){tp++} else if($1 == 1){fn++} else if($2 == 1){fp++}} END {print "'$pert'\t'$machine'\t" tp "\t" fp "\t" fn}'
done
done > unweighted_loo_stats.txt


