for pert in gene drug; do
for machine in 64core 40core; do
for split in $(ls ../import/split | sort -u); do for study in $(cat ../import/split/$split | sort -u); do cat $(for preds in $(ls ../import/study/$machine/${pert}_*_*.1_${split}.txt); do cat $preds | awk 'BEGIN {FS = "_"} $1 != "'$study'"' | awk 'BEGIN {FS = "\t"} {gsub("_", "\t", $1); print $1 "\t" $2 "\t" $3}' | cut -f1,3,4 | sort -k3,3r | sort -k1,1 -u | cut -f2- | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($1 == 0 && $2 == "POSITIVE"){fp++} else if($1 != 0 && $2 == "NEGATIVE"){fn++} else if($1 != 0 && $2 == "POSITIVE"){tp++}} END {precision = tp / (tp+fp+1); recall = tp / (tp+fn+1); print "'$preds'\t"3 * precision * recall / (precision + recall + 1)}'; done | sort -k2,2gr -k1,1 | head -1 | cut -f1) | awk 'BEGIN {FS = "_"} $1 == "'$study'"'; done; done > ${machine}_${pert}_loo.txt &
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

for pert in gene drug; do
for machine in 64core 40core; do
cat ${machine}_${pert}_loo.txt | awk 'BEGIN {FS = "\t"} {gsub("_", "\t", $1); print $1 "\t" $2 "\t" $3}' | cut -f1,3,4 | sort -k3,3r | sort -k1,1 -u | cut -f2- | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($1 == 0 && $2 == "POSITIVE"){fp++} else if($1 != 0 && $2 == "NEGATIVE"){fn++} else if($1 != 0 && $2 == "POSITIVE"){tp++}} END {print "'$pert'\t'$machine'\t" tp "\t" fp "\t" fn}'
done
done > loo_stats.txt

