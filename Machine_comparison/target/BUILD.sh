for pert in gene drug; do
if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then label=GENE
else label=CHEM; fi
for machine in 64core 40core; do
for split in $(ls ../import/split/); do for study in $(cat ../import/split/$split | sort -u); do
optimum=$(for model in distilbert biobert biomedbert; do for preds in $(ls ../import/target/$machine/${pert}_${model}_*.1_$split.txt); do
cat ../import/position_labels/${pert}_${model}/$split.txt | awk 'BEGIN {FS = "_"} {print $1 "\t" $0}' | grep -vwf <(cat ../import/sample_info.txt | cut -f-2 | awk '$1 == "'$study'"' | cut -f2 | sort -u) | cut -f2- | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u | comm - <(cat $preds | awk 'BEGIN {FS = "_"} {print $1 "\t" $0}' | grep -vwf <(cat ../import/sample_info.txt | cut -f-2 | awk '$1 == "'$study'"' | cut -f2 | sort -u) | cut -f2- | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u) | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($2 == ""){tp++} else if($1 == ""){fp++} else{fn++}} END {precision = tp / (tp + fp + 1); recall = tp / (tp + fn + 1); print "'$preds'\t" 3 * precision * recall / (precision + recall + 1)}'
done; done | sort -k2,2gr -k1,1 | head -1 | cut -f1)
cat $optimum | awk 'BEGIN {FS = "_"} {print $1 "\t" $0}' | sort -k1,1 | join -t$'\t' - <(cat ../import/sample_info.txt | cut -f-2 | awk '$1 == "'$study'"' | cut -f2 | sort -u) | cut -f2- | sort -u | awk '{print "'$optimum'\t" $0}' 
done
done > ${machine}_${pert}_loo.txt &
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

for pert in gene drug; do
if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then label=GENE
else label=CHEM; fi
for machine in 64core 40core; do
for model in distilbert biobert biomedbert; do
cat ../import/position_labels/${pert}_${model}/* | sort -k1,1 | join -t$'\t' - <(cat ${machine}_${pert}_loo.txt | awk 'BEGIN {FS = "_"} $2 == "'$model'"' | cut -f2 | sort -u) | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u | comm - <(cat ${machine}_${pert}_loo.txt | awk 'BEGIN {FS = "_"} $2 == "'$model'"' | cut -f2- | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u); done | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($2 == ""){tp++} else if($1 == ""){fp++} else{fn++}} END {precision = tp / (tp + fp + 1); recall = tp / (tp + fn + 1); print "'$pert'\t'$machine'\t" tp "\t" fp "\t" fn}'
done
done > loo_stats.txt

