for version in old new; do
(cat ../cleaned/${version}_perturbations_cleaned.txt | awk 'BEGIN {FS = "\t"} $2 == "KO" || $2 == "KD"' | cut -f4 | sed 's/;/\n/g' | sort -u | awk '{print $1 "\t1"}'; cat ../cleaned/${version}_perturbations_cleaned.txt | awk 'BEGIN {FS = "\t"} $2 == "KO" || $2 == "KD"' | cut -f4 | sed 's/;/\n/g' | sort -u | comm -13 - <(cat ../cleaned/${version}_perturbations_cleaned.txt | awk 'BEGIN {FS = "\t"} $2 == "KO" || $2 == "KD"' | cut -f1 | sort -u | join -t$'\t' - <(cat ../../import/sample_info.txt | cut -f-2 | sort -k1,1) | cut -f2 | sort -u) | awk '{print $1 "\t0"}') > ${version}_gene_labels.txt
(cat ../cleaned/${version}_perturbations_cleaned.txt | awk 'BEGIN {FS = "\t"} $2 == "D"' | cut -f4 | sed 's/;/\n/g' | sort -u | awk '{print $1 "\t1"}'; cat ../cleaned/${version}_perturbations_cleaned.txt | awk 'BEGIN {FS = "\t"} $2 == "D"' | cut -f4 | sed 's/;/\n/g' | sort -u | comm -13 - <(cat ../cleaned/${version}_perturbations_cleaned.txt | awk 'BEGIN {FS = "\t"} $2 == "D"' | cut -f1 | sort -u | join -t$'\t' - <(cat ../../import/sample_info.txt | cut -f-2 | sort -k1,1) | cut -f2 | sort -u) | awk '{print $1 "\t0"}') > ${version}_drug_labels.txt
done

for pert in gene drug; do
cat ../SNACKKSS_2C.txt | awk '$2 == "Ben"' | cut -f1 | sort -u | join -t$'\t' - <(cat ../../import/sample_info.txt | cut -f-2 | sort -k1,1) | cut -f2 | sort -u | join -t$'\t' - <(cat new_${pert}_labels.txt | sort -k1,1) > independent_${pert}_labels.txt
done

(echo $'Perturbation type\tAgreed\tTrials\tObserved agreement\tExpected agreement\tKappa'
for pert in gene drug; do
cat old_${pert}_labels.txt | sort -u | sort -k1,1 | join -t$'\t' - <(cat independent_${pert}_labels.txt | sort -u | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"; p1 = 0; p2 = 0; n1 = 0; n2 = 0; agreed = 0; disagreed = 0} {if($1 == 1){p1++} else{n1++} if($2 == 1){p2++} else{n2++} if($1 == $2){agreed++} else{disagreed++}} END {rate1 = p1/(p1+n1); rate2=p2/(p2+n2); pe = (rate1 * rate2) + ((1-rate1) * (1-rate2)); po = agreed / (agreed + disagreed); k = (po - pe) / (1 - pe); print "'$pert'\t" agreed "\t" (agreed + disagreed) "\t" po "\t" pe "\t" k}'
done) > agreement_stats.txt

