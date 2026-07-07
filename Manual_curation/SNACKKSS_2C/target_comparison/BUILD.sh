for version in new old; do
cat ../cleaned/${version}_perturbations_cleaned.txt | awk '$2 == "KO" || $2 == "KD"' | cut -f4- | awk 'BEGIN {FS = "\t"} {gsub(";", "\t", $1); print $2 "\t" $1}' | awk 'BEGIN {FS = "\t"} {for(i=2; i <= NF; i++){print $i "\t" $1}}' | sort -u | sed 's/;/\t/g' | sed 's/&/\t/g' | sort -k1,1 | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++){print $1 "\t" $i}}' | sort -u | sort -k1,1 | awk 'BEGIN {FS = "\t"; cur = ""; terms = ""} {if(cur != "" && cur != $1){print cur "\t" terms; terms = ""} cur = $1; terms = terms ";" $2} END {print cur "\t" terms}' | sed 's/\t;*/\t/g' | sed 's/;/:GENE;/g' | awk '{print $0 ":GENE"}' > ${version}_gene_target_names.txt
cat ../cleaned/${version}_perturbations_cleaned.txt | awk '$2 == "D"' | cut -f4- | awk 'BEGIN {FS = "\t"} {gsub(";", "\t", $1); print $2 "\t" $1}' | awk 'BEGIN {FS = "\t"} {for(i=2; i <= NF; i++){print $i "\t" $1}}' | sort -u | sed 's/;/\t/g' | sed 's/&/\t/g' | sort -k1,1 | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++){print $1 "\t" $i}}' | sort -u | sort -k1,1 | awk 'BEGIN {FS = "\t"; cur = ""; terms = ""} {if(cur != "" && cur != $1){print cur "\t" terms; terms = ""} cur = $1; terms = terms ";" $2} END {print cur "\t" terms}' | sed 's/\t;*/\t/g' | sed 's/;/:CHEM;/g' | awk '{print $0 ":CHEM"}' > ${version}_drug_target_names.txt
done

for version in new old; do
for pert in gene drug; do
cat ../../import/sample_info_protocols.txt | cut -f2- | sort -k1,1 | join -t$'\t' - <(cat ${version}_${pert}_target_names.txt | sort -k1,1) | awk 'BEGIN {FS = "\t"} {print $1 "\t" $3 "\t" $2}' | sort -u > ${version}_${pert}_name_labels.txt
python3 ../src/find_names.py <(cat ${version}_${pert}_name_labels.txt | sed 's/\r//g') | awk 'BEGIN {FS = "\t"} $2 != ""' > ${version}_${pert}_position_labels.txt
done
done

for pert in gene drug; do
cat ../SNACKKSS_2C.txt | awk '$2 == "Ben"' | cut -f1 | sort -u | join -t$'\t' - <(cat ../../import/sample_info.txt | cut -f-2 | sort -k1,1) | cut -f2 | sort -u | join -t$'\t' - <(cat new_${pert}_position_labels.txt | sort -k1,1) > independent_${pert}_position_labels.txt
done

#Get the disagreements
(for pert in gene drug; do
cat old_${pert}_position_labels.txt | sort -k1,1 | join -t$'\t' - <(cat independent_${pert}_position_labels.txt | cut -f1 | sort -u) | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i=2; i <= NF; i++){print $1 "\t" $i}}' | sort -u | comm - <(cat independent_${pert}_position_labels.txt | sort -k1,1 | join -t$'\t' - <(cat old_${pert}_position_labels.txt | cut -f1 | sort -u) | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i=2; i <= NF; i++){print $1 "\t" $i}}' | sort -u) | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($2 == ""){tp++} else if($1 == ""){fp++} else{fn++}} END {print "'$pert'\t" tp "\t" fp "\t" fn}'
done | awk 'BEGIN {print "Perturbation type\tShared\tOld only\tNew only"} {print}') > agreement_stats.txt

