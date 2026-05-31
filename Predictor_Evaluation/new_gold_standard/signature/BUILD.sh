#Do the required name/ID mapping
cat ../corpora/trrust/noreactome_mod_to_entrez.txt | sort -k2,2 | join -t$'\t' -1 1 -2 2 <(cat ../../import/ConnectivityMap/shrna_entrez_1to1.txt ../../import/ConnectivityMap/xpr_entrez_1to1.txt ../../import/ConnectivityMap/oe_entrez_1to1.txt | sort -k1,1) - | cut -f2- | sort -k2,2 | join -t$'\t' -1 1 -2 2 <(cat ../../import/ConnectivityMap/shrna_entrez_1to1.txt ../../import/ConnectivityMap/xpr_entrez_1to1.txt ../../import/ConnectivityMap/oe_entrez_1to1.txt | sort -k1,1) - | cut -f2- | sort -u > noreactome_names.txt
cat ../corpora/trrust/noreactome_mod_to_entrez.txt  | sort -k1,1 | join -t$'\t' <(cat ../../import/ConnectivityMap/shrna_entrez_1to1.txt ../../import/ConnectivityMap/xpr_entrez_1to1.txt ../../import/ConnectivityMap/oe_entrez_1to1.txt | sort -k1,1) - | cut -f2- | sort -u > noreactome_modnames.txt
cat ../corpora/nodgidb_aggregated_mod_to_entrez.txt  | sort -k2,2 | join -t$'\t' -1 1 -2 2 <(cat ../../import/ConnectivityMap/shrna_entrez_1to1.txt ../../import/ConnectivityMap/xpr_entrez_1to1.txt ../../import/ConnectivityMap/oe_entrez_1to1.txt | sort -k1,1) - | cut -f2- | sort -u | awk 'BEGIN {FS = "\t"} {print $2 "\t" $1 "\t" $3}' | sort -u > nodgidb_targnames.txt

#Match gene-disruption signatures
for db in SNACKKSS ConnectivityMap; do
for stage in f1_matches predictions; do
outfolder=f1_matches
if [ $(echo $stage | grep predictions | wc -l) -gt 0 ]; then outfolder=pre_predictions; fi
mkdir $db/$outfolder
for sign in pos neg; do
mkdir $db/$outfolder/$sign
thresh=$(cat ../../import/signature_based_optima.txt | awk '$1 == "'$db'" && $2 == "default" && $3 == "gene" && $4 == "'$stage'" && $5 == "'$sign'" {print $6}')
nestedz=../../import/${db}/gene_nested_z
for targ in $(ls $nestedz | cut -d. -f1 | sort -u); do
while [ $(jobs | grep Running | wc -l) -gt 30 ]; do jobs; sleep 1; done
groundtruth=../corpora/trrust/noreactome_mod_to_entrez.txt
if [ $(echo $db | grep ConnectivityMap | wc -l) -gt 0 ]; then groundtruth=noreactome_names.txt; fi
for mod in $(ls $nestedz | cut -d. -f1 | sort -u | comm -12 - <(cat $groundtruth | cut -f1 | sort -u)); do
cat $nestedz/$mod.txt | awk 'sqrt($2^2) > '$thresh | sort -k1,1 | join -t$'\t' -a 1 -a 2 -o auto -e 0 - <(cat $nestedz/$targ.txt | awk 'sqrt($2^2) > '$thresh | sort -k1,1) | awk 'BEGIN {FS = "\t"; ptp = 0; pfp = 0; pfn = 0; ntp = 0; nfp = 0; nfn = 0} {if(sqrt($3^2) <= '$thresh' && sqrt($2^2) <= '$thresh'){next} else if(sqrt($3^2) <= '$thresh'){pfp++; nfp++} else if(sqrt($2^2) <= '$thresh'){pfn++; nfn++} else if($2 * $3 > 0){ptp++; nfp++} else if($2 * $3 < 0){ntp++; pfp++}} END {pprec = ptp / (ptp+pfp+1); prec = ptp / (ptp+pfn+1); nprec = ntp / (ntp+nfp+1); nrec = ntp / (ntp+nfn+1); pf1 = 3 * pprec * prec / (pprec + prec+1); nf1 = 3 * nprec * nrec / (nprec + nrec+1); print "'$mod'\t'$targ'\t" ptp "\t" pfp "\t" pfn "\t" ntp "\t" nfp "\t" nfn "\t" pf1 "\t" nf1 "\t" (pf1-nf1)*2*sqrt((pf1-nf1)^2)/(pf1+nf1+1)}' | cut -f1,2,11 | awk '$3 != 0'
done > $db/$outfolder/$sign/$targ.txt &
done
done
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#Make fast reference files for easy ARCHS4-linking
ls SNACKKSS/pre_predictions/* | grep .txt | cut -d. -f1 | sort -u | join -t$'\t' -1 1 -2 2 - <(cat ../../import/archs4/human_column_gene_entrez.txt | cut -f2- | sort -t$'\t' -k2,2) | cut -f2 | sort -u | awk '{print "#GENE#" $1 "#GENE#"}' | sort -u | grep -wf - <(cat ../../import/archs4/human_correlation_table.txt | awk '{print "#GENE#" $1 "#GENE#\t" $0}') | cut -f2- | sort -u | sort -k1,1 > SNACKKSS/shortened_human_correlation_table.txt
ls ConnectivityMap/pre_predictions/* | rev | cut -d. -f2- | rev | sort -u | awk '{print "#GENE#" $0 "#GENE#"}' | sort -u | grep -wf - <(cat ../../import/archs4/human_correlation_table.txt | awk 'BEGIN {FS = "\t"} {print "#GENE#" $1 "#GENE#\t" $0}') | cut -f2- | sort -u | sort -k1,1 > ConnectivityMap/shortened_human_correlation_table.txt

#Gather the needed correlations from ARCHs4
for db in SNACKKSS ConnectivityMap; do
mkdir $db/correlations
for pair in $(cat ../../import/archs4/human_column_gene_entrez.txt | cut -f1,3 | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | cut -f2 | sort -u) | awk '{print $2 "_" $1}' | sort -u); do
column=$(echo $pair | cut -d_ -f1)
entrez=$(echo $pair | cut -d_ -f2)
while [ $(jobs | grep Running | wc -l) -gt 30 ]; do jobs; sleep 1; done
cat $db/shortened_human_correlation_table.txt | cut -f1,$column > $db/correlations/$entrez.txt &
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#SNACKKSS relations
for sign in pos neg; do
mkdir -p SNACKKSS/prepreds_quickref/$sign/
for mod in $(ls ../../import/SNACKKSS/gene_nested_z/ | cut -d. -f1 | sort -u); do
jobs
while [ $(jobs | grep Running | wc -l) -gt 30 ]; do jobs; sleep 1; done
cat SNACKKSS/pre_predictions/$sign/* | awk '$1 == "'$mod'"' | cut -f2- | sort -k1,1 | join -t$'\t' -1 2 -2 1 <(cat ../../import/archs4/human_column_gene_entrez.txt | cut -f2- | sort -t$'\t' -k2,2) - | cut -f2- | sort -u | sort -k1,1 > SNACKKSS/prepreds_quickref/$sign/$mod.txt &
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#CMap relations
for sign in pos neg; do
mkdir -p ConnectivityMap/prepreds_quickref/$sign/
for mod in $(ls ../../import/ConnectivityMap/gene_nested_z/ | cut -d. -f1 | sort -u); do
jobs
while [ $(jobs | grep Running | wc -l) -gt 30 ]; do jobs; sleep 0.1; done
cat ConnectivityMap/pre_predictions/$sign/* | awk '$1 == "'$mod'"' | cut -f2- | sort -u | sort -k1,1 | join -t$'\t' - <(cat ../../import/archs4/human_column_gene_entrez.txt | cut -f2 | sort -u) > ConnectivityMap/prepreds_quickref/$sign/$mod.txt &
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done
rm $(wc -l ConnectivityMap/prepreds_quickref/*/* | awk '$1 == 0 {print $2}')

#Calculate SA4 gene-gene relationship predictions
for db in SNACKKSS ConnectivityMap; do
groundtruth=../corpora/trrust/noreactome_mod_to_entrez.txt
if [ $(echo $db | grep ConnectivityMap | wc -l) -gt 0 ]; then groundtruth=noreactome_modnames.txt; fi
mkdir $db/predictions
for sign in pos neg; do
mkdir $db/predictions/$sign
for targ in $(ls $db/correlations | cut -d. -f1 | sort -u); do
jobs
while [ $(jobs | grep Running | wc -l) -gt 30 ]; do jobs; sleep 0.1; done
for mod in $(ls $db/prepreds_quickref/$sign/ | rev | cut -d. -f2- | rev | sort -u | comm -12 - <(cat $groundtruth | awk '$2 == "'$targ'"' | cut -f1 | sort -u)); do
cat $db/correlations/$targ.txt | join -t$'\t' - <(cat $db/prepreds_quickref/$sign/$mod.txt) | cut -f2- | awk 'BEGIN {FS = "\t"; p = 0; n = 0} $1 != 0 && $2 != 0 {L = $1 * $2 / (sqrt($1^2) + sqrt($2^2)); if(L > 0){p += L} else{n -= L}} END {if(p - n != 0){print "'$mod'\t'$targ'\t"(p - n) * sqrt((p-n)^2) / (p+n)}}'
done > $db/predictions/$sign/$targ.txt &
done
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#ID mapping
mkdir ConnectivityMap/id_f1_matches ConnectivityMap/id_predictions
for sign in pos neg; do
for targ in $(ls ConnectivityMap/f1_matches/$sign | cut -d. -f1 | sort -u | comm -12 - <(cat noreactome_names.txt | cut -f2 | sort -u)); do
jobs
cat ConnectivityMap/f1_matches/$sign/$targ.txt | sort -k2,2 | join -t$'\t' -1 2 -2 2 <(cat ../../import/ConnectivityMap/shrna_entrez_1to1.txt ../../import/ConnectivityMap/xpr_entrez_1to1.txt | sort -u | sort -k2,2) - | cut -f2- | sort -k2,2 | join -t$'\t' -1 2 -2 2 <(cat ../../import/ConnectivityMap/shrna_entrez_1to1.txt ../../import/ConnectivityMap/xpr_entrez_1to1.txt | sort -u | sort -k2,2) - | cut -f2- | awk '$3 != 0 {print $0 "\t" $3 / sqrt($3^2)}' | sort -t$'\t' -k4,4gr | sort -k1,1 -k2,2 -u | cut -f-3 | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | cut -f-2 | sed 's/\t/_/g' | sort -u) | sed 's/_/\t/g' | sort -u
done > ConnectivityMap/id_f1_matches/gene_$sign.txt &
for targ in $(ls ConnectivityMap/predictions/$sign | cut -d. -f1 | sort -u | comm -12 - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | cut -f2 | sort -u)); do
cat ConnectivityMap/predictions/$sign/$targ.txt | sort -k1,1 | join -t$'\t' -1 2 -2 1 <(cat ../../import/ConnectivityMap/shrna_entrez_1to1.txt ../../import/ConnectivityMap/xpr_entrez_1to1.txt | sort -u | sort -k2,2) - | cut -f2- | awk '$3 != 0 {print $0 "\t" $3 / sqrt($3^2)}' | sort -t$'\t' -k4,4gr | sort -k1,1 -k2,2 -u | cut -f-3 | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | cut -f-2 | sed 's/\t/_/g' | sort -u) | sed 's/_/\t/g' | sort -u
done > ConnectivityMap/id_predictions/gene_$sign.txt &
for targ in $(ls ../../import/ConnectivityMap/drug_f1_matches/$sign/ | rev | cut -d. -f2- | rev | sort -u | comm -12 - <(cat nodgidb_targnames.txt | cut -f2 | sort -u)); do
cat ../../import/ConnectivityMap/drug_f1_matches/$sign/$targ.txt | sort -k2,2 | join -t$'\t' -1 2 -2 2 <(cat ../../import/ConnectivityMap/shrna_entrez_1to1.txt ../../import/ConnectivityMap/xpr_entrez_1to1.txt | sort -u | sort -k2,2) - | cut -f2- | awk '{print $2 "_" $1 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat ../corpora/nodgidb_aggregated_mod_to_entrez.txt | cut -f-2 | sed 's/\t/_/g' | sort -u) | sed 's/_/\t/g' | sort -u
done > ConnectivityMap/id_f1_matches/drug_${sign}.txt &
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#The positive- and negative-optimized parameters must agree on the direction.
for pert in gene drug; do
for stage in f1_matches predictions; do
if [ $(echo $stage $pert | grep drug | grep predictions | wc -l) -gt 0 ]; then continue; fi
cat ConnectivityMap/id_${stage}/${pert}_pos.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' -a 1 -a 2 -o auto -e 0 - <(cat ConnectivityMap/id_${stage}/${pert}_neg.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {if($3 * $4 <= 0){next} if($3 < 0){print $1 "\t" $2 "\t" $3} else{print $1 "\t" $2 "\t" $4}}' | sort -u > ConnectivityMap/id_${stage}/${pert}_both.txt &
done
done
for targ in $(ls SNACKKSS/f1_matches/* | cut -d. -f1 | sort -u | comm -12 - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | cut -f2 | sort -u)); do
cat SNACKKSS/f1_matches/pos/$targ.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat SNACKKSS/f1_matches/neg/$targ.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {if($3 * $4 <= 0){next} if($3 < 0){print $1 "\t" $2 "\t" $3} else{print $1 "\t" $2 "\t" $4}}' | sort -u
done > SNACKKSS/f1_matches/gene_both.txt &
cat SNACKKSS/predictions/pos/* | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat SNACKKSS/predictions/neg/* | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {if($3 * $4 <= 0){next} if($3 < 0){print $1 "\t" $2 "\t" $3} else{print $1 "\t" $2 "\t" $4}}' | sort -u > SNACKKSS/predictions/gene_both.txt &
for targ in $(ls ../../import/SNACKKSS/drug_f1_matches/* | cut -d. -f1 | sort -u | comm -12 - <(cat ../corpora/nodgidb_aggregated_mod_to_entrez.txt | cut -f2 | sort -u)); do
cat ../../import/SNACKKSS/drug_f1_matches/pos/$targ.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' -a 1 -a 2 -o auto -e 0 - <(cat ../../import/SNACKKSS/drug_f1_matches/neg/$targ.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {if($3 * $4 <= 0){next} if($3 < 0){print $1 "\t" $2 "\t" $3} else{print $1 "\t" $2 "\t" $4}}' | sort -u
done > SNACKKSS/f1_matches/drug_both.txt &
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#Link the relationship predictions to accuracy estimates
for pert in gene drug; do
for db in SNACKKSS ConnectivityMap; do
for preds in f1_matches predictions; do
if [ $(echo $pert $preds | grep drug | grep predictions | wc -l) -gt 0 ]; then continue; fi
predspref=$preds
if [ $(echo $db | grep ConnectivityMap | wc -l) -gt 0 ]; then predspref=id_$preds; fi
for sign in pos neg; do
comp='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then comp='<'; fi
for targ in $(cat $db/$predspref/${pert}_both.txt | cut -f2 | sort -u); do cat $db/$predspref/${pert}_both.txt | awk 'BEGIN {FS = "\t"} $2 == "'$targ'" && $3 '$comp' 0 {print ".\t" $0}' | join -t$'\t' - <(cat ${pert}_scoretables/${db}_${preds}_${sign}.txt | awk 'BEGIN {FS = "\t"} {print ".\t" $0}') | cut -f2- | awk 'sqrt($3^2) > $4' | sort -k5,5gr | sort -k1,1 -k2,2 -u | cut -f1,2,3,5; done > $db/$predspref/${pert}_${sign}_estimated.txt &
done
done
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

