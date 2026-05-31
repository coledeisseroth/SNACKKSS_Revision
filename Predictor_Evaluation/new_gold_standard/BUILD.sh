cd corpora
bash BUILD.sh
cd ../signature
bash BUILD.sh
cd ../correlations
bash BUILD.sh
cd ../PARMESAN
bash BUILD.sh
cd ../PubTator3
bash BUILD.sh
cd ..

#Run the comparisons with the new drug-gene relations
mkdir drug_comparisons
for sign in 1 -1; do
signprefix=pos
if [ $sign -lt 0 ]; then signprefix=neg; fi
for db in PARMESAN PubTator3; do
cat ../import/$db/${signprefix}_consensus_preds.txt | awk '{print $1 "_" $2 "\t" $4}' | sort -k1,1 | join -t$'\t' - <(cat corpora/nodgidb_aggregated_mod_to_entrez.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | sort -u > drug_comparisons/${db}_consensus_${sign}.txt &
for preds in indirect pa4; do
predspref=$preds
if [ $(echo $db $preds | grep PubTator3 | grep pa4 | wc -l) -gt 0 ]; then predspref=p3a4; fi
for targ in $(cat corpora/nodgidb_aggregated_mod_to_entrez.txt | cut -f2 | sort -u | comm -12 - <(ls ../import/$db/${signprefix}_${predspref}_preds/ | cut -d. -f1 | sort -u)); do cat ../import/$db/${signprefix}_${predspref}_preds/$targ.txt | cut -f1,4 | sort -k1,1 | join -t$'\t' - <(cat corpora/nodgidb_aggregated_mod_to_entrez.txt | awk '$2 == "'$targ'" {print $1 "\t" $3}' | sort -k1,1) | awk '{print $1 "\t'$targ'\t" $2 "\t" $3}' | sort -u; done > drug_comparisons/${db}_${predspref}_${sign}.txt &
done
done
for db in SNACKKSS ConnectivityMap; do
predspref=f1_matches
if [ $(echo $db | grep ConnectivityMap | wc -l) -gt 0 ]; then predspref=id_f1_matches; fi
for targ in $(cat corpora/nodgidb_aggregated_mod_to_entrez.txt | cut -f2 | sort -u | comm -12 - <(cat signature/$db/$predspref/drug_${signprefix}_estimated.txt | cut -f2 | sort -u)); do cat signature/$db/$predspref/drug_${signprefix}_estimated.txt | awk '$2 == "'$targ'"' | cut -f1,4 | sort -k1,1 | join -t$'\t' - <(cat corpora/nodgidb_aggregated_mod_to_entrez.txt | awk '$2 == "'$targ'" {print $1 "\t" $3}' | sort -k1,1) | awk '{print $1 "\t'$targ'\t" $2 "\t" $3}' | sort -u; done > drug_comparisons/${db}_f1_matches_${sign}.txt &
for targ in $(cat corpora/nodgidb_aggregated_mod_to_entrez.txt | cut -f2 | sort -u | comm -12 - <(ls ../import/$db/drug_predictions_${signprefix}_estimated/ | cut -d. -f1 | sort -u)); do cat ../import/$db/drug_predictions_${signprefix}_estimated/$targ.txt | cut -f1,4 | sort -k1,1 | join -t$'\t' - <(cat corpora/nodgidb_aggregated_mod_to_entrez.txt | awk '$2 == "'$targ'" {print $1 "\t" $3}' | sort -k1,1) | awk '{print $1 "\t'$targ'\t" $2 "\t" $3}' | sort -u; done > drug_comparisons/${db}_predictions_${sign}.txt &
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#Run the comparisons with the new gene-gene relations
mkdir gene_comparisons
for sign in 1 -1; do
signprefix=pos
if [ $sign -lt 0 ]; then signprefix=neg; fi
for db in PARMESAN PubTator3; do
cat $db/${signprefix}_consensus_preds.txt | awk '{print $1 "_" $2 "\t" $4}' | sort -k1,1 | join -t$'\t' - <(cat corpora/trrust/noreactome_mod_to_entrez.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | sort -u > gene_comparisons/${db}_consensus_${sign}.txt &
for preds in indirect pa4; do
predspref=$preds
if [ $(echo $db $preds | grep PubTator3 | grep pa4 | wc -l) -gt 0 ]; then predspref=p3a4; fi
cat $db/${signprefix}_${predspref}_preds.txt | awk '{print $1 "_" $2 "\t" $4}' | sort -k1,1 | join -t$'\t' - <(cat corpora/trrust/noreactome_mod_to_entrez.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | sort -u > gene_comparisons/${db}_${predspref}_${sign}.txt &
done
done
for db in SNACKKSS ConnectivityMap; do
for preds in f1_matches predictions; do
predspref=$preds
if [ $(echo $db | grep ConnectivityMap | wc -l) -gt 0 ]; then predspref=id_$preds; fi
cat signature/$db/$predspref/${signprefix}_estimated.txt | awk '$3 * ('$sign') > 0 {print $1 "_" $2 "\t" $4}' | sort -k1,1 | join -t$'\t' - <(cat corpora/trrust/noreactome_mod_to_entrez.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | sort -u > gene_comparisons/${db}_${preds}_${sign}.txt &
done
done
for species in human mouse; do
cat correlations/${species}_${signprefix}_estimated.txt | awk '$3 * ('$sign') > 0 {print $1 "_" $2 "\t" $4}' | sort -k1,1 | join -t$'\t' - <(cat corpora/trrust/noreactome_mod_to_entrez.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | sort -u > gene_comparisons/${species}_a4c_${preds}_${sign}.txt &
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

for pert in gene drug; do
(ls ${pert}_comparisons/ | rev | cut -d_ -f2- | rev | sort -u | awk 'BEGIN {print} {print $0 "_pos"; print $0 "_neg"}' | paste -sd$'\t' | sed 's/_f1_matches/ signature-matching/g' | sed 's/ConnectivityMap_predictions/CMA4/g' | sed 's/SNACKKSS_predictions/SA4/g' | sed 's/_pos/ supportive/g' | sed 's/_neg/ inhibitory/g' | sed 's/_a4c_predictions/ A4C/g' | sed 's/human/Human/g' | sed 's/mouse/Mouse/g' | sed 's/PubTator3_p3a4/P3A4/g' | sed 's/PARMESAN_pa4/PA4/g' | sed 's/_indirect/ indirect/g' | sed 's/_consensus/ consensus/g' 
for i in $(seq 0 100 | awk '{print $1 / 100}'); do
for db in $(ls ${pert}_comparisons/ | rev | cut -d_ -f2- | rev | sort -u); do
for sign in 1 -1; do
cat ${pert}_comparisons/${db}_${sign}.txt | cut -f3- | awk 'BEGIN {FS = "\t"} {print ((($2 * ('$sign')) + 1) / 2) "\t" $1}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} {if(cur != "" && cur != $2 && p + n > 0){print cur "\t" p "\t" n "\t" p/(p+n)} cur = $2; if($1 == 1){p++} else{n++}} END {print cur "\t" p "\t" n "\t" p/(p+n)}' | cut -f2,4 | awk 'BEGIN {print 0} $2 > '$i' {print $1}' | sort -gr | head -1
done
done | paste -sd$'\t' | awk '{print "'$i'\t" $0}'
done) > ${pert}_unsmoothed_precision_count.txt &
(ls ${pert}_comparisons/ | rev | cut -d_ -f2- | rev | sort -u | awk 'BEGIN {print} {print $0 "_pos"; print $0 "_neg"}' | paste -sd$'\t' | sed 's/_f1_matches/ signature-matching/g' | sed 's/ConnectivityMap_predictions/CMA4/g' | sed 's/SNACKKSS_predictions/SA4/g' | sed 's/_pos/ supportive/g' | sed 's/_neg/ inhibitory/g' | sed 's/_a4c_predictions/ A4C/g' | sed 's/human/Human/g' | sed 's/mouse/Mouse/g' | sed 's/PubTator3_p3a4/P3A4/g' | sed 's/PARMESAN_pa4/PA4/g' | sed 's/_indirect/ indirect/g' | sed 's/_consensus/ consensus/g'
for i in $(seq 0 100 | awk '{print $1 / 100}'); do
for db in $(ls ${pert}_comparisons/ | rev | cut -d_ -f2- | rev | sort -u); do
for sign in 1 -1; do
cat ${pert}_comparisons/${db}_${sign}.txt | cut -f3- | awk 'BEGIN {FS = "\t"} {print ((($2 * ('$sign')) + 1) / 2) "\t" $1}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} {if(cur != "" && cur != $2){print cur "\t" p "\t" n "\t" p/(p+n+1)} cur = $2; if($1 == 1){p++} else{n++}} END {print cur "\t" p "\t" n "\t" p/(p+n+1)}' | cut -f2,4 | awk 'BEGIN {print 0} $2 > '$i' {print $1}' | sort -gr | head -1
done
done | paste -sd$'\t' | awk '{print "'$i'\t" $0}'
done) > ${pert}_smoothed_precision_count.txt &
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done


for pert in gene drug; do
for db in $(ls ${pert}_comparisons/ | rev | cut -d_ -f2- | rev | sort -u); do
for sign in 1 -1; do
python3 ../src/logrank.py <(cat ${pert}_comparisons/${db}_${sign}.txt | cut -f3- | awk 'BEGIN {FS = "\t"} {print ((($2 * ('$sign')) + 1) / 2) "\t" $1}') | awk '{print "'$pert'\t'$db'\t'$sign'\t" $0}'
done
done
done > predictor_logrank.txt

cd ablation_test
bash BUILD.sh
cd ..
