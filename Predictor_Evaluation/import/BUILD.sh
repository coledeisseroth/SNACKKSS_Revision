#Acquire the following files from a completed SNACKKSS_Eval run
cp ~/SNACKKSS_Eval/corpora/dgidb/mod_to_entrez.txt dgidb.txt
cp ~/SNACKKSS_Eval/corpora/reactome/reactome_entrez.txt reactome.txt
mkdir ConnectivityMap
cp ~/SNACKKSS_Eval/ConnectivityMap/id_mapping/*_entrez_1to1.txt ConnectivityMap/
cp -r ~/SNACKKSS_Eval/ConnectivityMap/default/kokd_nested_z_bytarget/ ConnectivityMap/gene_nested_z
cp -r ~/SNACKKSS_Eval/ensemble_predictor/ConnectivityMap/human/drug/f1_matches/ ConnectivityMap/drug_f1_matches
cp -r ~/SNACKKSS_Eval/ConnectivityMap/default/loo_prediction_accuracy_estimates/ ConnectivityMap
cp -r ~/SNACKKSS_Eval/corpora/pubchem/split_100K_clean/ pubchem_split_100K_clean
cp ~/SNACKKSS_Eval/corpora/entrez/human_aliases_tall_clean.txt entrez_human_aliases_tall_clean.txt

mkdir SNACKKSS
cp -r ~/SNACKKSS_Eval/SNACKKSS/default/gene/nested_z/ SNACKKSS/gene_nested_z
cp -r ~/SNACKKSS_Eval/ensemble_predictor/SNACKKSS/compound_nested_z/ SNACKKSS/drug_nested_z
cp -r ~/SNACKKSS_Eval/ensemble_predictor/SNACKKSS/human/drug/f1_matches/ SNACKKSS/drug_f1_matches
cp -r ~/SNACKKSS_Eval/SNACKKSS/output_stats/ SNACKKSS/output_stats
cp -r ~/SNACKKSS_Eval/SNACKKSS/default/loo_prediction_accuracy_estimates/ SNACKKSS/

for db in SNACKKSS ConnectivityMap; do
for sign in pos neg; do
cp -r ~/SNACKKSS_Eval/ensemble_predictor/$db/human/drug/predictions/${sign}_estimated/ $db/drug_predictions_${sign}_estimated
done
done

mkdir archs4
for species in human mouse; do
cp ~/SNACKKSS_Eval/correlations/${species}_archs4_loo_preds.txt archs4/
cp ~/SNACKKSS_Eval/correlations/${species}_archs4_correlation_preds.txt archs4/
cp ~/SNACKKSS_Eval/corpora/archs4/${species}_column_gene_entrez.txt archs4/
cp ~/SNACKKSS_Eval/corpora/archs4/${species}_correlation_table.txt archs4/
done

mkdir drugrepurposinghub
cd drugrepurposinghub
#Try running this--but you might have to get it manually from the website.
wget https://repo-hub.broadinstitute.org/public/data/repo-drug-annotation-20200324.txt
cd ..

cp ~/SNACKKSS_Eval/ensemble_predictor/signature_based_optima.txt .

for db in PARMESAN PubTator3; do
mkdir $db
cp -r ~/SNACKKSS_Eval/ensemble_predictor/$db/*_preds.txt $db/
cp ~/SNACKKSS_Eval/literature/$db/*loo*.txt $db/
cp -r ~/SNACKKSS_Eval/ensemble_predictor/$db/*_preds $db/
cp ~/SNACKKSS_Eval/literature/$db/gene_hypotheses.txt $db/
cp ~/SNACKKSS_Eval/literature/$db/gene_archs4_predictions.txt $db/
done
cp ~/SNACKKSS_Eval/corpora/PARMESAN/bulk_download/gene_consensus.txt PARMESAN/
cp ~/SNACKKSS_Eval/corpora/pubtator3/gene/consensus_directionality.txt PubTator3/gene_consensus.txt

for pert in gene drug; do
groundtruth=reactome.txt
if [ $(echo $pert | grep drug | wc -l) -gt 0 ]; then groundtruth=dgidb.txt; fi
for db in SNACKKSS ConnectivityMap; do
for preds in f1_matches predictions; do
predspref=$preds
if [ $(echo $db | grep ConnectivityMap | wc -l) -gt 0 ]; then predspref=id_$preds; fi
for sign in pos neg; do
sgn='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then sgn='<'; fi
champ=$(cat signature_based_optima.txt | awk '$1 == "'$db'" && $2 == "default" && $3 == "'$pert'" && $4 == "'$preds'" && $5 == "'$sign'" {print $6}')
cat ~/SNACKKSS_Eval/$db/default/$predspref/$pert/$champ.txt | awk '$3 '$sgn' 0 {print $1 "_" $2 "\t" $3}' | sort -t$'\t' -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -t$'\t' -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} {print 0.5 + (0.5 * $1 * $2 / sqrt($1^2)) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} {if(cur != "" && cur != $2){print cur "\t" p/(p+n+1)} cur = $2; if($1 == 1){p++} else{n++}} END {print "0\t" p / (p+n+1)}' | sort -k1,1g -k2,2gr | awk 'BEGIN {FS = "\t"; cur = 0} $2 > cur{print; cur = $2}' > $db/${pert}_${preds}_${sign}_scoretable.txt
done
done
done
done

cp ~/SNACKKSS_Eval/corpora/diopt/human_mouse_1to1.txt diopt_human_mouse_1to1.txt

#Get the average runtime from running DF1 on 30 cores
count=$(wc -l $(ls ~/SNACKKSS_Eval/ensemble_predictor/SNACKKSS/human/drug/f1_matches/*/* | head -1) | awk '{print $1}')
for sign in pos neg; do
(ls -lh ~/SNACKKSS_Eval/ensemble_predictor/SNACKKSS/human/drug/f1_matches/${sign}/ | head -32 | tail -31 | rev | cut -d' ' -f2 | rev | awk 'BEGIN {FS = ":"; t = 0} {t += ($1 * 60) + $2} END {print t / NR}'; ls -lh ~/SNACKKSS_Eval/ensemble_predictor/SNACKKSS/human/drug/f1_matches/${sign}/ | head -63 | tail -31 | rev | cut -d' ' -f2 | rev | awk 'BEGIN {FS = ":"; t = 0} {t += ($1 * 60) + $2} END {print t / NR}') | paste -sd$'\t'
done | awk 'BEGIN {t = 0} {t += ($2 - $1) / (62 * '$count')} END {print t}' > 30core_average_df1_runtime.txt


