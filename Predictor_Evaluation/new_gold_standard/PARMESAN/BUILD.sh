for sign in pos neg; do
comp='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then comp='<'; fi
cat ../../import/PARMESAN/gene_consensus.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat ../../import/reactome.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} $1 '$comp' 0 {print (0.5 + (0.5*$1*$2/sqrt($1^2))) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; cur = ""; p = 0; n = 0} {if(cur != "" && cur != int($2*10)/10){print cur "\t" p /(p+n+1)} cur = int($2*10)/10; if($1 > 0){p++} else{n++}} END {print "0\t" p/(p+n+1)}' | sort -k1,1g -k2,2gr | awk 'BEGIN {FS = "\t"; cur = 0} $2 > cur{print; cur = $2}' > ${sign}_consensus_scoretable.txt
cat ../../import/PARMESAN/gene_consensus.txt | cut -f-3 | awk 'NR > 1 && $3 '$comp' 0' | sort -k1,1 | join -t$'\t' - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | cut -f1 | sort -u) | awk '{print ".\t" $0}' | join -t$'\t' - <(cat ${sign}_consensus_scoretable.txt | awk '{print ".\t" $0}') | cut -f2- | awk 'sqrt($3^2) > $4' | cut -f1,2,3,5 | sort -k4,4gr | sort -k1,1 -k2,2 -u > ${sign}_consensus_preds.txt &
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

for sign in pos neg; do
comp='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then comp='<'; fi
cat ../../import/PARMESAN/gene_hypotheses.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat ../../import/reactome.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} $1 '$comp' 0 {print (0.5 + (0.5*$1*$2/sqrt($1^2))) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; cur = ""; p = 0; n = 0} {if(cur != "" && cur != int($2*10)/10){print cur "\t" p /(p+n+1)} cur = int($2*10)/10; if($1 > 0){p++} else{n++}} END {print "0\t" p/(p+n+1)}' | sort -k1,1g -k2,2gr | awk 'BEGIN {FS = "\t"; cur = 0} $2 > cur{print; cur = $2}' > ${sign}_indirect_scoretable.txt
done

for target in $(cat pos_consensus_preds.txt neg_consensus_preds.txt | cut -f2 | sort -u | comm -12 - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | cut -f2 | sort -u)); do
cat pos_consensus_preds.txt neg_consensus_preds.txt | cut -f-3 | sort -k1,1 | join -t$'\t' - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | awk '$2 == "'$target'"' | cut -f1 | sort -u) | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(cat pos_consensus_preds.txt neg_consensus_preds.txt | cut -f-3 | awk '$2 == "'$target'"' | cut -f1,3 | sort -u | sort -k1,1) | cut -f2- | sort -k1,1 | awk 'BEGIN {FS = "\t"; cur = ""; p = 0; n = 0} $2 * $3 != 0 {if(cur != "" && cur != $1){print cur "\t'$target'\t" (p-n)*sqrt((p-n)^2)/(p+n); p = 0; n = 0} cur = $1; x = $2 * $3; if(x > 0){p += x} else if(x < 0){n -= x}} END {print cur "\t'$target'\t" (p-n)*sqrt((p-n)^2)/(p+n)}'
done > indirect.txt

for sign in pos neg; do
comp='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then comp='<'; fi
cat indirect.txt | awk 'NR > 1 && $3 '$comp' 0 {print ".\t" $0}' | join -t$'\t' - <(cat ${sign}_indirect_scoretable.txt | awk '{print ".\t" $0}') | cut -f2- | awk 'sqrt($3^2) > $4' | cut -f1,2,3,5 | sort -k4,4gr | sort -k1,1 -k2,2 -u > ${sign}_indirect_preds.txt
done

for sign in pos neg; do
comp='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then comp='<'; fi
cat ../../import/PARMESAN/gene_archs4_predictions.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat ../../import/reactome.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} $1 '$comp' 0 {print (0.5 + (0.5*$1*$2/sqrt($1^2))) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; cur = ""; p = 0; n = 0} {if(cur != "" && cur != int($2*10)/10){print cur "\t" p /(p+n+1)} cur = int($2*10)/10; if($1 > 0){p++} else{n++}} END {print "0\t" p/(p+n+1)}' | sort -k1,1g -k2,2gr | awk 'BEGIN {FS = "\t"; cur = 0} $2 > cur{print; cur = $2}' > ${sign}_pa4_scoretable.txt
done

cat pos_consensus_preds.txt neg_consensus_preds.txt | cut -f2 | sort -u | join -t$'\t' -1 1 -2 2 - <(cat ../../import/archs4/human_column_gene_entrez.txt | cut -f2- | sort -t$'\t' -k2,2) | cut -f2 | sort -u | awk '{print "#GENE#" $1 "#GENE#"}' | sort -u | grep -wf - <(cat ../../import/archs4/human_correlation_table.txt | awk '{print "#GENE#" $1 "#GENE#\t" $0}') | cut -f2- | sort -u | sort -k1,1 > shortened_human_correlation_table.txt

mkdir pa4
for target in $(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | cut -f2 | sort -u | comm -12 - <(cat ../../import/archs4/human_column_gene_entrez.txt | cut -f3 | sort -u)); do
index=$(cat ../../import/archs4/human_column_gene_entrez.txt | awk '$3 == "'$target'"' | cut -f1 | sort -gu | head -1)
while [ $(jobs | grep Running | wc -l) -gt 30 ]; do jobs; sleep 1; done
cat pos_consensus_preds.txt neg_consensus_preds.txt | cut -f-3 | sort -k1,1 | join -t$'\t' - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | awk '$2 == "'$target'"' | cut -f1 | sort -u) | sort -k2,2 | join -t$'\t' -1 2 -2 2 <(cat ../../import/archs4/human_column_gene_entrez.txt | cut -f2- | sort -k2,2) - | cut -f2- | sort -k1,1 | join -t$'\t' - <(cat shortened_human_correlation_table.txt | cut -f1,$index | sort -k1,1) | cut -f2- | sort -k1,1 | awk 'BEGIN {FS = "\t"; cur = ""; p = 0; n = 0} $2 * $3 != 0 {if(cur != "" && cur != $1){print cur "\t'$target'\t" (p-n)*sqrt((p-n)^2)/(p+n); p = 0; n = 0} cur = $1; x = $2 * $3; if(x > 0){p += x} else if(x < 0){n -= x}} END {print cur "\t'$target'\t" (p-n)*sqrt((p-n)^2)/(p+n)}' > pa4/$target.txt &
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done
cat pa4/* | sort -u > pa4.txt
rm -r pa4 shortened_human_correlation_table.txt


for sign in pos neg; do
comp='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then comp='<'; fi
cat pa4.txt | awk 'NR > 1 && $3 '$comp' 0 {print ".\t" $0}' | join -t$'\t' - <(cat ${sign}_pa4_scoretable.txt | awk '{print ".\t" $0}') | cut -f2- | awk 'sqrt($3^2) > $4' | cut -f1,2,3,5 | sort -k4,4gr | sort -k1,1 -k2,2 -u > ${sign}_pa4_preds.txt
done


