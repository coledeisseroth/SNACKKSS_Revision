#Translate the coexpression matrix from ARCHS4 into a set of relationship predictions, for humans and for mice.
mkdir human_a4c_preds mouse_a4c_preds
for target in $(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | cut -f2 | sort -u | comm -12 - <(cat ../../import/archs4/human_column_gene_entrez.txt | cut -f3 | sort -u)); do
while [ $(jobs | grep Running | wc -l) -gt 30 ]; do jobs; sleep 1; done
index=$(cat ../../import/archs4/human_column_gene_entrez.txt | awk '$3 == '$target | cut -f1 | sort -gu | head -1)
cat ../../import/archs4/human_correlation_table.txt | cut -f1,$index | sort -k1,1 | join -t$'\t' - <(cat ../../import/archs4/human_column_gene_entrez.txt | cut -f2- | sort -k1,1) | cut -f2- | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | awk '$2 == "'$target'"' | cut -f1 | sort -u) | awk '{print $1 "\t'$target'\t" $2}' > human_a4c_preds/$target.txt &
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done
cat human_a4c_preds/* | sort -u > human_a4c_preds.txt
rm -r human_a4c_preds/

for target in $(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | cut -f2 | sort -u | join -t$'\t' - <(cat ../../import/diopt_human_mouse_1to1.txt | sort -k1,1) | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(cat ../../import/archs4/mouse_column_gene_entrez.txt | cut -f3 | sort -u) | sed 's/\t/_/g' | sort -u); do
while [ $(jobs | grep Running | wc -l) -gt 30 ]; do jobs; sleep 1; done
humantarg=$(echo $target | cut -d_ -f2)
mousetarg=$(echo $target | cut -d_ -f1)
index=$(cat ../../import/archs4/mouse_column_gene_entrez.txt | awk '$3 == '$mousetarg | cut -f1 | sort -gu | head -1)
cat ../../import/archs4/mouse_correlation_table.txt | cut -f1,$index | sort -k1,1 | join -t$'\t' - <(cat ../../import/archs4/mouse_column_gene_entrez.txt | cut -f2- | sort -k1,1) | cut -f2- | sort -k2,2 | join -t$'\t' -1 2 -2 2 <(cat ../../import/diopt_human_mouse_1to1.txt | sort -k2,2) - | cut -f2- | sort -k1,1 | join -t$'\t' - <(cat ../corpora/trrust/noreactome_mod_to_entrez.txt | awk '$2 == "'$humantarg'"' | cut -f1 | sort -u) | awk '{print $1 "\t'$humantarg'\t" $2}' > mouse_a4c_preds/$humantarg.txt &
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done
cat mouse_a4c_preds/* | sort -u > mouse_a4c_preds.txt
rm -r mouse_a4c_preds

for species in human mouse; do
for sign in pos neg; do
comp='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then comp='<'; fi
cat ../../import/archs4/${species}_archs4_correlation_preds.txt | awk '$3 '$comp' 0 {print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat ../../import/reactome.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} {print 0.5 + (0.5 * $1 * $2 / sqrt($1^2)) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} {if(cur != "" && cur != $2){print cur "\t" p/(p+n+1)} cur = $2; if($1 == 1){p++} else{n++}} END {print "0\t" p / (p+n+1)}' | sort -k1,1g -k2,2gr | awk 'BEGIN {FS = "\t"; cur = 0} $2 > cur{print; cur = $2}' > ${species}_${sign}_scoretable.txt
done
done

for species in human mouse; do
for sign in pos neg; do
comp='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then comp='<'; fi
for targ in $(cat ${species}_a4c_preds.txt | cut -f2 | sort -u); do cat ${species}_a4c_preds.txt | awk 'BEGIN {FS = "\t"} $2 == "'$targ'" && $3 '$comp' 0 {print ".\t" $0}' | join -t$'\t' - <(cat ${species}_${sign}_scoretable.txt | awk 'BEGIN {FS = "\t"} {print ".\t" $0}') | cut -f2- | awk 'sqrt($3^2) > $4' | sort -k5,5gr | sort -k1,1 -k2,2 -u | cut -f1,2,3,5; done > ${species}_${sign}_estimated.txt
done
done

