#Match signatures with Spearman correlations
for pair in $(cat ../import/dgidb.txt | cut -f-2 | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(ls ../import/SNACKKSS/gene_nested_z/ | cut -d. -f1 | sort -u) | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(ls ../import/SNACKKSS/drug_nested_z/ | cut -d. -f1 | sort -u) | sed 's/\t/_/g' | sort -u); do
mod=$(echo $pair | cut -d_ -f1)
targ=$(echo $pair | cut -d_ -f2)
for i in $(seq 0 9 | awk '{print $1 / 10}'); do
python3 ../src/spearman.py <(cat ../import/SNACKKSS/gene_nested_z/$targ.txt | awk 'sqrt($2^2) > '$i | sort -k1,1 | join -t$'\t' - <(cat ../import/SNACKKSS/drug_nested_z/$mod.txt | awk 'sqrt($2^2) > '$i | sort -k1,1) | cut -f2-) | awk 'BEGIN {FS = "\t"} {print "'$i'\t'$mod'\t'$targ'\t" $0}'
done
done > drug_spearman.txt &
for pair in $(cat ../import/reactome.txt | cut -f-2 | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(ls ../import/SNACKKSS/gene_nested_z/ | cut -d. -f1 | sort -u) | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(ls ../import/SNACKKSS/gene_nested_z/ | cut -d. -f1 | sort -u) | sed 's/\t/_/g' | sort -u); do
mod=$(echo $pair | cut -d_ -f1)
targ=$(echo $pair | cut -d_ -f2)
for i in $(seq 0 9 | awk '{print $1 / 10}'); do
python3 ../src/spearman.py <(cat ../import/SNACKKSS/gene_nested_z/$targ.txt | awk 'sqrt($2^2) > '$i | sort -k1,1 | join -t$'\t' - <(cat ../import/SNACKKSS/gene_nested_z/$mod.txt | awk 'sqrt($2^2) > '$i | sort -k1,1) | cut -f2-) | awk 'BEGIN {FS = "\t"} {print "'$i'\t'$mod'\t'$targ'\t" $0}'
done
done > gene_spearman.txt &
#Match signatures with DF1
for pair in $(cat ../import/dgidb.txt | cut -f-2 | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(ls ../import/SNACKKSS/gene_nested_z/ | cut -d. -f1 | sort -u) | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(ls ../import/SNACKKSS/drug_nested_z/ | cut -d. -f1 | sort -u) | sed 's/\t/_/g' | sort -u); do
mod=$(echo $pair | cut -d_ -f1)
targ=$(echo $pair | cut -d_ -f2)
for i in $(seq 0 9 | awk '{print $1 / 10}'); do
cat ../import/SNACKKSS/drug_nested_z/$mod.txt | awk 'sqrt($2^2) > '$i | sort -k1,1 | join -t$'\t' -a 1 -a 2 -o auto -e 0 - <(cat ../import/SNACKKSS/gene_nested_z/$targ.txt | awk 'sqrt($2^2) > '$i | sort -k1,1) | awk 'BEGIN {FS = "\t"; ptp = 0; pfp = 0; pfn = 0; ntp = 0; nfp = 0; nfn = 0} {if(sqrt($3^2) <= '$i' && sqrt($2^2) <= '$i'){next} else if(sqrt($3^2) <= '$i'){pfp++; nfp++} else if(sqrt($2^2) <= '$i'){pfn++; nfn++} else if($2 * $3 > 0){ptp++; nfp++} else if($2 * $3 < 0){ntp++; pfp++}} END {pprec = ptp / (ptp+pfp+1); prec = ptp / (ptp+pfn+1); nprec = ntp / (ntp+nfp+1); nrec = ntp / (ntp+nfn+1); pf1 = 3 * pprec * prec / (pprec + prec+1); nf1 = 3 * nprec * nrec / (nprec + nrec+1); print "'$i'\t'$mod'\t'$targ'\t" (pf1-nf1)*2*sqrt((pf1-nf1)^2)/(pf1+nf1+1)}'
done
done > drug_df1.txt &
for pair in $(cat ../import/reactome.txt | cut -f-2 | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(ls ../import/SNACKKSS/gene_nested_z/ | cut -d. -f1 | sort -u) | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(ls ../import/SNACKKSS/gene_nested_z/ | cut -d. -f1 | sort -u) | sed 's/\t/_/g' | sort -u); do
mod=$(echo $pair | cut -d_ -f1)
targ=$(echo $pair | cut -d_ -f2)
for i in $(seq 0 9 | awk '{print $1 / 10}'); do
cat ../import/SNACKKSS/gene_nested_z/$mod.txt | awk 'sqrt($2^2) > '$i | sort -k1,1 | join -t$'\t' -a 1 -a 2 -o auto -e 0 - <(cat ../import/SNACKKSS/gene_nested_z/$targ.txt | awk 'sqrt($2^2) > '$i | sort -k1,1) | awk 'BEGIN {FS = "\t"; ptp = 0; pfp = 0; pfn = 0; ntp = 0; nfp = 0; nfn = 0} {if(sqrt($3^2) <= '$i' && sqrt($2^2) <= '$i'){next} else if(sqrt($3^2) <= '$i'){pfp++; nfp++} else if(sqrt($2^2) <= '$i'){pfn++; nfn++} else if($2 * $3 > 0){ptp++; nfp++} else if($2 * $3 < 0){ntp++; pfp++}} END {pprec = ptp / (ptp+pfp+1); prec = ptp / (ptp+pfn+1); nprec = ntp / (ntp+nfp+1); nrec = ntp / (ntp+nfn+1); pf1 = 3 * pprec * prec / (pprec + prec+1); nf1 = 3 * nprec * nrec / (nprec + nrec+1); print "'$i'\t'$mod'\t'$targ'\t" (pf1-nf1)*2*sqrt((pf1-nf1)^2)/(pf1+nf1+1)}'
done
done > gene_df1.txt &
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#LOOCV with the signature-matchers
mkdir -p loo
for pert in gene drug; do
groundtruth=../import/dgidb.txt
compsgn='-'
if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then groundtruth=../import/reactome.txt; compsgn='+'; fi
for sign in 1 -1; do
for mod in $(cat ${pert}_spearman.txt | cut -f2 | sort -u); do
champ=$(for i in $(seq 0 9 | awk '{print $1 / 10}'); do
python3 ../src/logrank.py <(cat ${pert}_spearman.txt | awk '$1 == '$i' && $2 != "'$mod'" && $4 * ('$sign') > 0' | cut -f2- | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(cat $groundtruth | awk '$1 == "'$mod'"' | cut -f2 | sort -u | comm -13 - <(cat ${pert}_spearman.txt | cut -f3 | sort -u)) | awk '{print $2 "_" $1 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | awk 'BEGIN {FS = "\t"} {print (0.5 '$compsgn' ($2 * $3 * 0.5 / sqrt($2^2))) "\t" sqrt($2^2)}') | awk '{print "'$i'\t" $1}'
done | grep -v nan | sort -k2,2gr | head -1 | cut -f1)
cat ${pert}_spearman.txt | awk '$1 == "'$champ'" && $2 == "'$mod'" && $4 * ('$sign') > 0 {print "'$champ'\t" $0}' | cut -f2-
done > loo/${pert}_${sign}.txt &
for mod in $(cat ${pert}_df1.txt | cut -f2 | sort -u); do
champ=$(for i in $(seq 0 9 | awk '{print $1 / 10}'); do
python3 ../src/logrank.py <(cat ${pert}_df1.txt | awk '$1 == '$i' && $2 != "'$mod'" && $4 * ('$sign') > 0' | cut -f2- | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(cat $groundtruth | awk '$1 == "'$mod'"' | cut -f2 | sort -u | comm -13 - <(cat ${pert}_df1.txt | cut -f3 | sort -u)) | awk '{print $2 "_" $1 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | awk 'BEGIN {FS = "\t"} {print (0.5 '$compsgn' ($2 * $3 * 0.5 / sqrt($2^2))) "\t" sqrt($2^2)}') | awk '{print "'$i'\t" $1}'
done | grep -v nan | sort -k2,2gr | head -1 | cut -f1)
cat ${pert}_df1.txt | awk '$1 == "'$champ'" && $2 == "'$mod'" && $4 * ('$sign') > 0 {print "'$champ'\t" $0}' | cut -f2-
done > loo/df1_${pert}_${sign}.txt &
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#Calculate a larger set of matches for the SA4 procedure
mkdir linkable_matches
for pert in gene drug; do
groundtruth=../import/dgidb.txt
compsgn='-'
if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then groundtruth=../import/reactome.txt; compsgn='+'; fi
for sign in 1 -1; do
mkdir linkable_matches/${pert}_${sign}
mkdir linkable_matches/df1_${pert}_${sign}
for mod in $(ls ../import/SNACKKSS/${pert}_nested_z | cut -d. -f1 | sort -u | comm -12 - <(cat $groundtruth | cut -f1 | sort -u)); do
while [ $(jobs | grep Running | wc -l) -gt 10 ]; do jobs; sleep 1; done
champ=$(for i in $(seq 0 9 | awk '{print $1 / 10}'); do
python3 ../src/logrank.py <(cat ${pert}_spearman.txt | awk '$1 == '$i' && $2 != "'$mod'" && $4 * ('$sign') > 0' | cut -f2-4 | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(cat $groundtruth | awk '$1 == "'$mod'"' | cut -f2 | sort -u | comm -13 - <(cat ${pert}_spearman.txt | cut -f3 | sort -u)) | awk '{print $2 "_" $1 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | awk 'BEGIN {FS = "\t"} {print (0.5 '$compsgn' ($2 * $3 * 0.5 / sqrt($2^2))) "\t" sqrt($2^2)}') | awk '{print "'$i'\t" $1}'
done | grep -v nan | sort -k2,2gr | head -1 | cut -f1)
for intermediate in $(ls ../import/SNACKKSS/gene_nested_z/ | cut -d. -f1 | sort -u); do
if [ $(echo $pert $mod $intermediate | awk '$1 == "gene" && $2 == $3' | wc -l) -gt 0 ]; then continue; fi
python3 ../src/spearman.py <(cat ../import/SNACKKSS/gene_nested_z/$intermediate.txt | awk 'sqrt($2^2) > '$champ | sort -k1,1 | join -t$'\t' - <(cat ../import/SNACKKSS/${pert}_nested_z/$mod.txt | awk 'sqrt($2^2) > '$champ | sort -k1,1) | cut -f2-) | awk 'BEGIN {FS = "\t"} $1 != 0 {print "'$champ'\t'$mod'\t'$intermediate'\t" $0}'
done > linkable_matches/${pert}_${sign}/$mod.txt &
done
for mod in $(ls ../import/SNACKKSS/${pert}_nested_z | cut -d. -f1 | sort -u | comm -12 - <(cat $groundtruth | cut -f1 | sort -u)); do
while [ $(jobs | grep Running | wc -l) -gt 10 ]; do jobs; sleep 1; done
champ=$(for i in $(seq 0 9 | awk '{print $1 / 10}'); do
python3 ../src/logrank.py <(cat ${pert}_df1.txt | awk '$1 == '$i' && $2 != "'$mod'" && $4 * ('$sign') > 0' | cut -f2- | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(cat $groundtruth | awk '$1 == "'$mod'"' | cut -f2 | sort -u | comm -13 - <(cat ${pert}_df1.txt | cut -f3 | sort -u)) | awk '{print $2 "_" $1 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | awk 'BEGIN {FS = "\t"} {print (0.5 '$compsgn' ($2 * $3 * 0.5 / sqrt($2^2))) "\t" sqrt($2^2)}') | awk '{print "'$i'\t" $1}'
done | grep -v nan | sort -k2,2gr | head -1 | cut -f1)
for intermediate in $(ls ../import/SNACKKSS/gene_nested_z/ | cut -d. -f1 | sort -u); do
if [ $(echo $pert $mod $intermediate | awk '$1 == "gene" && $2 == $3' | wc -l) -gt 0 ]; then continue; fi
cat ../import/SNACKKSS/${pert}_nested_z/$mod.txt | awk 'sqrt($2^2) > '$champ | sort -k1,1 | join -t$'\t' -a 1 -a 2 -o auto -e 0 - <(cat ../import/SNACKKSS/gene_nested_z/$intermediate.txt | awk 'sqrt($2^2) > '$champ | sort -k1,1) | awk 'BEGIN {FS = "\t"; ptp = 0; pfp = 0; pfn = 0; ntp = 0; nfp = 0; nfn = 0} {if(sqrt($3^2) <= '$champ' && sqrt($2^2) <= '$champ'){next} else if(sqrt($3^2) <= '$champ'){pfp++; nfp++} else if(sqrt($2^2) <= '$champ'){pfn++; nfn++} else if($2 * $3 > 0){ptp++; nfp++} else if($2 * $3 < 0){ntp++; pfp++}} END {pprec = ptp / (ptp+pfp+1); prec = ptp / (ptp+pfn+1); nprec = ntp / (ntp+nfp+1); nrec = ntp / (ntp+nfn+1); pf1 = 3 * pprec * prec / (pprec + prec+1); nf1 = 3 * nprec * nrec / (nprec + nrec+1); print "'$champ'\t'$mod'\t'$intermediate'\t" (pf1-nf1)*2*sqrt((pf1-nf1)^2)/(pf1+nf1+1)}'
done > linkable_matches/df1_${pert}_${sign}/$mod.txt &
done
done
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#Gather the expression correlations needed to run SA4
mkdir a4c_confirmable
for pair in $(cat ../import/archs4/human_column_gene_entrez.txt | cut -f1,3 | sort -k2,2 | join -t$'\t' -1 2 -2 1 - <(cat ../import/reactome.txt ../import/dgidb.txt | cut -f2 | sort -u | comm -23 - <(ls a4c_confirmable/ | cut -d. -f1 | sort -u)) | sed 's/\t/_/g' | sort -u); do
while [ $(jobs | grep Running | wc -l) -gt 40 ]; do jobs; sleep 1; done
gene=$(echo $pair | cut -d_ -f1)
column=$(echo $pair | cut -d_ -f2)
cat ../import/archs4/human_correlation_table.txt | cut -f1,$column | sort -t$'\t' -k1,1 | join -t$'\t' <(cat ../import/archs4/human_column_gene_entrez.txt | cut -f2- | sort -t$'\t' -k1,1) - | cut -f2- | sort -u | sort -k1,1 | join -t$'\t' - <(ls ../import/SNACKKSS/gene_nested_z | cut -d. -f1 | sort -u) | sort -t$'\t' -k1,1 > a4c_confirmable/$gene.txt &
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#Run SA4
mkdir a4c_linked_predictions
for pert in gene drug; do
groundtruth=../import/dgidb.txt
if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then groundtruth=../import/reactome.txt; fi
for mod in $(ls linkable_matches/${pert}_* | cut -d. -f1 | sort -u); do
for targ in $(cat $groundtruth | awk '$1 == "'$mod'"' | cut -f2 | sort -u | comm -12 - <(ls a4c_confirmable | cut -d. -f1 | sort -u)); do
cat linkable_matches/${pert}_1/$mod.txt | sed 's/_/\t/g' | awk '$4 > 0 {print $2 "_" $3 "\t" $4}' | sort -t$'\t' -k1,1 | join -t$'\t' -a 1 -a 2 -o auto -e 0 - <(cat linkable_matches/${pert}_-1/$mod.txt | sed 's/_/\t/g' | awk '$4 < 0 {print $2 "_" $3 "\t" $4}' | sort -t$'\t' -k1,1) | awk '$2 == 0 || $3 == 0 {print $1 "\t" $2 + $3}' | sed 's/_/\t/g' | cut -f2- | sort -t$'\t' -k1,1 | join -t$'\t' - <(cat a4c_confirmable/$targ.txt | sort -t$'\t' -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"; p = 0; n = 0} $1 * $2 != 0 {link = $1 * $2 / (sqrt($1^2) + sqrt($2^2)); if(link > 0){p += link} else{n -= link}} END {print "'$mod'\t'$targ'\t" (p - n) * sqrt((p-n)^2) / (p+n)}'
done
done > a4c_linked_predictions/$pert.txt
for mod in $(ls linkable_matches/df1_${pert}_* | cut -d. -f1 | sort -u); do
for targ in $(cat $groundtruth | awk '$1 == "'$mod'"' | cut -f2 | sort -u | comm -12 - <(ls a4c_confirmable | cut -d. -f1 | sort -u)); do
cat linkable_matches/df1_${pert}_1/$mod.txt | sed 's/_/\t/g' | awk '$4 > 0 {print $2 "_" $3 "\t" $4}' | sort -t$'\t' -k1,1 | join -t$'\t' -a 1 -a 2 -o auto -e 0 - <(cat linkable_matches/df1_${pert}_-1/$mod.txt | sed 's/_/\t/g' | awk '$4 < 0 {print $2 "_" $3 "\t" $4}' | sort -t$'\t' -k1,1) | awk '$2 == 0 || $3 == 0 {print $1 "\t" $2 + $3}' | sed 's/_/\t/g' | cut -f2- | sort -t$'\t' -k1,1 | join -t$'\t' - <(cat a4c_confirmable/$targ.txt | sort -t$'\t' -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"; p = 0; n = 0} $1 * $2 != 0 {link = $1 * $2 / (sqrt($1^2) + sqrt($2^2)); if(link > 0){p += link} else{n -= link}} END {print "'$mod'\t'$targ'\t" (p - n) * sqrt((p-n)^2) / (p+n)}'
done
done > a4c_linked_predictions/df1_$pert.txt &
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

#Run SA4, only accepting statistically significant Spearman correlations
for pert in gene drug; do
groundtruth=../import/dgidb.txt
if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then groundtruth=../import/reactome.txt; fi
for mod in $(ls linkable_matches/${pert}_* | cut -d. -f1 | sort -u); do
for targ in $(cat $groundtruth | awk '$1 == "'$mod'"' | cut -f2 | sort -u | comm -12 - <(ls a4c_confirmable | cut -d. -f1 | sort -u)); do
cat linkable_matches/${pert}_1/$mod.txt | sed 's/_/\t/g' | awk '$4 > 0 && $5 < 0.05 {print $2 "_" $3 "\t" $4}' | sort -t$'\t' -k1,1 | join -t$'\t' -a 1 -a 2 -o auto -e 0 - <(cat linkable_matches/${pert}_-1/$mod.txt | sed 's/_/\t/g' | awk '$4 < 0 && $5 < 0.05 {print $2 "_" $3 "\t" $4}' | sort -t$'\t' -k1,1) | awk '$2 == 0 || $3 == 0 {print $1 "\t" $2 + $3}' | sed 's/_/\t/g' | cut -f2- | sort -t$'\t' -k1,1 | join -t$'\t' - <(cat a4c_confirmable/$targ.txt | sort -t$'\t' -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"; p = 0; n = 0} $1 * $2 != 0 {link = $1 * $2 / (sqrt($1^2) + sqrt($2^2)); if(link > 0){p += link} else{n -= link}} END {print "'$mod'\t'$targ'\t" (p - n) * sqrt((p-n)^2) / (p+n)}'
done
done > a4c_linked_predictions/sig_$pert.txt &
done
while [ $(jobs | grep Running | wc -l) gt 0 ]; do jobs; sleep 1; done

mkdir output_stats
#How long, on average, does it take to calculate one spearman correlation, 11 in parallel?
count=$(wc -l $(ls linkable_matches/gene_1/* | head -1; ls linkable_matches/drug_1/* | head -1) | tail -1 | awk '{print $1}')
for pert in gene drug; do
for sign in 1 -1; do
echo $pert $sign
(ls -lh linkable_matches/${pert}_${sign}/ | head -12 | tail -11 | rev | cut -d' ' -f2 | rev | awk 'BEGIN {FS = ":"; t = 0} {t += ($1 * 60) + $2} END {print t / NR}'; ls -lh linkable_matches/${pert}_${sign}/ | head -23 | tail -11 | rev | cut -d' ' -f2 | rev | awk 'BEGIN {FS = ":"; t = 0} {t += ($1 * 60) + $2} END {print t / NR}') | paste -sd$'\t'
done
done | awk 'BEGIN {t = 0}{t += ($2 - $1) / (22 * '$count')} END {print t}' > output_stats/average_spearman_runtime.txt
#One DF1 score, 11 in parallel?
count=$(wc -l $(ls linkable_matches/df1_gene_${sign}/* | head -1; ls linkable_matches/df1_drug_${sign}/* | head -1) | tail -1 | awk '{print $1}')
for pert in gene drug; do
for sign in 1 -1; do
echo $pert $sign
(ls -lh linkable_matches/df1_${pert}_${sign}/ | head -12 | tail -11 | rev | cut -d' ' -f2 | rev | awk 'BEGIN {FS = ":"; t = 0} {t += ($1 * 60) + $2} END {print t / NR}'; ls -lh linkable_matches/df1_${pert}_${sign}/ | head -23 | tail -11 | rev | cut -d' ' -f2 | rev | awk 'BEGIN {FS = ":"; t = 0} {t += ($1 * 60) + $2} END {print t / NR}') | paste -sd$'\t'
done
done | awk 'BEGIN {t = 0}{t += ($2 - $1) / (22 * '$count')} END {print t}' > output_stats/average_df1_runtime.txt

#How did the predictors do?
for pert in gene drug; do
pertsgn='+'
groundtruth=../import/reactome.txt
if [ $(echo $pert | grep drug | wc -l) -gt 0 ]; then groundtruth=../import/dgidb.txt; pertsgn='-'; fi
for sgn in 1 -1; do
python3 ../src/logrank.py <(cat loo/${pert}_${sgn}.txt | awk 'BEGIN {FS = "\t"} $4 * ('$sgn') > 0 && $5 < 0.05 {print $2 "_" $3 "\t" sqrt($4^2)}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} {print int(0.5 '$pertsgn' ($1*$2*0.5)/sqrt($1^2)) "\t" sqrt($1^2)}') | awk '{print "Signature-matching\t'$pert'\tSpearman significant\t'$sgn'\t" $0}'
python3 ../src/logrank.py <(cat a4c_linked_predictions/sig_${pert}.txt | awk 'BEGIN {FS = "\t"} $3 * ('$sgn') > 0 {print $1 "_" $2 "\t" sqrt($3^2)}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} {print int(0.5 '$pertsgn' ($1*$2*0.5)/sqrt($1^2)) "\t" sqrt($1^2)}') | awk '{print "SA4\t'$pert'\tSpearman significant\t'$sgn'\t" $0}'
for db in $pert df1_${pert}; do
dbname=Spearman
if [ $(echo $db | grep df1 | wc -l) -gt 0 ]; then dbname=DF1; fi
python3 ../src/logrank.py <(cat loo/${db}_${sgn}.txt | awk 'BEGIN {FS = "\t"} $4 * ('$sgn') > 0 {print $2 "_" $3 "\t" sqrt($4^2)}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} {print int(0.5 '$pertsgn' ($1*$2*0.5)/sqrt($1^2)) "\t" sqrt($1^2)}') | awk '{print "Signature-matching\t'$pert'\t'$dbname'\t'$sgn'\t" $0}'
python3 ../src/logrank.py <(cat a4c_linked_predictions/$db.txt | awk 'BEGIN {FS = "\t"} $3 * ('$sgn') > 0 {print $1 "_" $2 "\t" sqrt($3^2)}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} {print int(0.5 '$pertsgn' ($1*$2*0.5)/sqrt($1^2)) "\t" sqrt($1^2)}') | awk '{print "SA4\t'$pert'\t'$dbname'\t'$sgn'\t" $0}'
done
done
done > output_stats/performance_logrank.txt

for i in $(seq 0 100 | awk '{print $1 / 100}'); do
for pert in gene drug; do
groundtruth=../import/reactome.txt
pertsgn='+'
if [ $(echo $pert | grep drug | wc -l) -gt 0 ]; then pertsgn='-'; groundtruth=../import/dgidb.txt; fi
for db in $pert df1_${pert} sig_${pert}; do
for sgn in '>' '<'; do
cat a4c_linked_predictions/$db.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} $1 '$sgn' 0 {print int(0.5 '$pertsgn' ($1 * $2 * 0.5 / sqrt($1^2 * $2^2))) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} {if(cur != "" && cur != $2 && p + n > 0){print cur "\t" p "\t" n "\t" p/(p+n)} cur = $2; p += $1; n += 1; n -= $1} END {print cur "\t" p "\t" n "\t" p/(p+n)}' | cut -f2,4 | awk 'BEGIN {print 0} $2 > '$i' {print $1}' | sort -gr | head -1
done
done
done | paste -sd$'\t' | awk '{print "'$i'\t" $0}'
done | awk 'BEGIN {FS = "\t"; print "\tGene Spearman positive\tGene Spearman negative\tGene DF1 positive\tGene DF1 negative\tGene Spearman significant positive\tGene Spearman significant negative\tDrug Spearman negative\tDrug Spearman positive\tDrug DF1 negative\tDrug DF1 posiive\tDrug Spearman significant negative\tDrug Spearman significant positive"} {print}' > output_stats/sa4_accuracy.txt

for i in $(seq 0 100 | awk '{print $1 / 100}'); do
for pert in gene drug; do
groundtruth=../import/reactome.txt
pertsgn='+'
if [ $(echo $pert | grep drug | wc -l) -gt 0 ]; then pertsgn='-'; groundtruth=../import/dgidb.txt; fi
for db in $pert df1_${pert} sig_${pert}; do
for sgn in '>' '<'; do
cat a4c_linked_predictions/$db.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} $1 '$sgn' 0 {print int(0.5 '$pertsgn' ($1 * $2 * 0.5 / sqrt($1^2 * $2^2))) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} {if(cur != "" && cur != $2){print cur "\t" p "\t" n "\t" p/(p+n+1)} cur = $2; p += $1; n += 1; n -= $1} END {print cur "\t" p "\t" n "\t" p/(p+n+1)}' | cut -f2,4 | awk 'BEGIN {print 0} $2 > '$i' {print $1}' | sort -gr | head -1
done
done
done | paste -sd$'\t' | awk '{print "'$i'\t" $0}'
done | awk 'BEGIN {FS = "\t"; print "\tGene Spearman positive\tGene Spearman negative\tGene DF1 positive\tGene DF1 negative\tGene Spearman significant positive\tGene Spearman significant negative\tDrug Spearman negative\tDrug Spearman positive\tDrug DF1 negative\tDrug DF1 posiive\tDrug Spearman significant negative\tDrug Spearman significant positive"} {print}' > output_stats/sa4_accuracy_smoothed.txt

for i in $(seq 0 100 | awk '{print $1 / 100}'); do
for pert in gene drug; do
groundtruth=../import/reactome.txt
pertsgn='+'
if [ $(echo $pert | grep drug | wc -l) -gt 0 ]; then pertsgn='-'; groundtruth=../import/dgidb.txt; fi
for db in $pert df1_${pert}; do
for sgn in 1 -1; do
cat loo/${db}_${sgn}.txt | awk '{print $2 "_" $3 "\t" $4}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} $1 * ('$sgn') > 0 {print int(0.5 '$pertsgn' ($1 * $2 * 0.5 / sqrt($1^2 * $2^2))) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} {if(cur != "" && cur != $2 && p + n > 0){print cur "\t" p "\t" n "\t" p/(p+n)} cur = $2; p += $1; n += 1; n -= $1} END {print cur "\t" p "\t" n "\t" p/(p+n)}' | cut -f2,4 | awk 'BEGIN {print 0} $2 > '$i' {print $1}' | sort -gr | head -1
done
done
for sgn in 1 -1; do
cat loo/${pert}_${sgn}.txt | awk '$5 < 0.05 {print $2 "_" $3 "\t" $4}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} $1 * ('$sgn') > 0 {print int(0.5 '$pertsgn' ($1 * $2 * 0.5 / sqrt($1^2 * $2^2))) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} {if(cur != "" && cur != $2 && p + n > 0){print cur "\t" p "\t" n "\t" p/(p+n)} cur = $2; p += $1; n += 1; n -= $1} END {print cur "\t" p "\t" n "\t" p/(p+n)}' | cut -f2,4 | awk 'BEGIN {print 0} $2 > '$i' {print $1}' | sort -gr | head -1
done
done | paste -sd$'\t' | awk '{print "'$i'\t" $0}'
done | awk 'BEGIN {FS = "\t"; print "\tGene Spearman positive\tGene Spearman negative\tGene DF1 positive\tGene DF1 negative\tGene Spearman significant positive\tGene Spearman significant negative\tDrug Spearman negative\tDrug Spearman positive\tDrug DF1 negative\tDrug DF1 posiive\tDrug Spearman significant negative\tDrug Spearman significant positive"} {print}' > output_stats/signature_match_accuracy.txt

for i in $(seq 0 100 | awk '{print $1 / 100}'); do
for pert in gene drug; do
groundtruth=../import/reactome.txt
pertsgn='+'
if [ $(echo $pert | grep drug | wc -l) -gt 0 ]; then pertsgn='-'; groundtruth=../import/dgidb.txt; fi
for db in $pert df1_${pert}; do
for sgn in 1 -1; do
cat loo/${db}_${sgn}.txt | awk '{print $2 "_" $3 "\t" $4}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} $1 * ('$sgn') > 0 {print int(0.5 '$pertsgn' ($1 * $2 * 0.5 / sqrt($1^2 * $2^2))) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} {if(cur != "" && cur != $2){print cur "\t" p "\t" n "\t" p/(p+n+1)} cur = $2; p += $1; n += 1; n -= $1} END {print cur "\t" p "\t" n "\t" p/(p+n+1)}' | cut -f2,4 | awk 'BEGIN {print 0} $2 > '$i' {print $1}' | sort -gr | head -1
done
done
for sgn in 1 -1; do
cat loo/${pert}_${sgn}.txt | awk '$5 < 0.05 {print $2 "_" $3 "\t" $4}' | sort -k1,1 | join -t$'\t' - <(cat $groundtruth | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} $1 * ('$sgn') > 0 {print int(0.5 '$pertsgn' ($1 * $2 * 0.5 / sqrt($1^2 * $2^2))) "\t" sqrt($1^2)}' | sort -k2,2gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} {if(cur != "" && cur != $2){print cur "\t" p "\t" n "\t" p/(p+n+1)} cur = $2; p += $1; n += 1; n -= $1} END {print cur "\t" p "\t" n "\t" p/(p+n+1)}' | cut -f2,4 | awk 'BEGIN {print 0} $2 > '$i' {print $1}' | sort -gr | head -1
done
done | paste -sd$'\t' | awk '{print "'$i'\t" $0}'
done | awk 'BEGIN {FS = "\t"; print "\tGene Spearman positive\tGene Spearman negative\tGene DF1 positive\tGene DF1 negative\tGene Spearman significant positive\tGene Spearman significant negative\tDrug Spearman negative\tDrug Spearman positive\tDrug DF1 negative\tDrug DF1 posiive\tDrug Spearman significant negative\tDrug Spearman significant positive"} {print}' > output_stats/signature_match_accuracy_smoothed.txt

