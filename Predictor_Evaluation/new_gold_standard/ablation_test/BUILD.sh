for smoothing in smoothed unsmoothed; do
increment=1
if [ $(echo $smoothing | grep unsmoothed | wc -l) -gt 0 ]; then increment=0; fi
for pert in gene drug; do
for sign in 1 -1; do
(ls ../${pert}_loo_comparisons/ | grep -v f1_matches | rev | cut -d_ -f2- | rev | sort -u | paste -sd$'\t' | sed 's/ConnectivityMap_predictions/CMA4/g' | sed 's/SNACKKSS_predictions/SA4/g' | sed 's/_pos/ supportive/g' | sed 's/_neg/ inhibitory/g' | sed 's/_a4c_predictions/ A4C/g' | sed 's/human/Human/g' | sed 's/mouse/Mouse/g' | sed 's/PubTator3_p3a4/P3A4/g' | sed 's/PARMESAN_pa4/PA4/g' | sed 's/_indirect/ indirect/g' | sed 's/_consensus/ consensus/g' | awk '{print "\tNothing\t" $0}'
for i in $(seq 0 99 | awk '{print $1 / 100}'); do echo -n $i$'\t'; for excluded in everything $(ls ../${pert}_loo_comparisons/* | grep -v f1_matches | rev | cut -d'/' -f1 | cut -d_ -f2- | rev | sort -u); do for sgn in 1 -1; do cat $(ls ../${pert}_loo_comparisons/*_${sgn}.txt | grep -v f1_matches | grep -v $excluded | paste -sd' ') | awk '{print $1 "\t" $2 "\t" $3 "\t'$sgn'\t" $4}'; done | sort -k3,3gr -k4,4gr | sort -k1,1 -k2,2 -u | sort -k3,3gr | awk 'BEGIN {FS = "\t"; cur = ""; p = 0; n = 0} $4 * ('$sign') > 0 {if(cur != "" && cur != $3 && p+n>0){print cur "\t" p "\t" n "\t" p/(p+n+'$increment')} cur = $3; if(('$sign') * $5 > 0){p++} else{n++}} END {print cur "\t" p "\t" n "\t" p/(p+n+'$increment'); print "1\t0\t0\t1"}' | awk '$4 > '$i | sort -k2,2gr | head -1 | cut -f2; done | paste -sd$'\t'
done) > ${pert}_${sign}_${smoothing}_ablation_test.txt &
done
done
#Ablation statistical test: Count the number of accuracy thresholds where the coverage improved after adding each predictor
(echo $'Perturbation type\tDirection\tExcluded database\tThresholds improved\tThresholds not improved'
for pert in gene drug; do
for sign in 1 -1; do
predfiles=$(ls ../${pert}_loo_comparisons/* | grep -v f1_matches)
for excluded in nothing $(echo $predfiles | sed 's/ /\n/g' | rev | cut -d'/' -f1 | cut -d_ -f2- | rev | sort -u); do
for sgn in 1 -1; do cat $(ls $predfiles | grep _${sgn}.txt) | awk '{print $1 "\t" $2 "\t" $3 "\t'$sgn'\t" $4}'; done | sort -k3,3gr -k4,4gr | sort -k1,1 -k2,2 -u | sort -k3,3gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} $4 * ('$sign') > 0 {if(cur != "" && cur != $3){print cur "\t" p "\t" p/(p+n+'$increment')} cur = $3; if($4 * $5 > 0){p++} else{n++}} END {print cur "\t" p "\t" p/(p+n+'$increment')}' | cut -f2- | sort -k1,1gr -k2,2gr | awk 'BEGIN {prev = 0} $2 > prev {print ".\t" $0; prev = $2}' | sort -k1,1 | join -t$'\t' - <(for sgn in 1 -1; do cat $(echo $predfiles | sed 's/ /\n/g' | grep _${sgn}.txt | grep -v $excluded | paste -sd' ') | awk '{print $1 "\t" $2 "\t" $3 "\t'$sgn'\t" $4}'; done | sort -k3,3gr -k4,4gr | sort -k1,1 -k2,2 -u | sort -k3,3gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""} $4 * ('$sign') > 0 {if(cur != "" && cur != $3){print cur "\t" p "\t" p/(p+n+'$increment')} cur = $3; if($4 * $5 > 0){p++} else{n++}} END {print cur "\t" p "\t" p/(p+n+'$increment')}' | cut -f2- | sort -k1,1gr -k2,2gr | awk 'BEGIN {prev = 0} $2 > prev {print ".\t" $0; prev = $2}' | sort -k1,1) | cut -f2- | awk '$2 == 0 || $2 > $4' | sort -k1,1gr | sort -k4,4 -u | sort -k4,4gr -k1,1gr | awk 'BEGIN {FS = "\t"; cur = 0; p = 0; n = 0} {if($1 > cur){cur = $1} if(cur > $3){p++} else{n++}} END {print "'$pert'\t'$sign'\t'$excluded'\t" p "\t" n}' | sed 's/\///g' | sed 's/\.\.//g' | sed 's/gene_loo_comparisons//g' | sed 's/drug_loo_comparisons//g' | sed 's/_1.txt//g' | sed 's/_-1.txt//g' | sed 's/ConnectivityMap_predictions/CMA4/g' | sed 's/SNACKKSS_predictions/SA4/g' | sed 's/_pos/ supportive/g' | sed 's/_neg/ inhibitory/g' | sed 's/_a4c_predictions/ A4C/g' | sed 's/human/Human/g' | sed 's/mouse/Mouse/g' | sed 's/PubTator3_p3a4/P3A4/g' | sed 's/PARMESAN_pa4/PA4/g' | sed 's/_indirect/ indirect/g' | sed 's/_consensus/ consensus/g'
done
done
done) > ${smoothing}_ablation_statistics.txt &
done
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done

for pert in gene drug; do
for sign in 1 -1; do
python3 ../../src/logrank.py <(for sgn in 1 -1; do cat $(ls ../${pert}_loo_comparisons/*_${sgn}.txt | grep -v f1_matches | paste -sd' ') | awk '{print $1 "\t" $2 "\t" $3 "\t'$sgn'\t" $4}'; done | sort -k3,3gr -k4,4gr | sort -k1,1 -k2,2 -u | awk 'BEGIN {FS = "\t"} $4 == '$sign' {print (0.5+($4*$5/2)) "\t" $3}') | awk '{print "'$pert'\t'$sign'\t" $0}'
done
done

