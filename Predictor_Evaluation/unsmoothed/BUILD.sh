#Plot accuracy against coverage for all predictors
(for file in $(ls ../import/SNACKKSS/output_stats | grep _accuracy.txt); do for col in $(seq $(cat ../import/SNACKKSS/output_stats/$file | head -1 | awk 'BEGIN {FS = "\t"}{print int(NF / 2)}') | awk '{print $1*2}'); do head -1 ../import/SNACKKSS/output_stats/$file | cut -f$col | awk '{print "'$file' "$0}' | sed 's/_accuracy.txt//g' | sed 's/ true//g' | sed 's/Gene/gene/g' | sed 's/Drug/drug/g' | sed 's/_/ /g'; done; done | paste -sd$'\t' | awk '{print "\t" $0}'
for acc in $(seq 0 999 | awk '{print $1 / 1000}'); do (echo $acc; for file in $(ls ../import/SNACKKSS/output_stats | grep _accuracy.txt); do for col in $(seq $(cat ../import/SNACKKSS/output_stats/$file | head -1 | awk 'BEGIN {FS = "\t"}{print int(NF / 2)}') | awk '{print $1*2}'); do cat ../import/SNACKKSS/output_stats/$file | cut -f$col,$(echo $col | awk '{print $1 + 1}') | awk 'BEGIN {print 0} {if($1 + $2 == 0){next} else if($1 / ($1+$2) > '$acc'){print}}' | cut -f1 | sort -gr | head -1; done; done) | paste -sd$'\t'; done) > accuracy_vs_coverage.txt

#Ablation test
mkdir ablation
genepredfiles="../import/PARMESAN/gene_loo_consensus.txt ../import/PARMESAN/gene_loo_hypotheses.txt ../import/PARMESAN/gene_loo_archs4_predictions.txt ../import/PubTator3/gene_loo_consensus.txt ../import/PubTator3/gene_loo_hypotheses.txt ../import/PubTator3/gene_loo_archs4_predictions.txt ../import/archs4/human_archs4_loo_preds.txt ../import/archs4/mouse_archs4_loo_preds.txt ../import/ConnectivityMap/loo_prediction_accuracy_estimates/gene.txt ../import/SNACKKSS/loo_prediction_accuracy_estimates/gene.txt"
drugpredfiles=$(echo $genepredfiles | sed 's/ /\n/g' | grep -v archs4_loo_preds | sed 's/gene/drug/g' | paste -sd' ')
for pert in gene drug; do
predfiles=$genepredfiles
if [ $(echo $pert | grep drug | wc -l) -gt 0 ]; then predfiles=$drugpredfiles; fi
for sign in pos neg; do
sgn='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then sgn='<'; fi
(if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then echo $'Mimimum accuracy\tNothing\tPARMESAN consensus\tPARMESAN hypotheses\tPARMESAN ARCHS4-linked hypotheses\tPubTator3 consensuses\tPubTator3 hypotheses\tPubTator3 ARCHS4-linked hypotheses\tHuman ARCHS4 coexpression\tMouse ARCHS4 coexpression\tCMap ARCHS4-linked hypotheses\tSNACKKSS ARCHS4-linked hypotheses'
else echo $'Mimimum accuracy\tNothing\tPARMESAN consensus\tPARMESAN hypotheses\tPARMESAN ARCHS4-linked hypotheses\tPubTator3 consensuses\tPubTator3 hypotheses\tPubTator3 ARCHS4-linked hypotheses\tCMap ARCHS4-linked hypotheses\tSNACKKSS ARCHS4-linked hypotheses'
fi
for i in $(seq 0 999 | awk '{print $1 / 1000}'); do echo -n $i$'\t'; for excluded in everything $predfiles; do cat $(echo $predfiles | sed 's/ /\n/g' | grep -vw $excluded | paste -sd' ') | sort -k3,3gr -k4,4gr | sort -k1,1 -k2,2 -u | awk '$4 '$sgn' 0' | sort -k3,3gr | awk 'BEGIN {FS = "\t"; cur = ""; p = 0; n = 0} {if(cur != "" && cur != $3){if(p+n == 0){acc = 0}else{acc = p/(p+n)} print cur "\t" p "\t" n "\t" acc} cur = $3; if($4 * $5 > 0){p++} else{n++}} END {if(p+n == 0){acc = 0}else{acc = p/(p+n)} print cur "\t" p "\t" n "\t" acc; print "1\t0\t0\t1"}' | awk '$4 > '$i | sort -k2,2gr | head -1 | cut -f2; done | paste -sd$'\t'
done) > ablation/${pert}_${sign}_ablation_test.txt &
done
done
#Ablation-esque LOO test using each predictor alone
genepredfiles="../import/PARMESAN/gene_loo_consensus.txt ../import/PARMESAN/gene_loo_hypotheses.txt ../import/PARMESAN/gene_loo_archs4_predictions.txt ../import/PubTator3/gene_loo_consensus.txt ../import/PubTator3/gene_loo_hypotheses.txt ../import/PubTator3/gene_loo_archs4_predictions.txt ../import/archs4/human_archs4_loo_preds.txt ../import/archs4/mouse_archs4_loo_preds.txt ../import/ConnectivityMap/loo_prediction_accuracy_estimates/gene.txt ../import/SNACKKSS/loo_prediction_accuracy_estimates/gene.txt"
drugpredfiles=$(echo $genepredfiles | sed 's/ /\n/g' | grep -v archs4_loo_preds | sed 's/gene/drug/g' | paste -sd' ')
for pert in gene drug; do
predfiles=$genepredfiles
if [ $(echo $pert | grep drug | wc -l) -gt 0 ]; then predfiles=$drugpredfiles; fi
for sign in pos neg; do
sgn='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then sgn='<'; fi
(if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then echo $'Mimimum accuracy\tPARMESAN consensus\tPARMESAN hypotheses\tPARMESAN ARCHS4-linked hypotheses\tPubTator3 consensuses\tPubTator3 hypotheses\tPubTator3 ARCHS4-linked hypotheses\tHuman ARCHS4 coexpression\tMouse ARCHS4 coexpression\tCMap ARCHS4-linked hypotheses\tSNACKKSS ARCHS4-linked hypotheses'
else echo $'Mimimum accuracy\tPARMESAN consensus\tPARMESAN hypotheses\tPARMESAN ARCHS4-linked hypotheses\tPubTator3 consensuses\tPubTator3 hypotheses\tPubTator3 ARCHS4-linked hypotheses\tCMap ARCHS4-linked hypotheses\tSNACKKSS ARCHS4-linked hypotheses'
fi
for i in $(seq 0 999 | awk '{print $1 / 1000}'); do echo -n $i$'\t'; for preds in $predfiles; do cat $preds | sort -k3,3gr -k4,4gr | sort -k1,1 -k2,2 -u | awk '$4 '$sgn' 0' | sort -k3,3gr | awk 'BEGIN {FS = "\t"; cur = ""; p = 0; n = 0} {if(cur != "" && cur != $3){if(p+n == 0){acc = 0}else{acc = p/(p+n)} print cur "\t" p "\t" n "\t" acc} cur = $3; if($4 * $5 > 0){p++} else{n++}} END {if(p+n == 0){acc = 0}else{acc = p/(p+n)} print cur "\t" p "\t" n "\t" acc; print "1\t0\t0\t1"}' | awk '$4 > '$i | sort -k2,2gr | head -1 | cut -f2; done | paste -sd$'\t'
done) > ablation/${pert}_${sign}_loo_test.txt &
done
done
#Ablation test with a focus on number of targets instead of number of relations
for i in $(seq 0 999 | awk '{print $1 / 1000}'); do echo -n $i$'\t'; for excluded in nothing '../import/SNACKKSS/loo_prediction_accuracy_estimates/drug.txt'; do
cat $(echo $drugpredfiles | sed 's/ /\n/g' | grep -vw $excluded | paste -sd' ') | sort -k3,3gr -k4,4gr | sort -k1,1 -k2,2 -u | awk '$4 < 0' | sort -k3,3gr | awk 'BEGIN {FS = "\t"; cur = ""; p = 0; n = 0} {if(p+n>0){print $2 "\t" $3 "\t" p/(p+n)} else{print $2 "\t" $3 "\t" 0} if($4 * $5 > 0){p++} else{n++}}' | sort -k2,2g -k3,3gr | awk 'BEGIN {FS = "\t"; cur = ""} {if(cur < $3){cur = $3} if(cur > '$i'){print $1 "\t" cur}}' | cut -f1 | sort -u | wc -l
done | paste -sd$'\t'
done > ablation/SA4_target_coverage_ablation_test.txt &
while [ $(jobs | grep Running | wc -l) -gt 0 ]; do jobs; sleep 1; done


#Get the exact accuracy caps after ablating each predictor
(echo $'Perturbation type\tDirection\tExcluded database\tAccuracy cap'
for pert in gene drug; do
genepredfiles="../import/PARMESAN/gene_loo_consensus.txt ../import/PARMESAN/gene_loo_hypotheses.txt ../import/PARMESAN/gene_loo_archs4_predictions.txt ../import/PubTator3/gene_loo_consensus.txt ../import/PubTator3/gene_loo_hypotheses.txt ../import/PubTator3/gene_loo_archs4_predictions.txt ../import/archs4/human_archs4_loo_preds.txt ../import/archs4/mouse_archs4_loo_preds.txt ../import/ConnectivityMap/loo_prediction_accuracy_estimates/gene.txt ../import/SNACKKSS/loo_prediction_accuracy_estimates/gene.txt"
drugpredfiles=$(echo $genepredfiles | sed 's/ /\n/g' | grep -v archs4_loo_preds | sed 's/gene/drug/g' | paste -sd' ')
predfiles=$genepredfiles
if [ $(echo $pert | grep drug | wc -l) -gt 0 ]; then predfiles=$drugpredfiles; fi
for sign in pos neg; do
sgn='>'
if [ $(echo $sign | grep neg | wc -l) -gt 0 ]; then sgn='<'; fi
for excluded in nothing $predfiles; do
cat $(echo $predfiles | sed 's/ /\n/g' | grep -vw $excluded | paste -sd' ') | sed 's/_/\t/g' | sort -k3,3gr -k4,4gr | sort -k1,1 -k2,2 -u | sort -k3,3gr | awk 'BEGIN {FS = "\t"; p = 0; n = 0; cur = ""; champ=0} $4 '$sgn' 0 {if(cur != "" && cur != $3 && p + n > 0){acc = p/(p+n); if(acc > champ){champ = acc}} cur = $3; if($4 * $5 > 0){p++} else{n++}} END {print "'$pert'\t'$sign'\t'$excluded'\t" champ}'
done
done
done) > ablation_accuracy_caps.txt


