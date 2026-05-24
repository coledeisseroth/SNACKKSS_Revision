#!/bin/bash

#Train on CREEDS
for pert in gene drug; do
mkdir $pert/CREEDS_models $pert/CREEDS_checkpoint_predictions $pert/creedstrain_smctest_checkpoint_predictions
for split in $(ls ../import/target/$pert/CREEDS_training_datasets | cut -d. -f1 | sort -u); do
python3 ../src/target_finetune.py ../import/target/$pert/CREEDS_training_datasets/$split.json $pert/CREEDS_models/$split microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext 15
for ckpt in $(ls $pert/CREEDS_models/$split | grep checkpoint); do
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/CREEDS_testing_datasets/$split.txt | cut -f1,3 | sort -u) $pert/CREEDS_models/$split/$ckpt) > $pert/CREEDS_checkpoint_predictions/${split}_$ckpt.txt
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/SNACKKSS_MC_testing_datasets/$split.txt | cut -f1,3 | sort -u) $pert/CREEDS_models/$split/$ckpt) > $pert/creedstrain_smctest_checkpoint_predictions/${split}_$ckpt.txt
done
done
done
rm -r */CREEDS_models/*/checkpoint-*

#Train on SNACKKSS-MC
for pert in gene drug; do
mkdir $pert/SNACKKSS_MC_models $pert/SNACKKSS_MC_checkpoint_predictions
for split in $(ls $pert/CREEDS_models/ | sort -u); do
python3 ../src/target_finetune.py ../import/target/$pert/SNACKKSS_MC_training_datasets/$split.json $pert/SNACKKSS_MC_models/$split $pert/CREEDS_models/$split 17
for ckpt in $(ls $pert/SNACKKSS_MC_models/$split | grep checkpoint); do
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/SNACKKSS_MC_testing_datasets/$split.txt | cut -f1,3 | sort -u) ${pert}/SNACKKSS_MC_models/$split/$ckpt) > $pert/SNACKKSS_MC_checkpoint_predictions/${split}_$ckpt.txt
done
done
done
rm -r */SNACKKSS_MC_models/*/checkpoint-*

#Do a formal comparison.
#Performance while training on CREEDS
for pert in gene drug; do
if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then label=GENE
else label=CHEM; fi
for db in CREEDS creedstrain_smctest SNACKKSS_MC; do
testdb=$(echo $db | sed 's/creedstrain_smctest/SNACKKSS_MC/g')
for ckpt in $(ls $pert/${db}_checkpoint_predictions/ | cut -d- -f2- | cut -d. -f1 | sort -gu); do
for split in $(ls $pert/${db}_checkpoint_predictions/ | cut -d_ -f1 | sort -u); do
latest=$(ls $pert/${db}_checkpoint_predictions/ | grep ${split}_ | cut -d'-' -f2 | cut -d. -f1 | awk '$1 <= '$ckpt | sort -gr | head -1)
cat ../import/target/$pert/${testdb}_testing_datasets/$split.txt | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u | comm - <(cat $pert/${db}_checkpoint_predictions/${split}_checkpoint-${latest}.txt | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u)
done | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($2 == ""){tp++} else if($1 == ""){fp++} else{fn++}} END {precision = tp / (tp + fp + 1); recall = tp / (tp + fn + 1); print tp "\t" fp "\t" fn}' | awk '{print "'$db'\t'$ckpt'\t" $0}'
done
done | awk 'BEGIN {print "Predictions\tCheckpoint\tTrue positives\tFalse positives\tFalse negatives"} {print}' > $pert/saturation.txt
done



