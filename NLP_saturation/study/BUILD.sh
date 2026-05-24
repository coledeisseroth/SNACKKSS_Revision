#!/bin/bash

rm -rf gene drug

#Train on CREEDS
for pert in gene drug; do
mkdir -p $pert
hfmodel='microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext'
mkdir $pert/CREEDS_models $pert/CREEDS_checkpoint_predictions $pert/creedstrain_smctest_checkpoint_predictions
for split in $(ls ../import/study/$pert/CREEDS_training_datasets | cut -d. -f1 | sort -u); do
python3 ../src/text_classification_finetune.py ../import/study/$pert/CREDS_training_datasets/$split.json $pert/CREEDS_models/$split $hfmodel 392
for ckpt in $(ls $pert/CREEDS_models/$split | grep checkpoint); do
python3 ../src/text_classification_predict.py $pert/CREEDS_models/$split/$ckpt ../import/study/$pert/CREEDS_testing_datasets/$split.txt > $pert/CREEDS_checkpoint_predictions/${split}_$ckpt.txt
python3 ../src/text_classification_predict.py $pert/CREEDS_models/$split/$ckpt ../import/study/$pert/SNACKKSS_MC_testing_datasets/$split.txt > $pert/creedstrain_smctest_checkpoint_predictions/${split}_$ckpt.txt
done
done
done
rm -r */CREEDS_models/*/checkpoint-*

#Train on SNACKKSS_MC
for pert in gene drug; do
mkdir $pert/SNACKKSS_MC_models $pert/SNACKKSS_MC_checkpoint_predictions
for split in $(ls $pert/CREEDS_models | sort -u); do
python3 ../src/text_classification_finetune.py ../import/study/$pert/SNACKKSS_MC_training_datasets/$split.json $pert/SNACKKSS_MC_models/$split $pert/CREEDS_models/$split 144
for ckpt in $(ls $pert/SNACKKSS_MC_models/$split | grep checkpoint); do
python3 ../src/text_classification_predict.py $pert/SNACKKSS_MC_models/$split/$ckpt ../import/study/$pert/SNACKKSS_MC_testing_datasets/$split.txt > $pert/SNACKKSS_MC_checkpoint_predictions/${split}_$ckpt.txt
done
done
done
rm -r */SNACKKSS_MC_models/*/checkpoint-*

#Do a formal comparison.
for pert in gene drug; do
for db in CREEDS SNACKKSS_MC creedstrain_smctest; do
for ckpt in $(ls $pert/${db}_checkpoint_predictions/ | cut -d_ -f2- | cut -d- -f2 | cut -d. -f1 | sort -gu); do
for split in $(ls $pert/${db}_checkpoint_predictions | cut -d_ -f1 | sort -u); do
latest=$(ls $pert/${db}_checkpoint_predictions/ | grep ${split}_ | cut -d'-' -f2 | cut -d. -f1 | awk '$1 <= '$ckpt | sort -gr | head -1)
cat $pert/${db}_checkpoint_predictions/${split}_checkpoint-$latest.txt | awk 'BEGIN {FS = "\t"} {gsub("_", "\t", $1); print $1 "\t" $2 "\t" $3}' | cut -f1,3,4 | sort -k3,3r | sort -k1,1 -u | cut -f2-
done | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($1 == 0 && $2 == "POSITIVE"){fp++} else if($1 != 0 && $2 == "NEGATIVE"){fn++} else if($1 != 0 && $2 == "POSITIVE"){tp++}} END {print "'$db'\t'$ckpt'\t" tp "\t" fp "\t" fn}'
done
done | awk 'BEGIN {print "Predictions\tCheckpoint\tTrue positives\tFalse positives\tFalse negatives"} {print}' > $pert/saturation.txt
done

