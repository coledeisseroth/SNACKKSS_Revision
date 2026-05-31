#!/bin/bash

rm -rf gene drug

#Train and evaluate one model
for pert in gene drug; do
mkdir -p $pert
mkdir $pert/CREEDS_models $pert/CREEDS_checkpoint_predictions $pert/creedstrain_smctest_checkpoint_predictions
for split in $(ls ../import/control/$pert/CREEDS_training_datasets | cut -d. -f1 | sort -u); do
python3 ../src/text_classification_finetune.py ../import/control/$pert/CREEDS_training_datasets/$split.json ${pert}/CREEDS_models/$split microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext 1634
for ckpt in $(ls $pert/CREEDS_models/$split | grep checkpoint); do
python3 ../src/text_classification_predict.py $pert/CREEDS_models/$split/$ckpt ../import/control/$pert/CREEDS_testing_datasets/$split.txt > $pert/CREEDS_checkpoint_predictions/${split}_$ckpt.txt
python3 ../src/text_classification_predict.py $pert/CREEDS_models/$split/$ckpt ../import/control/$pert/SNACKKSS_MC_testing_datasets/$split.txt > $pert/creedstrain_smctest_checkpoint_predictions/${split}_$ckpt.txt
done
done
done
rm -r */CREEDS_models/*/checkpoint-*

for pert in gene drug; do
mkdir $pert/SNACKKSS_MC_models $pert/SNACKKSS_MC_checkpoint_predictions
for split in $(ls ../import/control/$pert/SNACKKSS_MC_training_datasets | cut -d. -f1 | sort -u); do
python3 ../src/text_classification_finetune.py ../import/control/$pert/SNACKKSS_MC_training_datasets/$split.json $pert/SNACKKSS_MC_models/$split $pert/CREEDS_models/$split 1244
for ckpt in $(ls $pert/SNACKKSS_MC_models/$split | grep checkpoint); do
python3 ../src/text_classification_predict.py $pert/SNACKKSS_MC_models/$split/$ckpt ../import/control/$pert/SNACKKSS_MC_testing_datasets/$split.txt > $pert/SNACKKSS_MC_checkpoint_predictions/${split}_$ckpt.txt
done
done
done
rm -r */SNACKKSS_MC_models/*/checkpoint-*


#Performance while training on CREEDS
for pert in gene drug; do
for db in CREEDS SNACKKSS_MC creedstrain_smctest; do
testdb=$(echo $db | sed 's/creedstrain_smctest/SNACKKSS_MC/g')
for ckpt in $(ls $pert/${db}_checkpoint_predictions/ | cut -d- -f2- | cut -d. -f1 | sort -gu); do
for split in $(ls $pert/${db}_checkpoint_predictions | cut -d_ -f1 | sort -u); do
latest=$(ls $pert/${db}_checkpoint_predictions/ | grep ${split}_ | cut -d'-' -f2 | cut -d. -f1 | awk '$1 <= '$ckpt | sort -gr | head -1)
cat $pert/${db}_checkpoint_predictions/${split}_checkpoint-$latest.txt | cut -f-2 | awk 'BEGIN {FS = "_"} {print $2 "\t" $0}' | sort -k1,1 | join -t$'\t' - <(cat $pert/${db}_checkpoint_predictions/${split}_checkpoint-$latest.txt | cut -d_ -f2 | sort | uniq -c | awk '{print $2 "\t" 1 / $1}' | sort -u | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "_" $4 "\t" $7 "\t" $8}' | sort -k1,1 | join -t$'\t' <(cat ../import/control/$pert/${testdb}_testing_datasets/$split.txt | cut -f1 | awk 'BEGIN {FS = "_"} {print $1 "_" $2 "_" $3 "\t" $5}' | sort -k1,1) - | sed 's/POSITIVE/1/g' | sed 's/NEGATIVE/0/g' | awk 'BEGIN {FS = "\t"} {if($3 == 1){$3 = 1} else{$3 = 0} print $2 "\t" $3 "\t" $4}'
done | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($1 == 1 && $2 == 1){tp += $3} else if($1 == 1){fn += $3} else if($2 == 1){fp += $3}} END {print "'$db'\t'$ckpt'\t" tp "\t" fp "\t" fn}'
done
done | awk 'BEGIN {print "Predictions\tCheckpoint\tTrue positives\tFalse positives\tFalse negatives"} {print}' > $pert/saturation.txt
done


