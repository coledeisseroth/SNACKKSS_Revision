#!/bin/bash

#rm -rf gene drug

#Train on CREEDS
hfmodel=microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext
for pert in gene drug; do
mkdir $pert
mkdir $pert/CREEDS_models $pert/CREEDS_checkpoint_predictions $pert/creedstrain_smctest_checkpoint_predictions
mkdir $pert/CREEDS_2epoch_batch1_models $pert/CREEDS_2epoch_batch1_predictions $pert/creedstrain_smctest_2epoch_batch1_predictions
for split in $(ls ../import/target/$pert/CREEDS_training_datasets | cut -d. -f1 | sort -u); do
python3 ../src/target_finetune.py ../import/target/$pert/CREEDS_training_datasets/$split.json $pert/CREEDS_2epoch_batch1_models/$split $hfmodel 100000 2 1
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/CREEDS_testing_datasets/$split.txt | cut -f1,3 | sort -u) $pert/CREEDS_2epoch_batch1_models/$split) > $pert/CREEDS_2epoch_batch1_predictions/$split.txt
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/SNACKKSS_MC_testing_datasets/$split.txt | cut -f1,3 | sort -u) $pert/CREEDS_2epoch_batch1_models/$split) > $pert/creedstrain_smctest_2epoch_batch1_predictions/$split.txt
done
for batchsize in 1 16; do
mkdir -p $pert/CREEDS_models/$batchsize $pert/CREEDS_checkpoint_predictions/$batchsize $pert/creedstrain_smctest_checkpoint_predictions/$batchsize $pert/CREEDS_1epoch_predictions/$batchsize $pert/creedstrain_smctest_1epoch_predictions/$batchsize
for split in $(ls ../import/target/$pert/CREEDS_training_datasets | cut -d. -f1 | sort -u); do
steps=$(cat ../import/target/$pert/CREEDS_training_datasets/$split.json | wc -l | awk '{print int($1 / (11 * '$batchsize'))}')
python3 ../src/target_finetune.py ../import/target/$pert/CREEDS_training_datasets/$split.json $pert/CREEDS_models/$batchsize/$split $hfmodel $steps 1 $batchsize
for ckpt in $(ls $pert/CREEDS_models/$batchsize/$split | grep checkpoint); do
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/CREEDS_testing_datasets/$split.txt | cut -f1,3 | sort -u) $pert/CREEDS_models/$batchsize/$split/$ckpt) > $pert/CREEDS_checkpoint_predictions/$batchsize/${split}_$ckpt.txt
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/SNACKKSS_MC_testing_datasets/$split.txt | cut -f1,3 | sort -u) $pert/CREEDS_models/$batchsize/$split/$ckpt) > $pert/creedstrain_smctest_checkpoint_predictions/$batchsize/${split}_$ckpt.txt
done
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/CREEDS_testing_datasets/$split.txt | cut -f1,3 | sort -u) $pert/CREEDS_models/$batchsize/$split) > $pert/CREEDS_1epoch_predictions/$batchsize/$split.txt
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/SNACKKSS_MC_testing_datasets/$split.txt | cut -f1,3 | sort -u) $pert/CREEDS_models/$batchsize/$split) > $pert/creedstrain_smctest_1epoch_predictions/$batchsize/$split.txt
done
done
done
rm -r */CREEDS_models/*/*/checkpoint-*

#Train on SNACKKSS-MC
for pert in gene drug; do
mkdir $pert/SNACKKSS_MC_models $pert/SNACKKSS_MC_checkpoint_predictions
mkdir $pert/SNACKKSS_MC_2epoch_batch1_models $pert/SNACKKSS_MC_2epoch_batch1_predictions
for split in $(ls ../import/target/$pert/CREEDS_training_datasets | cut -d. -f1 | sort -u); do
python3 ../src/target_finetune.py ../import/target/$pert/SNACKKSS_MC_training_datasets/$split.json $pert/SNACKKSS_MC_2epoch_batch1_models/$split $pert/CREEDS_2epoch_batch1_models/$split 100000 2 1
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/SNACKKSS_MC_testing_datasets/$split.txt | cut -f1,3 | sort -u) $pert/SNACKKSS_MC_2epoch_batch1_models/$split) > $pert/SNACKKSS_MC_2epoch_batch1_predictions/$split.txt
done
for batchsize in 1 16; do
creedsmodeldir=../import/target/$pert/CREEDS_models/
if [ $batchsize -eq 16 ]; then creedsmodeldir=$pert/CREEDS_2epoch_batch1_models/; fi
mkdir -p $pert/SNACKKSS_MC_models/$batchsize $pert/SNACKKSS_MC_checkpoint_predictions/$batchsize $pert/SNACKKSS_MC_1epoch_predictions/$batchsize
for split in $(ls $pert/CREEDS_models/$batchsize | sort -u); do
steps=$(cat ../import/target/$pert/SNACKKSS_MC_training_datasets/$split.json | wc -l | awk '{print int($1 / (11 * '$batchsize'))}')
python3 ../src/target_finetune.py ../import/target/$pert/SNACKKSS_MC_training_datasets/$split.json $pert/SNACKKSS_MC_models/$batchsize/$split $creedsmodeldir/$split $steps 1 $batchsize
for ckpt in $(ls $pert/SNACKKSS_MC_models/$batchsize/$split | grep checkpoint); do
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/SNACKKSS_MC_testing_datasets/$split.txt | cut -f1,3 | sort -u) ${pert}/SNACKKSS_MC_models/$batchsize/$split/$ckpt) > $pert/SNACKKSS_MC_checkpoint_predictions/$batchsize/${split}_$ckpt.txt
done
python3 ../src/merge_entities.py <(python3 ../src/target_predict.py <(cat ../import/target/$pert/SNACKKSS_MC_testing_datasets/$split.txt | cut -f1,3 | sort -u) $pert/SNACKKSS_MC_models/$batchsize/$split) > $pert/SNACKKSS_MC_1epoch_predictions/$batchsize/$split.txt
done
done
done
rm -r */SNACKKSS_MC_models/*/*/checkpoint-*

#Do a formal comparison.
for pert in gene drug; do
for db in CREEDS SNACKKSS_MC creedstrain_smctest; do
testdb=$(echo $db | sed 's/creedstrain_smctest/SNACKKSS_MC/g')
if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then label=GENE
else label=CHEM; fi
for batchsize in 1 16; do
cat ../import/target/$pert/${testdb}_testing_datasets/*.txt | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u | comm - <(cat $pert/${db}_1epoch_predictions/$batchsize/*.txt | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u) | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($2 == ""){tp++} else if($1 == ""){fp++} else{fn++}} END {precision = tp / (tp + fp + 1); recall = tp / (tp + fn + 1); print tp "\t" fp "\t" fn}' | awk '{print "'$pert'\t'$db'\t'$batchsize'\t" $0}'
done
done
done | awk 'BEGIN {print "Perturbation type\tTrial\tBatch size\tTrue positives\tFalse positives\tFalse negatives"} {print}' > 1epoch_performance.txt

for pert in gene drug; do
if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then label=GENE
else label=CHEM; fi
for db in CREEDS creedstrain_smctest SNACKKSS_MC; do
traindb=SNACKKSS_MC
if [ $(echo $db | grep SNACKKSS_MC | wc -l) -lt 1 ]; then traindb=CREEDS; fi
testdb=$(echo $db | sed 's/creedstrain_smctest/SNACKKSS_MC/g')
for batchsize in 1 16; do
for progress in $(seq 0 9 | awk '{print $1 / 10}'); do
for split in $(ls $pert/${db}_checkpoint_predictions/$batchsize | cut -d_ -f1 | sort -u); do
steps=$(cat ../import/target/$pert/${traindb}_training_datasets/$split.json | wc -l | awk '{print int($1 / '$batchsize')}')
soonest=$(ls $pert/${db}_checkpoint_predictions/$batchsize/ | grep $split | cut -d- -f2 | cut -d. -f1 | awk '$1 / '$steps' > '$progress'' | sort -gu | head -1)
cat ../import/target/$pert/${testdb}_testing_datasets/$split.txt | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u | comm - <(cat $pert/${db}_checkpoint_predictions/$batchsize/${split}_checkpoint-${soonest}.txt | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u)
done | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($2 == ""){tp++} else if($1 == ""){fp++} else{fn++}} END {precision = tp / (tp + fp + 1); recall = tp / (tp + fn + 1); print tp "\t" fp "\t" fn}' | awk '{print "'$db'\t'$batchsize'\t'$progress'\t" $0}'
done
cat 1epoch_performance.txt | awk '$1 == "'$pert'" && $2 == "'$db'" && $3 == "'$batchsize'"' | cut -f2-
done
done | awk 'BEGIN {print "Predictions\tBatch size\tProgress\tTrue positives\tFalse positives\tFalse negatives"} {print}' > $pert/saturation.txt
done

for pert in gene drug; do
for db in CREEDS SNACKKSS_MC creedstrain_smctest; do
testdb=$(echo $db | sed 's/creedstrain_smctest/SNACKKSS_MC/g')
if [ $(echo $pert | grep gene | wc -l) -gt 0 ]; then label=GENE
else label=CHEM; fi
cat ../import/target/$pert/${testdb}_testing_datasets/*.txt | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u | comm - <(cat $pert/${db}_2epoch_batch1_predictions/*.txt | cut -f-2 | sed 's/;/\t/g' | awk 'BEGIN {FS = "\t"} {for(i = 2; i <= NF; i++) {print $1 "\t" $i}}' | grep ":"$label | cut -d: -f1 | awk '$2 != ""' | sort -u) | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($2 == ""){tp++} else if($1 == ""){fp++} else{fn++}} END {precision = tp / (tp + fp + 1); recall = tp / (tp + fn + 1); print tp "\t" fp "\t" fn}' | awk '{print "'$pert'\t'$db'\t" $0}'
done
done | awk 'BEGIN {print "Perturbation type\tTrial\tTrue positives\tFalse positives\tFalse negatives"} {print}' > batch1_performance.txt


