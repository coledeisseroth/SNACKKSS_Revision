#!/bin/bash

rm -rf gene drug

#Train on CREEDS
for pert in gene drug; do
mkdir -p $pert
hfmodel='microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext'
mkdir $pert/CREEDS_models $pert/CREEDS_checkpoint_predictions $pert/creedstrain_smctest_checkpoint_predictions
mkdir $pert/CREEDS_2epoch_batch16_models $pert/CREEDS_2epoch_batch16_predictions $pert/creedstrain_smctest_2epoch_batch16_predictions
for split in $(ls ../import/control/$pert/CREEDS_training_datasets | cut -d. -f1 | sort -u); do
python3 ../src/text_classification_finetune.py ../import/control/$pert/CREEDS_training_datasets/$split.json $pert/CREEDS_2epoch_batch16_models/$split $hfmodel 100000 2 16
python3 ../src/text_classification_predict.py $pert/CREEDS_2epoch_batch16_models/$split ../import/control/$pert/CREEDS_testing_datasets/$split.txt > $pert/CREEDS_2epoch_batch16_predictions/$split.txt
python3 ../src/text_classification_predict.py $pert/CREEDS_2epoch_batch16_models/$split ../import/control/$pert/SNACKKSS_MC_testing_datasets/$split.txt > $pert/creedstrain_smctest_2epoch_batch16_predictions/$split.txt
done
for batchsize in 1 16; do
mkdir -p $pert/CREEDS_models/$batchsize $pert/CREEDS_checkpoint_predictions/$batchsize $pert/creedstrain_smctest_checkpoint_predictions/$batchsize $pert/CREEDS_1epoch_predictions/$batchsize $pert/creedstrain_smctest_1epoch_predictions/$batchsize
for split in $(ls ../import/control/$pert/CREEDS_training_datasets | cut -d. -f1 | sort -u); do
steps=$(cat ../import/control/$pert/CREEDS_training_datasets/$split.json | wc -l | awk '{print int($1 / (11 * '$batchsize'))}')
python3 ../src/text_classification_finetune.py ../import/control/$pert/CREEDS_training_datasets/$split.json $pert/CREEDS_models/$batchsize/$split $hfmodel $steps 1 $batchsize
for ckpt in $(ls $pert/CREEDS_models/$batchsize/$split | grep checkpoint); do
python3 ../src/text_classification_predict.py $pert/CREEDS_models/$batchsize/$split/$ckpt ../import/control/$pert/CREEDS_testing_datasets/$split.txt > $pert/CREEDS_checkpoint_predictions/$batchsize/${split}_$ckpt.txt
python3 ../src/text_classification_predict.py $pert/CREEDS_models/$batchsize/$split/$ckpt ../import/control/$pert/SNACKKSS_MC_testing_datasets/$split.txt > $pert/creedstrain_smctest_checkpoint_predictions/$batchsize/${split}_$ckpt.txt
done
python3 ../src/text_classification_predict.py $pert/CREEDS_models/$batchsize/$split ../import/control/$pert/CREEDS_testing_datasets/$split.txt > $pert/CREEDS_1epoch_predictions/$batchsize/$split.txt
python3 ../src/text_classification_predict.py $pert/CREEDS_models/$batchsize/$split ../import/control/$pert/SNACKKSS_MC_testing_datasets/$split.txt > $pert/creedstrain_smctest_1epoch_predictions/$batchsize/$split.txt
done
done
done
rm -r */CREEDS_models/*/*/checkpoint-*

#Train on SNACKKSS_MC
for pert in gene drug; do
mkdir $pert/SNACKKSS_MC_models $pert/SNACKKSS_MC_checkpoint_predictions
mkdir $pert/SNACKKSS_MC_2epoch_batch16_models $pert/SNACKKSS_MC_2epoch_batch16_predictions
for split in $(ls $pert/CREEDS_2epoch_batch16_models | sort -u); do
python3 ../src/text_classification_finetune.py ../import/control/$pert/SNACKKSS_MC_training_datasets/$split.json $pert/SNACKKSS_MC_2epoch_batch16_models/$split $pert/CREEDS_2epoch_batch16_models/$split 100000 2 16
python3 ../src/text_classification_predict.py $pert/SNACKKSS_MC_2epoch_batch16_models/$split ../import/control/$pert/SNACKKSS_MC_testing_datasets/$split.txt > $pert/SNACKKSS_MC_2epoch_batch16_predictions/$split.txt
done
for batchsize in 1 16; do
creedsmodeldir=../import/control/$pert/CREEDS_models/
if [ $batchsize -eq 16 ]; then creedsmodeldir=$pert/CREEDS_2epoch_batch16_models/; fi
mkdir -p $pert/SNACKKSS_MC_models/$batchsize $pert/SNACKKSS_MC_checkpoint_predictions/$batchsize $pert/SNACKKSS_MC_1epoch_predictions/$batchsize
for split in $(ls $pert/CREEDS_models/$batchsize | sort -u); do
steps=$(cat ../import/control/$pert/SNACKKSS_MC_training_datasets/$split.json | wc -l | awk '{print int($1 / (11 * '$batchsize'))}')
python3 ../src/text_classification_finetune.py ../import/control/$pert/SNACKKSS_MC_training_datasets/$split.json $pert/SNACKKSS_MC_models/$batchsize/$split $creedsmodeldir/$split $steps 1 $batchsize
for ckpt in $(ls $pert/SNACKKSS_MC_models/$batchsize/$split | grep checkpoint); do
python3 ../src/text_classification_predict.py $pert/SNACKKSS_MC_models/$batchsize/$split/$ckpt ../import/control/$pert/SNACKKSS_MC_testing_datasets/$split.txt > $pert/SNACKKSS_MC_checkpoint_predictions/$batchsize/${split}_$ckpt.txt
done
python3 ../src/text_classification_predict.py $pert/SNACKKSS_MC_models/$batchsize/$split ../import/control/$pert/SNACKKSS_MC_testing_datasets/$split.txt > $pert/SNACKKSS_MC_1epoch_predictions/$batchsize/$split.txt
done
done
done
rm -r */SNACKKSS_MC_models/*/*/checkpoint-*

#Do a formal comparison.
for pert in gene drug; do
for db in CREEDS SNACKKSS_MC creedstrain_smctest; do
testdb=$(echo $db | sed 's/creedstrain_smctest/SNACKKSS_MC/g')
for batchsize in 1 16; do
cat $pert/${db}_1epoch_predictions/$batchsize/*.txt | cut -f-2 | awk 'BEGIN {FS = "_"} {print $2 "\t" $0}' | sort -k1,1 | join -t$'\t' - <(cat $pert/${db}_1epoch_predictions/$batchsize/*.txt | cut -d_ -f2 | sort | uniq -c | awk '{print $2 "\t" 1 / $1}' | sort -u | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "_" $4 "\t" $7 "\t" $8}' | sort -k1,1 | join -t$'\t' <(cat ../import/control/$pert/${testdb}_testing_datasets/*.txt | cut -f1 | awk 'BEGIN {FS = "_"} {print $1 "_" $2 "_" $3 "\t" $5}' | sort -k1,1) - | sed 's/POSITIVE/1/g' | sed 's/NEGATIVE/0/g' | awk 'BEGIN {FS = "\t"} {if($3 == 1){$3 = 1} else{$3 = 0} print $2 "\t" $3 "\t" $4}' | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($1 == 1 && $2 == 1){tp += $3} else if($1 == 1){fn += $3} else if($2 == 1){fp += $3}} END {print "'$pert'\t'$db'\t'$batchsize'\t1\t" tp "\t" fp "\t" fn}'
done
done
done | awk 'BEGIN {print "Perturbation type\tPredictions\tBatch Size\tProgress\tTrue positives\tFalse positives\tFalse negatives"} {print}' > 1epoch_performance.txt

for pert in gene drug; do
for db in CREEDS SNACKKSS_MC creedstrain_smctest; do
testdb=$(echo $db | sed 's/creedstrain_smctest/SNACKKSS_MC/g')
traindb=SNACKKSS_MC
if [ $(echo $db | grep SNACKKSS_MC | wc -l) -lt 1 ]; then traindb=CREEDS; fi
for batchsize in 1 16; do
for progress in $(seq 0 9 | awk '{print $1 / 10}'); do
for split in $(ls $pert/${db}_checkpoint_predictions/$batchsize/ | cut -d_ -f1 | sort -u); do
steps=$(cat ../import/control/$pert/${traindb}_training_datasets/$split.json | wc -l | awk '{print int($1 / '$batchsize')}')
soonest=$(ls $pert/${db}_checkpoint_predictions/$batchsize/ | grep $split | cut -d- -f2 | cut -d. -f1 | awk '$1 / '$steps' > '$progress'' | sort -gu | head -1)
cat $pert/${db}_checkpoint_predictions/$batchsize/${split}_checkpoint-$soonest.txt | cut -f-2 | awk 'BEGIN {FS = "_"} {print $2 "\t" $0}' | sort -k1,1 | join -t$'\t' - <(cat $pert/${db}_checkpoint_predictions/$batchsize/${split}_checkpoint-$soonest.txt | cut -d_ -f2 | sort | uniq -c | awk '{print $2 "\t" 1 / $1}' | sort -u | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "_" $4 "\t" $7 "\t" $8}' | sort -k1,1 | join -t$'\t' <(cat ../import/control/$pert/${testdb}_testing_datasets/$split.txt | cut -f1 | awk 'BEGIN {FS = "_"} {print $1 "_" $2 "_" $3 "\t" $5}' | sort -k1,1) - | sed 's/POSITIVE/1/g' | sed 's/NEGATIVE/0/g' | awk 'BEGIN {FS = "\t"} {if($3 == 1){$3 = 1} else{$3 = 0} print $2 "\t" $3 "\t" $4}'
done | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($1 == 1 && $2 == 1){tp += $3} else if($1 == 1){fn += $3} else if($2 == 1){fp += $3}} END {print "'$db'\t'$batchsize'\t'$progress'\t" tp "\t" fp "\t" fn}'
done
cat 1epoch_performance.txt | awk '$1 == "'$pert'" && $2 == "'$db'" && $3 == "'$batchsize'"' | cut -f2-
done
done | awk 'BEGIN {print "Predictions\tBatch size\tProgress\tTrue positives\tFalse positives\tFalse negatives"} {print}' > $pert/saturation.txt
done

for pert in gene drug; do
for db in CREEDS SNACKKSS_MC creedstrain_smctest; do
testdb=$(echo $db | sed 's/creedstrain_smctest/SNACKKSS_MC/g')
cat $pert/${db}_2epoch_batch16_predictions/*.txt | cut -f-2 | awk 'BEGIN {FS = "_"} {print $2 "\t" $0}' | sort -k1,1 | join -t$'\t' - <(cat $pert/${db}_2epoch_batch16_predictions/*.txt | cut -d_ -f2 | sort | uniq -c | awk '{print $2 "\t" 1 / $1}' | sort -u | sort -k1,1) | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "_" $4 "\t" $7 "\t" $8}' | sort -k1,1 | join -t$'\t' <(cat ../import/control/$pert/${testdb}_testing_datasets/*.txt | cut -f1 | awk 'BEGIN {FS = "_"} {print $1 "_" $2 "_" $3 "\t" $5}' | sort -k1,1) - | sed 's/POSITIVE/1/g' | sed 's/NEGATIVE/0/g' | awk 'BEGIN {FS = "\t"} {if($3 == 1){$3 = 1} else{$3 = 0} print $2 "\t" $3 "\t" $4}' | awk 'BEGIN {FS = "\t"; tp = 0; fp = 0; fn = 0} {if($1 == 1 && $2 == 1){tp += $3} else if($1 == 1){fn += $3} else if($2 == 1){fp += $3}} END {print "'$pert'\t'$db'\t16\t2\t" tp "\t" fp "\t" fn}'
done
done | awk 'BEGIN {print "Perturbation type\tPredictions\tBatch Size\tProgress\tTrue positives\tFalse positives\tFalse negatives"} {print}' > batch16_performance.txt


