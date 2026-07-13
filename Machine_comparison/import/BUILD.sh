#This script assumes that you have a completed SNACKKSS-NLP run from both cpus, each contained within a directory titled "64core" or "40core"
for task in study sample target control; do
mkdir $task
for core in 64core 40core
mkdir $task/$core
for pert in gene drug; do
for model in distilbert biobert biomedbert; do
for preds in 1.1 1.2 12.1 12.2 21.1 21.2 2.1 2.2; do
for split in $(ls ./$core/SNACKKSS_NLP/study/$pert/combo/$model/predictions$preds/); do
cp ./$core/SNACKKSS_NLP/study/$pert/combo/$model/predictions$preds/$split > $task/$core/${pert}_${model}_${preds}_${split}
done
done
done
done
done

cp ./64core/SNACKKSS_NLP/metadata/SNACKKSS_MC/split/ .
cp ./64core/SNACKKSS_NLP/metadata/SNACKKSS_MC/sample_info.txt .

mkdir position_labels
for pert in gene drug; do
for model in distilbert biobert biomedbert; do
cp -r ./64core/SNACKKSS_NLP/target/$pert/SNACKKSS_MC/$model/position_labels_test/ position_labels/${pert}_${model}
done
done

mkdir labeled_pairs
for pert in gene drug; do 
cp ./64core/SNACKKSS_NLP/control/$pert/SNACKKSS_MC/labeled_pairs.txt labeled_pairs/$pert.txt
done

