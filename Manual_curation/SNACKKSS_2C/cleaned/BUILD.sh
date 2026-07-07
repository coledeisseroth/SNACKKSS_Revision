#Clean the annotations
python3 ../src/2c_gsm_fetch.py <(cat ../SNACKKSS_2C.txt | sed 's/;\t/\t/g' | awk 'NF > 1') <(cat ../../import/sample_info.txt | awk '{print $2 "\t" $0}' | sort -k1,1 | join -t$'\t' - <(cat ../../import/sample_library_strategy.txt | awk '$2 == "RNA-Seq"' | cut -f1 | sort -u) | cut -f2- | sort -u) | awk 'index($1 $2 $3 $4, "#") == 0' > new_perturbations_cleaned.txt

python3 ../src/mc_gsm_fetch.py <(cat ../../import/corrected_curated_dataset.txt | sed 's/;\t/\t/g' | awk 'NF > 1') <(cat ../../import/sample_info.txt | awk '{print $2 "\t" $0}' | sort -k1,1 | join -t$'\t' - <(cat ../../import/sample_library_strategy.txt | awk '$2 == "RNA-Seq"' | cut -f1 | sort -u) | cut -f2- | sort -u) | awk 'index($1 $2 $3 $4, "#") == 0' > old_perturbations_cleaned.txt

