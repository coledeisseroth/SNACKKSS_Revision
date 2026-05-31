cat ../../../import/drugrepurposinghub/repo-drug-annotation-20200324.txt | cut -f1,3,4 | awk 'BEGIN {FS = "\t"} $3 != "" && $2 != ""' | sort -u | awk 'BEGIN {FS = "\t"} {gsub(" \\| ", "\t", $3); print $1 "\t" $2 "\t" $3}' | awk 'BEGIN {FS = "\t"} {for(i=3; i <= NF; i++){print $1 "\t" $2 "\t" $i}}' | sort -u > interactions_tall.txt

for i in 1 -1; do cat interactions_tall.txt | grep -wf <(cat raw/action_direction.txt | awk '$2 == '$i | cut -f1 | sort -u) | cut -f1,3 | sort -u | awk '{print $0 "\t'$i'"}'; done | sort -u > named_interactions.txt

python3 ../../../src/clean_column.py <(python3 ../../../src/clean_column.py named_interactions.txt 0) 1 | sort -u > named_interactions_clean.txt

for split in $(ls ../../../import/pubchem_split_100K_clean/); do
python3 ../../../src/resolve_synonyms.py <(cat named_interactions_clean.txt | sort -t$'\t' -k2,2 | join -t$'\t' -1 2 -2 2 - <(cat ../../../import/entrez_human_aliases_tall_clean.txt | sort -t$'\t' -k2,2) | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "\t" $1 "\t" $4}' | sort -t$'\t' -k1,1) | cut -f1,3 | sed 's/_/\t/g' | sort -t$'\t' -k1,1 | join -t$'\t' -1 1 -2 2 - <(cat ../../../import/pubchem_split_100K_clean/$split | sort -t$'\t' -k2,2) | awk 'BEGIN {FS = "\t"} {print $3 "_" $2 "\t" $1 "\t" $4}'
done > options.txt
python3 ../../../src/resolve_synonyms.py <(cat options.txt | sort -u | sort -t$'\t' -k1,1) | cut -f1,3 | sed 's/_/\t/g' | awk '{print $3 "\t" $1 "\t" $2}' | sort -u > mod_to_entrez.txt
rm options.txt

#Exclude anything present in DGIdb
cat mod_to_entrez.txt | cut -f-2 | sort -u | comm -23 - <(cat ../../../import/dgidb.txt | cut -f-2 | sort -u) | sed 's/\t/_/g' | sort -u | join -t$'\t' - <(cat mod_to_entrez.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | sort -u > nodgidb_mod_to_entrez.txt

