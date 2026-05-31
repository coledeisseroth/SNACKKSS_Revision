if [ $(ls | grep trrust_rawdata.human.tsv | wc -l) -lt 1 ]; then wget https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv; fi

python3 ../../../src/resolve_synonyms.py <(python3 ../../../src/resolve_synonyms.py <(python3 ../../../src/clean_column.py <(python3 ../../../src/clean_column.py <(cat trrust_rawdata.human.tsv | cut -f-3 | awk 'BEGIN {FS = "\t"} {if($3 == "Activation"){print $1 "\t" $2 "\t1"} else if($3 == "Repression"){print $1 "\t" $2 "\t-1"}}') 1) 0 | sort -u | sort -k1,1 | join -t$'\t' -1 1 -2 2 - <(cat ../../../import/entrez_human_aliases_tall_clean.txt | sort -k2,2) | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "\t" $1 "\t" $4}' | sort -k1,1 -k2,2 -k3,3) | cut -f1,3 | sort -u | sed 's/_/\t/g' | sort -t$'\t' -k1,1 | join -t$'\t' -1 1 -2 2 - <(cat ../../../import/entrez_human_aliases_tall_clean.txt | sort -k2,2) | awk 'BEGIN {FS = "\t"} {print $3 "_" $2 "\t" $1 "\t" $4}' | sort -k1,1 -k2,2 -k3,3) | cut -f1,3 | sort -u | sed 's/_/\t/g' | awk 'BEGIN {FS = "\t"} {print $1 "\t" $3 "\t" $2}' | sort -u > mod_to_entrez.txt

#See how much it has in common with Reactome's FI dataset
cat ../../../import/reactome.txt | cut -f-2 | sort -u | comm -13 - <(cat mod_to_entrez.txt | cut -f-2 | sort -u) | sed 's/\t/_/g' | sort -u | join -t$'\t' - <(cat mod_to_entrez.txt | awk '{print $1 "_" $2 "\t" $0}' | sort -k1,1) | cut -f2- | sort -u > noreactome_mod_to_entrez.txt

