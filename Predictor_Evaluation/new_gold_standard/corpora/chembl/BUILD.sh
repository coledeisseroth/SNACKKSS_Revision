for file in chembl_36.0_moa.ttl.gz chembl_36.0_targetrel.ttl.gz chembl_36.0_activity.ttl.gz chembl_36.0_target.ttl.gz chembl_36.0_molecule.ttl.gz; do
if [ $(ls | grep $file | wc -l) -lt 1 ]; then wget https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBL-RDF/latest/$file; fi
done

zcat chembl_36.0_moa.ttl.gz | grep chembl_molecule: | grep hasMechanism | cut -d: -f2- | sed 's/ .*chembl_moa:/\t/g' | cut -d' ' -f1 | sort -u > molecule_to_mechanism.txt

zcat chembl_36.0_moa.ttl.gz | grep chembl_target | grep isTargetForMechanism | cut -d: -f2- | sed 's/ .*chembl_moa:/\t/g' | cut -d' ' -f1 | sort -u > target_to_mechanism.txt

zcat chembl_36.0_moa.ttl.gz | grep 'mechanismDescription\|mechanismActionType' | paste -sd$'\t' | sed 's/chembl_moa:/\n/g' | sed 's/ cco:mechanismDescription /\t/g' | sed 's/ ;\t\tcco:mechanismActionType /\t/g' | sed 's/" \./"/g' | sed 's/"//g' | sort -u > mechanism_description.txt

cat molecule_to_mechanism.txt | sort -k2,2 | join -t$'\t' -1 2 -2 2 - <(cat target_to_mechanism.txt | sort -k2,2) | join -t$'\t' -1 1 -2 2 - <(cat mechanism_description.txt | cut -f1,3 | sort -t$'\t' -k2,2 | join -t$'\t' -1 2 -2 1 - <(cat raw/actiontype_direction.txt | awk 'BEGIN {FS = "\t"} $2 != 0' | sort -t$'\t' -k1,1) | sort -t$'\t' -k2,2) | cut -f2- | sort -u > mod_target_action_direction.txt

zcat chembl_36.0_molecule.ttl.gz | grep 'chembl_molecule:\|rdfs:label' | awk 'BEGIN {FS = "\t"; cur = ""; names = ""} {if(substr($0, 0, 16) == "chembl_molecule:"){if(cur != "" && names != ""){print cur "\t" names} names = ""; cur = $1; gsub("chembl_molecule:", "", cur); gsub(" .*", "", cur)} else if(index($2, "rdfs:label") > 0){gsub("rdfs:label \"", "", $2); gsub("\".*", "", $2); if(index($2, cur) == 0){names = names "\t" $2}}} END {if(cur != "" && names != ""){print cur "\t" names}}' | sed 's/\t\t/\t/g' > drug_names.txt

zcat chembl_36.0_target.ttl.gz | awk 'BEGIN {FS = "\t"; cur = ""; species = ""; names = ""} {if(index($1, "chembl_target:") > 0){if(cur != "" && names != "" && species == "Homo sapiens"){print cur "\t" names} names = ""; species = ""; cur = $1; gsub("chembl_target:", "", cur); gsub(" .*", "", cur)} else if(index($2, "cco:organismName") > 0){species = $2; gsub("cco:organismName \"", "", species); gsub("\".*", "", species)} else if(index($2, "rdfs:label") > 0){gsub("rdfs:label \"", "", $2); gsub("\".*", "", $2); names = names "\t" $2}} END {if(cur != "" && names != "" && species == "Homo sapiens"){print cur "\t" names}}' | sed 's/\t\t/\t/g' | sort -u > protein_names.txt 

python3 ../../../src/clean_column.py <(python3 ../../../src/clean_column.py <(cat mod_target_action_direction.txt | cut -f1,2,4 | sort -t$'\t' -k1,1 | join -t$'\t' - <(cat drug_names.txt | sort -t$'\t' -k1,1) | cut -f2- | sort -t$'\t' -k1,1 | join -t$'\t' - <(cat protein_names.txt | sort -t$'\t' -k1,1) | cut -f2- | awk 'BEGIN {FS = "\t"} {print $2 "\t" $3 "\t" $1}') 1) 0 | sort -u > named_interactions_clean.txt

for split in $(ls ../../../import/pubchem_split_100K_clean/); do
python3 ../../../src/resolve_synonyms.py <(cat named_interactions_clean.txt | sort -t$'\t' -k2,2 | join -t$'\t' -1 2 -2 2 - <(cat ../../../import/entrez_human_aliases_tall_clean.txt | sort -t$'\t' -k2,2) | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "\t" $1 "\t" $4}' | sort -t$'\t' -k1,1) | cut -f1,3 | sed 's/_/\t/g' | sort -t$'\t' -k1,1 | join -t$'\t' -1 1 -2 2 - <(cat ../../../import/pubchem_split_100K_clean/$split | sort -t$'\t' -k2,2) | awk 'BEGIN {FS = "\t"} {print $3 "_" $2 "\t" $1 "\t" $4}'
done > options.txt
python3 ../../../src/resolve_synonyms.py <(cat options.txt | sort -u | sort -t$'\t' -k1,1) | cut -f1,3 | sed 's/_/\t/g' | awk '{print $3 "\t" $1 "\t" $2}' | sort -u > mod_to_entrez.txt
rm options.txt

cat mod_to_entrez.txt | cut -f-2 | sort -u | comm -23 - <(cat ../../../import/dgidb.txt | cut -f-2 | sort -u) | sed 's/\t/_/g' | sort -u | join -t$'\t' - <(cat mod_to_entrez.txt | awk '{print $1 "_" $2 "\t" $3}' | sort -k1,1) | sed 's/_/\t/g' | sort -u > nodgidb_mod_to_entrez.txt

