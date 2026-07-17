#Download page here: https://www.guidetopharmacology.org/download.jsp#database
for file in approved_drug_detailed_interactions.tsv interactions.tsv; do
if [ $(ls | grep -w $file | wc -l) -lt 1 ]; then wget https://www.guidetopharmacology.org/DATA/$file; fi
done

#Make sure nothing in DGIdb is here
wget https://dgidb.org/data/2024-Dec/interactions.tsv -O dgidb_interactions.tsv
#Find any IDs in DGIdb that it claims were in GTP, but aren't currently
cat dgidb_interactions.tsv | awk 'BEGIN {FS = "\t"} $10 == "GuideToPharmacology"' | cut -f4 | cut -d: -f2 | sort -u | comm -23 - <(cat interactions.tsv | cut -f15 | sed 's/"//g' | sort -u)
#There were two renamed interactions. We had to annotate those manually
#Find any RELATIONS in DGIdb that it claims were in GTP, but aren't currently
comm -23 <(cat dgidb_interactions.tsv | awk 'BEGIN {FS = "\t"} $10 == "GuideToPharmacology"' | cut -f1,4 | sed 's/NCBIGENE://g' | sed 's/IUPHAR.LIGAND://g' | awk '{print $2 "\t" $1}' | sort -u) <(cat interactions.tsv | cut -f2,15 | sed 's/"//g' | awk 'BEGIN {FS = "\t"} NR > 2 && NF > 1 && $2 != "" && $1 != "" {print $2 "\t" $1}' | sort -k1,1 | join -t$'\t' -1 1 -2 2 -a 1 -o auto -e '.' - <(cat raw/renamed_ligands.txt | sort -k2,2) | awk 'BEGIN {FS = "\t"} {if($3 == "."){print $1 "\t" $2} else{print $3 "\t" $2}}' | sort -u)
#These were the relations that came up:
#12689	2388: Ligand renamed
#2092	360: Ligand renamed
#4942	1806: Relation no longer present (GM-CSF --> DPYD), and it was never directed to begin with
#8407	2446: Relation no longer present (Etrolizumab --> FRA11H), and it was never directed to begin with

(head -2 interactions.tsv | tail -1; cat dgidb_interactions.tsv | awk 'BEGIN {FS = "\t"} $10 == "GuideToPharmacology"' | cut -f1,4 | sed 's/NCBIGENE://g' | sed 's/IUPHAR.LIGAND://g' | awk '{print $2 "\t" $1}' | sort -u | sort -k1,1 | join -t$'\t' -a 1 -o auto -e '.' - <(cat raw/renamed_ligands.txt | sort -k1,1) | awk 'BEGIN {FS = "\t"} {if($3 == "."){print $1 "\t" $2} else{print $3 "\t" $2}}' | sort -u | comm -13 - <(cat interactions.tsv | cut -f2,15 | sed 's/"//g' | awk 'BEGIN {FS = "\t"} NR > 2 && NF > 1 && $2 != "" && $1 != "" {print $2 "\t" $1}' | sort -u) | sed 's/\t/_/g' | sort -u | join -t$'\t' - <(cat interactions.tsv | awk 'BEGIN {FS = "\t"} {x = $15 "_" $2; gsub("\"", "", x); print x "\t" $0}' | sort -k1,1) | cut -f2-)  > interactions_minus_dgidb_relations.tsv

cat interactions_minus_dgidb_relations.tsv | cut -f2,21,23,24 | sed 's/"//g' | awk 'BEGIN {FS = "\t"} NR > 2 && $1 != "" && $2 != "" && $3 != ""' | sort -u > nogtpdgidb_entrez_pcid_actions.txt

python3 ../../../src/clean_column.py <(cat interactions_minus_dgidb_relations.tsv | cut -f2,16,23,24 | sed 's/"//g' | awk 'BEGIN {FS = "\t"} NR > 2 && $1 != "" && $2 != "" && $3 != ""') 1 | sort -u > nogtpdgidb_entrez_drugname_actions.txt

for file in $(ls ../../../import/pubchem_split_100K_clean/); do cat nogtpdgidb_entrez_drugname_actions.txt | sort -t$'\t' -k2,2 | join -t$'\t' -1 2 -2 2 - <(cat ../../../import/pubchem_split_100K_clean/$file | sort -t$'\t' -k2,2) | awk 'BEGIN {FS = "\t"} {print $2 "_" $3 "_" $4 "\t" $1 "\t" $5}' | sort -u; done > options.txt

python3 ../../../src/resolve_synonyms.py <(cat options.txt | sort -t$'\t' -k1,1 -k2,2 -k3,3) | uniq | sort -u | awk 'BEGIN {FS = "\t"} {gsub("_", "\t", $1); print $3 "\t" $1}' > pcid_entrez_actions.txt
(cat pcid_entrez_actions.txt | cut -f-3 | sort -t$'\t' -k3,3 | join -t$'\t' -1 3 -2 1 - <(cat raw/type_direction.txt | sort -k1,1) | cut -f2-
cat pcid_entrez_actions.txt | cut -f1,2,4 | sort -t$'\t' -k3,3 | join -t$'\t' -1 3 -2 1 - <(cat raw/type_direction.txt | cut -f1 | sort -u) | cut -f2- | sort -u | comm -13 - <(cat pcid_entrez_actions.txt | cut -f-2 | sort -u) | sed 's/\t/_/g' | sort -u | join -t$'\t' - <(cat pcid_entrez_actions.txt | awk 'BEGIN {FS = "\t"} {print $1 "_" $2 "\t" $0}' | sort -k1,1) | cut -f2,3,5 | sort -t$'\t' -k3,3 | join -t$'\t' -1 3 -2 1 - <(cat raw/action_direction.txt | sort -k1,1) | cut -f2-) | sort -u > resolved_mod_to_entrez.txt

cat ../../../import/dgidb.txt | cut -f-2 | sort -u | comm -13 - <(cat resolved_mod_to_entrez.txt | cut -f-2 | sort -u) | sed 's/\t/_/g' | sort -u | join -t$'\t' - <(cat resolved_mod_to_entrez.txt | awk '{print $1 "_" $2 "\t" $0}' | sort -k1,1) | cut -f2- | sort -u > resolved_nodgidb_mod_to_entrez.txt

(cat nogtpdgidb_entrez_pcid_actions.txt | cut -f-3 | sort -t$'\t' -k3,3 | join -t$'\t' -1 3 -2 1 - <(cat raw/type_direction.txt | sort -k1,1) | cut -f2-
cat nogtpdgidb_entrez_pcid_actions.txt | cut -f1,2,4 | sort -t$'\t' -k3,3 | join -t$'\t' -1 3 -2 1 - <(cat raw/type_direction.txt | cut -f1 | sort -u) | cut -f2- | sort -u | comm -13 - <(cat nogtpdgidb_entrez_pcid_actions.txt | cut -f-2 | sort -u) | sed 's/\t/_/g' | sort -u | join -t$'\t' - <(cat nogtpdgidb_entrez_pcid_actions.txt | awk 'BEGIN {FS = "\t"} {print $1 "_" $2 "\t" $0}' | sort -k1,1) | cut -f2,3,5 | sort -t$'\t' -k3,3 | join -t$'\t' -1 3 -2 1 - <(cat raw/action_direction.txt | sort -k1,1) | cut -f2-) | awk 'BEGIN {FS = "\t"} {print $2 "\t" $1 "\t" $3}' | sort -u > nogtpdgidb_mod_to_entrez.txt


cat ../../../import/dgidb.txt | cut -f-2 | sort -u | comm -13 - <(cat nogtpdgidb_mod_to_entrez.txt | cut -f-2 | sort -u) | sed 's/\t/_/g' | sort -u | join -t$'\t' - <(cat nogtpdgidb_mod_to_entrez.txt | awk '{print $1 "_" $2 "\t" $0}' | sort -k1,1) | cut -f2- | sort -u > nodgidb_mod_to_entrez.txt

