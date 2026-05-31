cd chembl
bash BUILD.sh
cd ../drugrepurposinghub
bash BUILD.sh
cd ../guidetopharmacology
bash BUILD.sh
cd ../pharmgkb
bash BUILD.sh
cd ../trrust
bash BUILD.sh
cd ..

cat chembl/nodgidb_mod_to_entrez.txt drugrepurposinghub/nodgidb_mod_to_entrez.txt guidetopharmacology/nodgidb_mod_to_entrez.txt pharmgkb/nodgidb_mod_to_entrez.txt | sort -u > nodgidb_aggregated_mod_to_entrez.txt

