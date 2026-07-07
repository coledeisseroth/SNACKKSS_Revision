(
echo -n $'Lines that were unchanged\t'
comm -12 <(cat ../import/corrected_curated_dataset.txt | awk 'NR > 1' | sort -u) <(cat ../import/curated_dataset.txt | awk 'NR > 1' | sort -u) | wc -l

echo -n $'Lines that were changed\t'
comm -13 <(cat ../import/corrected_curated_dataset.txt | awk 'NR > 1' | sort -u) <(cat ../import/curated_dataset.txt | awk 'NR > 1' | sort -u) > remainder.txt
cat remainder.txt | wc -l

echo -n $'Lines with a stray newline character\t'
cat remainder.txt | grep -v GSE | wc -l
cat remainder.txt | grep GSE > tmp; mv tmp remainder.txt

echo -n $'Lines with a stray TAB character\t'
cat remainder.txt | grep $'\t"\t' | wc -l
cat remainder.txt | grep -v $'\t"\t' | awk 'BEGIN {FS = "\t"} NF == 12' > tmp; mv tmp remainder.txt

echo -n $'Lines where the experiment was not using RNA-Seq\t'
cat ../import/corrected_curated_dataset.txt | grep NOT_RNA | awk 'BEGIN {FS = "\t"} {print $2 "_" $4 "_" $6 "_" $7 "_" $11}' | sort -t$'\t' -k1,1 | join -t$'\t' - <(cat remainder.txt | awk 'BEGIN {FS = "\t"} {print $2 "_" $4 "_" $6 "_" $7 "_" $11 "\t" $0}' | sort -t$'\t' -k1,1) | cut -f2- | sort -u | wc -l
cat ../import/corrected_curated_dataset.txt | grep NOT_RNA | awk 'BEGIN {FS = "\t"} {print $2 "_" $4 "_" $6 "_" $7 "_" $11}' | sort | comm -13 - <(cat remainder.txt | awk 'BEGIN {FS = "\t"} {print $2 "_" $4 "_" $6 "_" $7 "_" $11}' | sort) | sort | join -t$'\t' - <(cat remainder.txt | awk 'BEGIN {FS = "\t"} {print $2 "_" $4 "_" $6 "_" $7 "_" $11 "\t" $0}' | sort -t$'\t' -k1,1) | cut -f2- | sort -u > tmp; mv tmp remainder.txt

echo -n $'Lines with the outdated GD experiment label\t'
cat remainder.txt | awk 'BEGIN {FS = "\t"} $4 == "GD"' | wc -l
cat remainder.txt | awk 'BEGIN {FS = "\t"} $4 != "GD"' > tmp; mv tmp remainder.txt

echo -n $'Lines where the comments were changed\t'
cat remainder.txt | cut -f-2,4-11 | sort | comm -12 - <(cat ../import/corrected_curated_dataset.txt | cut -f-2,4-11 | sort) | wc -l
cat remainder.txt | cut -f-2,4-11 | sort | comm -23 - <(cat ../import/corrected_curated_dataset.txt | cut -f-2,4-11 | sort) > tmp; mv tmp remainder.txt

echo -n $'Lines where the perturbation type was changed\t'
cat remainder.txt | awk 'BEGIN {FS = "\t"} {print $2 "_" $5 "_" $6 "\t" $3}' | sort -t$'\t' -k1,1 | join -t$'\t' - <(cat ../import/corrected_curated_dataset.txt | awk 'BEGIN {FS = "\t"} {print $2 "_" $6 "_" $7 "\t" $4}' | sort -t$'\t' -k1,1) | awk '$2 != $3' | wc -l
cat remainder.txt | awk 'BEGIN {FS = "\t"} {print $2 "_" $5 "_" $6 "\t" $3}' | sort -t$'\t' -k1,1 | join -t$'\t' - <(cat ../import/corrected_curated_dataset.txt | awk 'BEGIN {FS = "\t"} {print $2 "_" $6 "_" $7 "\t" $4}' | sort -t$'\t' -k1,1) | awk '$2 != $3' | cut -f1 | sort -u | comm -13 - <(cat remainder.txt | awk 'BEGIN {FS = "\t"} {print $2 "_" $5 "_" $6}' | sort -u) | join -t$'\t' - <(cat remainder.txt | awk 'BEGIN {FS = "\t"} {print $2 "_" $5 "_" $6 "\t" $0}' | sort -t$'\t' -k1,1) | cut -f2- > tmp; mv tmp remainder.txt

echo -n $'Lines where the labels were changed\t'
cat remainder.txt | wc -l
) > correction_statistics.txt

#Gather the remaining differing lines
#diff <(cat remainder.txt | cut -f1 | sort -u | join -t$'\t' - <(cat curated_dataset.txt | sort -k1,1) | sort) <(cat remainder.txt | cut -f1 | sort -u | join -t$'\t' - <(cat corrected_curated_dataset.txt | sort -k1,1) | sort)


