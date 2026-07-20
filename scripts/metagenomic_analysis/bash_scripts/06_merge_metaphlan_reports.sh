module load MetaPhlAn/4.1.0-gimkl-2022a-Python-3.10.5

### MERGE SAMPLE TABLE BASED ON SIZE, AND OTHER DATA CATEGORIES ###

merge_metaphlan_tables.py reports_negatives/*.metaphlan_profile.txt > merged_negatives_table.txt
merge_metaphlan_tables.py reports_jun23_database/MS*.metaphlan_profile.txt > merged_palaeofaeces_table.txt
merge_metaphlan_tables.py reports_jun23_database/SR*.metaphlan_profile.txt > merged_modern_table.txt

#create species only table
grep -E "s__|MS" merged_palaeofaeces_table.txt \
| grep -v "t__" \
| sed "s/^.*|//g" \
| sed "s/.metaphlan//g" \
> merged_palaeofaeces_metaphlan_table_species.txt

grep -E "s__|MS" merged_negatives_table.txt \
| grep -v "t__" \
| sed "s/^.*|//g" \
| sed "s/.metaphlan//g" \
> merged_negatives_metaphlan_table_species.txt

grep -E "s__|SR" merged_modern_table.txt \
| grep -v "t__" \
| sed "s/^.*|//g" \
| sed "s/_4.1.0.metaphlan//g" \
> merged_modern_table_species.txt
