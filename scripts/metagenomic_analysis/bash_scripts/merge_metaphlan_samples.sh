module load MetaPhlAn/4.1.0-gimkl-2022a-Python-3.10.5

merge_metaphlan_tables.py *4.1.0.metaphlan_profile.txt > merged_abundance_table.txt

#create species only table
grep -E "s__|MS" merged_abundance_table.txt \
| grep -v "t__" \
| sed "s/^.*|//g" \
| sed "s/_4.1.0.metaphlan//g" \
> merged_abundance_table_species.txt

#create genus only table
grep -E "g__|MS" merged_abundance_table.txt \
| grep -v "t__" \
| sed "s/_4.1.0.metaphlan//g" \
| sed "s/^.*|//g" \
> merged_abundance_table_genus.txt

#create phylum only table
grep -E "p__|MS" merged_abundance_table.txt \
| grep -v "t__" \
| sed "s/^.*|//g" \
| sed "s/_4.1.0.metaphlan//g" \
> merged_abundance_table_phylum.txt
