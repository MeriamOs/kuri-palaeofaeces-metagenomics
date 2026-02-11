module purge
module load Bracken/2.7-GCC-11.3.0
module load Python/3.11.3-gimkl-2022a

DB='/nesi/nobackup/uoo02328/meriam/kraken/databases/2024-standard'
prefix=2024-standard

for sample in MS10790 MS10902 MS10903 MS10904 MS11102 MS11103 MS11107 MS11108 \
              Blank1_WH Blank2_WH KH_blank_1 KH_blank_2 LB_blank_1 LB_blank_2 \
              MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
              MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
do 
python filter_kraken_minimizers.py reports/${sample}_kraken_${prefix}_conf0.50.txt \
    reports/${sample}_kraken_${prefix}_conf0.50_minimizer.txt
done

for sample in MS10790 MS10902 MS10903 MS10904 MS11102 MS11103 MS11107 MS11108 \
              Blank1_WH Blank2_WH KH_blank_1 KH_blank_2 LB_blank_1 LB_blank_2;
do
bracken -d $DB -i reports/${sample}_kraken_${prefix}_conf0.50_minimizer.txt \
    -o bracken_decontam/${sample}_minimizer.bracken -l S
python ../KrakenTools/filter_bracken.out.py -i bracken_decontam/${sample}_minimizer.bracken \
    --exclude 9606 -o bracken_decontam/${sample}_minimizer_filtered.bracken
done

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
              MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
do
bracken -d $DB -i reports/${sample}_kraken_${prefix}_conf0.50_filtered.txt \
        -o bracken_species/${sample}_minimizer.bracken -l S
python ../KrakenTools/filter_bracken.out.py -i bracken_species/${sample}_minimizer.bracken \
        --exclude 9606 -o bracken_species/${sample}_minimizer_filtered.bracken
done

python combine_bracken_outputs.py \
--files bracken_decontam/*minimizer_filtered.bracken \
-o negatives_combined_${prefix}_minimizer_bracken_species.bracken

python combine_bracken_outputs.py \
--files bracken_species/*minimizer_filtered.bracken \
-o palaeofaeces_combined_${prefix}_minimizer_bracken_species.bracken

#python tsv_to_csv.py < input_file.tsv > output_file.csv
#python tsv_to_csv.py < palaeofaeces_modern_combined_bracken_phylum.bracken > palaeofaeces_modern_combined_bracken_phylum.csv
#python tsv_to_csv.py < palaeofaeces_modern_combined_bracken_genus.bracken > palaeofaeces_modern_combined_bracken_genus.csv

python tsv_to_csv.py < palaeofaeces_combined_${prefix}_minimizer_bracken_species.bracken \
        > palaeofaeces_combined_${prefix}_minimizer_bracken_species.csv
python tsv_to_csv.py < negatives_combined_${prefix}_minimizer_bracken_species.bracken \
        > negatives_combined_${prefix}_minimizer_bracken_species.csv

## csv files then get further processed in RStudio
