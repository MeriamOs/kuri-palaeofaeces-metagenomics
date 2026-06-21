#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J krk2
#SBATCH --time 4:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --mem=90G
#SBATCH --mail-type=ALL
#SBATCH --output krk2_palaeofaeces.%j.out # CHANGE map1 part each run
#SBATCH --error krk2_palaeofaeces.%j.err # CHANGE map1 part each run

module purge

module load Kraken2
module load BLAST/2.13.0-GCC-11.3.0

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis

DATADIR='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/eager_CanFam3.1/metagenomic_complexity_filter/fastp'
DATADIR2='/nesi/project/uoo02328/meriam/data/decontam_data'

#########################
## Reference Database ###
#########################

## standard database obtained by:
# wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz
# tar -xvzf k2_standard_20240904.tar.gz 

#################
## Run kraken2 ##
#################

ref='/nesi/nobackup/uoo02328/meriam/kraken/databases/2024-standard'

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 
			MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775 \
			WH_blank_1 WH_blank_2 KH_blank_1 KH_blank_2 LB_blank_1 LB_blank_2;
	do
	echo ${sample}
	kraken2 --db $ref \
		--threads 32 \
		--classified-out kraken/kraken-2024-standard/${sample}_dedup_kraken_2024-standard_conf0.50.fasta \
		--report kraken/kraken-2024-standard/${sample}_dedup_kraken_2024-standard_conf0.50.txt \
		--output kraken/kraken-2024-standard/${sample}_dedup_kraken_2024-standard_conf0.50 \
		--use-names \
		--report-minimizer-data \
		--gzip-compressed \
		--minimum-base-quality 30 \
		--confidence 0.50 \
		${DATADIR}/${sample}.unmapped.fastq.gz_lowcomplexityremoved.dedup.fq.gz
done

for sample in MS10790 MS10902 MS10903 MS10904 MS11102 MS11103 MS11107 MS11108;
	do
	echo ${sample}
	kraken2 --db $ref \
		--threads 32 \
		--classified-out kraken/kraken-2024-standard/${sample}_dedup_kraken_2024-standard_conf0.50.fasta \
		--report kraken/kraken-2024-standard/${sample}_dedup_kraken_2024-standard_conf0.50.txt \
		--output kraken/kraken-2024-standard/${sample}_dedup_kraken_2024-standard_conf0.50 \
		--use-names \
		--report-minimizer-data \
		--gzip-compressed \
		--minimum-base-quality 30 \
		--confidence 0.50 \
		${DATADIR2}/${sample}.collapsed.dedup.fq.gz
done
