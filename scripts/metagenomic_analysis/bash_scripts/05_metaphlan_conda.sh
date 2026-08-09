#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J metaphlan
#SBATCH --time 24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output metaphlan_palaeofaeces.%j.out # CHANGE map1 part each run
#SBATCH --error metaphlan_palaeofaeces.%j.err # CHANGE map1 part each runmodule purge

module purge 
module load Miniconda3

source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#conda create -p /nesi/nobackup/uoo02328/meriam/conda_environments/mpa -c conda-forge -c bioconda python=3.7 metaphlan=4.1.2

conda activate /nesi/nobackup/uoo02328/meriam/conda_environments/mpa

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/metaphlan

#metaphlan --install --index mpa_vJun23_CHOCOPhlAnSGB_202403 --bowtie2db metaphlan_databases

DATADIR='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/eager_CanFam3.1/metagenomic_complexity_filter/fastp'
INDEXDIR='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/metaphlan/metaphlan_databases/metaphlan_jun23/'
INDEX='mpa_vJun23_CHOCOPhlAnSGB_202403'

#module load SeqKit/2.4.0

### PALAEOFAECES ###

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
              MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
    do
    metaphlan ${DATADIR}/${sample}.unmapped.lowcomplexityremoved.dedup.fq.gz \
        --input_type fastq \
        --bowtie2out reports_jun23_database/${sample}_jun23_dedup.bt2.out  \
        --nproc 16 \
        --bowtie2db metaphlan_databases \
        --read_min_len 30 \
        -t rel_ab_w_read_stats \
        --bt2_ps very-sensitive \
        --index $INDEX \
        --offline \
        > reports_jun23_database/${sample}_jun23_nreads.metaphlan_profile.txt
done

### BLANKS ###

for sample in KH_blank_1 KH_blank_2 LB_blank_1 LB_blank_2 WH_blank_1 WH_blank_2;
    do
    metaphlan ${DATADIR}/${sample}.unmapped.lowcomplexityremoved.dedup.fq.gz \
        --input_type fastq \
        --bowtie2out reports_jun23_database/${sample}_jun23_dedup.bt2.out  \
        --nproc 16 \
        --bowtie2db metaphlan_databases \
        --read_min_len 30 \
        -t rel_ab_w_read_stats \
        --bt2_ps very-sensitive \
        --index $INDEX \
        --offline \
        > reports_jun23_database/${sample}_jun23_nreads.metaphlan_profile.txt
done

### ENVIRONMENTAL CONTROLS ###

DATADIR2='/nesi/project/uoo02328/meriam/data/decontam_data'
    for sample in MS10790 MS10902 MS10903 MS10904 MS11102 MS11103 MS11107 MS11108;
    do
    metaphlan ${DATADIR2}/${sample}.collapsed.dedup.fq.gz \
        --input_type fastq \
        --bowtie2out reports_jun23_database/${sample}_jun23_dedup.bt2.out  \
        --nproc 4 \
        --bowtie2db metaphlan_databases \
        --read_min_len 30 \
        -t rel_ab_w_read_stats \
        --bt2_ps very-sensitive \
        --index $INDEX \
        --offline \
        > reports_jun23_database/${sample}_jun23_nreads.metaphlan_profile.txt
done
