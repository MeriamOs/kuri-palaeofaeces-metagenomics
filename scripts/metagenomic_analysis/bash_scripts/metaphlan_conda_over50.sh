#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J metaphlan
#SBATCH --time 48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output metaphlan_palaeofaeces.%j.out # CHANGE map1 part each run
#SBATCH --error metaphlan_palaeofaeces.%j.err # CHANGE map1 part each runmodule purge

module purge 
module load Miniconda3

source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#conda create --name mpa -c conda-forge -c bioconda python=3.7 metaphlan

conda activate mpa

DATADIR='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/eager/metagenomic_complexity_filter'

#metaphlan --install --bowtie2db metaphlan/metaphlan_databases

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/metaphlan

INDEXDIR='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/metaphlan/metaphlan_databases/metaphlan_jun23/'
INDEX='mpa_vJun23_CHOCOPhlAnSGB_202403'

module load SeqKit/2.4.0

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
              MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
do
seqkit seq ${DATADIR}/${sample}.unmapped.fastq.gz_lowcomplexityremoved.fq.gz \
    --min-len 50 > bowtie2/${sample}_50bp_filtered.fq
bowtie2 -x ${INDEXDIR}/${INDEX} \
    -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
    -p 8 \
    -U bowtie2/${sample}_50bp_filtered.fq \
    -S bowtie2/${sample}_over50bp.sam
done

MS11669_nread=23143727
MS11670_nread=12166638
MS11673_nread=20452543
MS11674_nread=7041742
MS11675_nread=26515180
MS11676_nread=16224865
MS11677_nread=15346193
MS11678_nread=19996252
MS11679_nread=25568474
MS11683_nread=28024789
MS11684_nread=12793039
MS11686_nread=2778827
MS11770_nread=15437518
MS11771_nread=19579440
MS11774_nread=22172485
MS11775_nread=21080790

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
    MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
do
value=$(eval echo \${${sample}_nread})
metaphlan bowtie2/${sample}_over50bp.sam \
    --input_type sam \
    --nproc 4 \
    -t rel_ab_w_read_stats \
    --bowtie2db metaphlan_databases/metaphlan_jun23 \
    --offline \
    --nreads $value \
    -o reports_over50bp/${sample}_jun23_over50bp.metaphlan_profile.txt
done

for sample in Blank1_WH Blank2_WH KH_blank_1 KH_blank_2 LB_blank_1 LB_blank_2 ;
do
metaphlan ${DATADIR}/${sample}.unmapped.fastq.gz_lowcomplexityremoved.fq.gz \
    --input_type fastq \
    --bowtie2out reports_jun23_database/${sample}_jun23_nreads.bt2.out  \
    --nproc 4 \
    --bowtie2db metaphlan_databases/metaphlan_jun23 \
    --read_min_len 30 \
    -t rel_ab_w_read_stats \
    --index $INDEX \
    --offline \
    -t rel_ab_w_read_stats \
    > reports_jun23_database/${sample}_jun23_nreads.metaphlan_profile.txt
done

DATADIR2='/nesi/project/uoo02328/meriam/data/decontam_data'
for sample in MS10790 MS10902 MS10903 MS10904 MS11102 MS11103 MS11107 MS11108;
do
metaphlan ${DATADIR2}/${sample}.collapsed.gz \
    --input_type fastq \
    --bowtie2out reports_jun23_database/${sample}_jun23_nreads.bt2.out  \
    --nproc 4 \
    --bowtie2db metaphlan_databases/metaphlan_jun23 \
    --read_min_len 30 \
    -t rel_ab_w_read_stats \
    --bt2_ps very-sensitive \
    --index $INDEX \
    --offline \
    -t rel_ab_w_read_stats \
    > reports_jun23_database/${sample}_jun23_nreads.metaphlan_profile.txt
done

#Laos
DATADIR3='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/laos_modern/trimmed_data'
for sample in SRR14842426 SRR14842338 SRR14842349 SRR14842385 SRR14842352 \
        SRR14842396 SRR14842341 SRR14842328 SRR14842325 SRR14842339;
do
metaphlan ${DATADIR3}/${sample}_trimmed/${sample}.collapsed.gz \
    --input_type fastq \
    --bowtie2out reports_jun23_database/${sample}_jun23_nreads.bt2.out  \
    --nproc 4 \
    --bowtie2db metaphlan_databases/metaphlan_jun23 \
    --read_min_len 30 \
    -t rel_ab_w_read_stats \
    --index $INDEX \
    --offline \
    -t rel_ab_w_read_stats \
    > reports_jun23_database/${sample}_jun23_nreads.metaphlan_profile.txt
done

#Inda
DATADIR4='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/yarlagadda_modern'
for sample in SRR14842364 SRR14842376 SRR14842401 SRR14842403 SRR14842389 \
        SRR14842361 SRR14842368 SRR14842383 SRR14842399 SRR14842362;
do
metaphlan ${DATADIR4}/${sample}_trimmed/${sample}.collapsed.gz \
    --input_type fastq \
    --bowtie2out reports_jun23_database/${sample}_jun23_nreads.bt2.out  \
    --nproc 4 \
    --bowtie2db metaphlan_databases/metaphlan_jun23 \
    --read_min_len 30 \
    -t rel_ab_w_read_stats \
    --index $INDEX \
    --offline \
    -t rel_ab_w_read_stats \
    > reports_jun23_database/${sample}_jun23_nreads.metaphlan_profile.txt
done
