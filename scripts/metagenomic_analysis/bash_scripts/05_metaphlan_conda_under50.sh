#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J metaphlan
#SBATCH --time 6:00:00
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

#MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775 ERR5863536 SRR12455959 ERR3761407 ERR3761412

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/metaphlan

INDEXDIR='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/metaphlan/metaphlan_databases/metaphlan_jun23/'
INDEX='mpa_vJun23_CHOCOPhlAnSGB_202403'

#module load SeqKit/2.4.0

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
              MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
do
seqkit seq ${DATADIR}/${sample}.unmapped.fastq.gz_lowcomplexityremoved.fq.gz \
    --max-len 50 > bowtie2/${sample}_under50bp_filtered.fq
bowtie2 -x ${INDEXDIR}/${INDEX} \
    -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 \
    -p 16 \
    -U bowtie2/${sample}_under50bp_filtered.fq \
    -S bowtie2/${sample}_under50bp.sam
done


MS11669_nread=46541918
MS11670_nread=21260786
MS11673_nread=33907031
MS11674_nread=29401446
MS11675_nread=48303933
MS11676_nread=44088490
MS11677_nread=40862354
MS11678_nread=39290741
MS11679_nread=34964342
MS11683_nread=46908140
MS11684_nread=23459306
MS11686_nread=25848236
MS11770_nread=27043303
MS11771_nread=33468829
MS11774_nread=45864468
MS11775_nread=40356435


#MS11669 MS11670 
for sample in MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
    MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
do
value=$(eval echo \${${sample}_nread})
metaphlan bowtie2/${sample}_under50bp.sam \
    --input_type sam \
    --nproc 4 \
    -t rel_ab_w_read_stats \
    --bowtie2db metaphlan_databases/metaphlan_jun23 \
    --offline \
    --nreads $value \
    -o reports_under50bp/${sample}_jun23_under50bp.metaphlan_profile.txt
done
