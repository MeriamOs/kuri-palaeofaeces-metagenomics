#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=eager_screening
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output eager_screening.%j.out # CHANGE map1 part each run
#SBATCH --error eager_screening.%j.err # CHANGE map1 part each runmodule purge 

########################################
### PREPARATION OF GENOME REFERENCES ###
########################################

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/references

module load datasets/18.23.0

datasets download genome accession GCF_000002285.3 \
    --include genome \
    --filename GCF_000002285.3.zip

unzip GCF_000002285.3.zip

datasets download genome accession GCF_000001405.40 \
    --include genome \
    --filename GCF_000001405.40.zip 

unzip GCF_000001405.40.zip 

datasets download genome accession GCF_014441545.1 --chromosomes Y \
    --include genome \
    --filename GCF_014441545.1.zip

unzip GCF_014441545.1.zip

cp ncbi_dataset/data/GCF_014441545.1/chrY.fna GCF_014441545.1_chrY.fna
cp ncbi_dataset/data/GCF_000002285.3/GCF_000002285.3_CanFam3.1_genomic.fna .
cp ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna .

cat GCF_000002285.3_CanFam3.1_genomic.fna GCF_014441545.1_chrY.fna > GCF_000002285.3_CanFam3.1_chrY_genomic.fna

module load BWA/0.7.18-GCC-12.3.0

bwa index GCF_000002285.3_CanFam3.1_chrY_genomic.fna
bwa index GCF_000001405.40_GRCh38.p14_genomic.fna

#####################
### GET SEQUENCES ###
#####################

/home/mvanos/programs/bs project list

/home/mvanos/programs/bs download project -i <ID_number> --extension=fastq.gz

/home/mvanos/programs/bs download project -i 399153801 --extension=fastq.gz
/home/mvanos/programs/bs download project -i 418160226 --extension=fastq.gz

bash rename_screening_files.sh 

module purge

module load Nextflow/22.04.3

export SINGULARITY_TMPDIR=/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/cache/cache
export SINGULARITY_CACHEDIR=$SINGULARITY_TMPDIR
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_TMPDIR
setfacl -b "$SINGULARITY_TMPDIR"

#nextflow run nf-core/eager -profile test,singularity

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis

###################
### DOG MAPPING ###
###################

nextflow run nf-core/eager \
    -r 2.5.0 \
    -c 02_eager.config \
    -profile singularity \
    --outdir '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/screening/eager_CanFam3.1_screening/' \
    --input '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/03-data/NZ_coprolites/raw/nanoruns/screening_data/*_R{1,2}_001.fastq.gz' \
    --fasta '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/references/GCF_000002285.3_CanFam3.1_chrY_genomic.fna' \
    --bwa_index '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/references/' \
    --single_stranded \
    --bwaalnn 0.03 \
    --dedupper 'dedup' \
    --mergedonly \
    --skip_damage_calculation \
    --run_bam_filtering \
    --bam_mapping_quality_threshold 20

#####################
### HUMAN MAPPING ###
#####################

nextflow run nf-core/eager \
    -r 2.5.0 \
    -c 02_eager.config \
    -profile singularity \
    --outdir '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/screening/eager_GRCh38_screening/' \
    --input '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/screening/eager_CanFam3.1_screening/adapterremoval/output/*_L001_R1_001.fastq_L0.pe.combined.fq.gz' \
    --fasta '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/references/GCF_000001405.40_GRCh38.p14_genomic.fna' \
    --bwa_index '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/references/' \
    --single_end \
    --single_stranded \
    --skip_fastqc --skip_adapterremoval --skip_preseq \
    --bwaalnn 0.03 \
    --dedupper 'dedup' \
    --mergedonly \
    --skip_damage_calculation \
    --run_bam_filtering \
    --bam_mapping_quality_threshold 20 

#clean up intermediate files
#nextflow clean -f -k

