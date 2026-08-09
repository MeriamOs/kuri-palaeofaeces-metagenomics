#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=eager_canFam3.1
#SBATCH --time=128:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output eager_canFam3.1.%j.out # CHANGE map1 part each run
#SBATCH --error eager_canFam3.1.%j.err # CHANGE map1 part each runmodule purge 

module purge

module load Nextflow/22.04.3

export SINGULARITY_TMPDIR=/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/cache/cache
export SINGULARITY_CACHEDIR=$SINGULARITY_TMPDIR
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_TMPDIR
setfacl -b "$SINGULARITY_TMPDIR"

#nextflow run nf-core/eager -profile test,singularity

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis

nextflow run nf-core/eager \
    -r 2.5.0 \
    -c 02_eager.config \
    -profile singularity \
    --outdir '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/eager_CanFam3.1/' \
    --input '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/03-data/NZ_coprolites/*_R{1,2}_001.fastq.gz' \
    --fasta '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/references/GCF_000002285.3_CanFam3.1_chrY_genomic.fna' \
    --bwa_index '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/references/' \
    --single_stranded \
    --bwaalnn 0.03 \
    --dedupper 'dedup' \
    --mergedonly \
    --damage_calculation_tool 'mapdamage' --mapdamage_downsample 100000 \
    --run_bam_filtering --bam_mapping_quality_threshold 20 --bam_unmapped_type 'fastq' --metagenomic_complexity_filter \
    --run_bcftools_stats \
    -resume 

#clean up intermediate files
#nextflow clean -f -k
