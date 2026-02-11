#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=eager_Meriam
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output eager_mtDNA.%j.out # CHANGE map1 part each run
#SBATCH --error eager_mtDNA.%j.err # CHANGE map1 part each runmodule purge 

module purge

module load Nextflow/22.04.3
module load Singularity/3.11.3

export SINGULARITY_TMPDIR=/nesi/nobackup/uoo02328/meriam/container-cache
export SINGULARITY_CACHEDIR=$SINGULARITY_TMPDIR
export NXF_SINGULARITY_CACHEDIR=$SINGULARITY_TMPDIR
setfacl -b "$SINGULARITY_TMPDIR"

#nextflow run nf-core/eager -r 2.5.0 -profile test,singularity

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis

nextflow run nf-core/eager \
-r 2.5.0 \
-c eager.config \
-profile singularity \
--outdir '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/canis_mtDNA/' \
--input '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/03-data/NZ_coprolites/*_R{1,2}_001.fastq.gz' \
--fasta '/nesi/project/uoo02328/references/domesticates_mtDNA/Canis_lupus_familiaris_mtDNA_NC_002008.fasta' \
--bwaalnn 0.03 --dedupper dedup --mergedonly \
--damage_calculation_tool 'mapdamage' --mapdamage_downsample 100000 \
--run_genotyping --genotyping_tool 'ug' \
--gatk_ploidy 1 --gatk_ug_out_mode 'EMIT_ALL_SITES' --gatk_ug_genotype_model 'SNP' \
--run_vcf2genome --vcf2genome_minc 3 --vcf2genome_minfreq 0.66 \
-resume 

#clean up intermediate files
#nextflow clean -f -k

