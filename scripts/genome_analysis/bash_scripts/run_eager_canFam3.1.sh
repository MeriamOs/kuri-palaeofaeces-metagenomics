#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=eager_canFam3.1
#SBATCH --qos=debug
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=ALL
#SBATCH --output eager_canFam3.1.%j.out # CHANGE map1 part each run
#SBATCH --error eager_canFam3.1.%j.err # CHANGE map1 part each runmodule purge 

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
--outdir '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/eager_CanFam3.1/' \
--input '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/canis_mtDNA/adapterremoval/output/MS*.fq.gz' \
--fasta '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/references/Canis_lupus_familiaris.CanFam3.1.dna.toplevel.fa' \
--bwa_index '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/references/' \
--single_end \
--single_stranded \
--skip_fastqc --skip_adapterremoval --skip_preseq \
--run_trim_bam --bamutils_clip_single_stranded_none_udg_left 2 --bamutils_clip_single_stranded_none_udg_right 2 --bamutils_softclip \
--bwaalnn 0.01 \
--dedupper 'dedup' \
--run_mtnucratio \
--mergedonly \
--damage_calculation_tool 'mapdamage' --mapdamage_downsample 100000 \
--run_bam_filtering --bam_unmapped_type 'fastq' --metagenomic_complexity_filter \
--run_genotyping --genotyping_tool 'pileupcaller' --genotyping_source 'trimmed' \
--pileupcaller_bedfile '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/genotyping/souilmi_files/souilmi_positions.txt' \
--pileupcaller_snpfile '/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/genotyping/souilmi_files/souilmi2024.snp' \
--pileupcaller_method 'randomHaploid' \
--run_bcftools_stats \
-resume 

#clean up intermediate files
#nextflow clean -f -k

