#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=bacteria_map
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --output bacterial_mapping.%j.out # CHANGE map1 part each run
#SBATCH --error bacterial_mapping.%j.err # CHANGE map1 part each run

##slurm job to run aDNA_trimmer.sh and aDNA_mapper.sh to map to reference genome, with same parameters as ancientDNA_trimQC.sh

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis

#bash /nesi/project/uoo02328/programs/microbiome_capture/scripts/aDNA_trimmer.sh trim_QC_parameters.txt

bash aDNA_mapper.sh mapping_parameters-bacterial.txt




