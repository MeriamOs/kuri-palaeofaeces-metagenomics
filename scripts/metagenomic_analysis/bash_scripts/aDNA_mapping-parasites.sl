#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=dogs_map
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-user=meriamvanos96@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output parasites.%j.out # CHANGE map1 part each run
#SBATCH --error parasites.%j.err # CHANGE map1 part each run


##slurm job to run aDNA_trimmer.sh and aDNA_mapper.sh to map to reference genome, with same parameters as ancientDNA_trimQC.sh

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis

#bash /nesi/project/uoo02328/programs/microbiome_capture/scripts/aDNA_trimmer.sh trim_QC_parameters.txt

bash aDNA_mapper.sh mapping_parameters-parasites.txt




