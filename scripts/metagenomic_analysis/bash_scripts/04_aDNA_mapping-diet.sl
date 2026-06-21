#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=dogs_map
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-user=meriamvanos96@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output plants.%j.out # CHANGE map1 part each run
#SBATCH --error plants.%j.err # CHANGE map1 part each run


##slurm job to run aDNA_trimmer.sh and aDNA_mapper.sh to map to reference genome, with same parameters as ancientDNA_trimQC.sh

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis

bash aDNA_mapper.sh 04A_mapping_parameters-diet_mtDNA.txtbash 
aDNA_mapper.sh 04B_mapping_parameters-diet_COI.txt
bash aDNA_mapper.sh 04C_mapping_parameters-diet_plants.txt
