#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=dogs_map
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-user=meriamvanos96@hotmail.com
#SBATCH --mail-type=ALL
#SBATCH --output laos.%j.out # CHANGE map1 part each run
#SBATCH --error laos.%j.err # CHANGE map1 part each run


#cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/laos_modern

bash /nesi/project/uoo02328/programs/microbiome_capture/scripts/aDNA_trimmer.sh mapping_parameters.txt

bash modern_mapper.sh mapping_parameters_modern.txt


