#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=bamplotter
#SBATCH --time=00:10:00
#SBATCH --mem=50G
#SBATCH --mail-type=ALL
#SBATCH --output bamplotter.%j.out # CHANGE map1 part each run
#SBATCH --error bamplotter.%j.err # CHANGE map1 part each runmodule purge 

#module purge
module load Miniconda3

source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#conda create -p /nesi/nobackup/uoo02328/meriam/conda_environments/pydamage -c bioconda -c conda-forge pydamage 

conda activate /nesi/nobackup/uoo02328/meriam/conda_environments/pydamage

MARKER='mitochondrion'
    
for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775
	do
	echo ${sample}
	python BAMPlotter_pydamage_10readsonly.py \
	-b ${MARKER}/${sample}_maponly.bam \
    -d ${MARKER}/${sample}_${MARKER}_pydamage_results.csv \
	-o ${MARKER}/BAMPlotter/${sample}_BAMPlotter_${MARKER}_pydamage_10+reads.pdf
	done 

MARKER='COI_region'

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775
	do
	echo ${sample}
	python BAMPlotter_pydamage_COI.py \
	-b ${MARKER}/${sample}_maponly.bam \
    -d ${MARKER}/${sample}_${MARKER}_pydamage_results.csv \
	-o ${MARKER}/BAMPlotter/${sample}_BAMPlotter_${MARKER}_pydamage_10+reads.pdf
	done 

conda deactivate	
