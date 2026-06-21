#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=bamplotter
#SBATCH --time=2:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-type=ALL
#SBATCH --output bamplotter.%j.out # CHANGE map1 part each run
#SBATCH --error bamplotter.%j.err # CHANGE map1 part each runmodule purge 

###################
## MAPPING STATS ##
###################

module load SAMtools/1.12-GCC-9.2.0

for ref in $(grep '>' ../../references/bacteria_all.fasta | awk '{print substr($1,2)}'); 
    do
        echo "$ref"
        cat header.txt | sed "s/#rname/$ref/" > ${ref}_stats.txt
        for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
                        MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
            do 
            samtools coverage ${sample}_maponly.bam | grep $ref | sed "s/${ref}/$sample/" >> ${ref}_stats.txt
        done
done

####################
## ACTIVATE CONDA ##
####################

module purge
module load Miniconda3

source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

##################
## Run pydamage ##
##################

#conda create -p /nesi/nobackup/uoo02328/meriam/conda_environments/pydamage -c bioconda -c conda-forge pydamage 

conda activate /nesi/nobackup/uoo02328/meriam/conda_environments/pydamage
MARKER='all_bacteria'

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
            MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775
    do 
    pydamage analyze -f --no_ga bam_files/${sample}_maponly.bam
    mv pydamage_results/pydamage_results.csv ${sample}_${MARKER}_pydamage_results.csv
done

conda deactivate

####################
## Run BAMPlotter ##
####################

#conda create -p /nesi/nobackup/uoo02328/meriam/conda_environments/pysamstats_env \
#  -c bioconda -c conda-forge \
#  python=3.7 pysam=0.15.4 pip setuptools wheel

conda activate /nesi/nobackup/uoo02328/meriam/conda_environments/pysamstats_env
#pip install pysamstats --no-build-isolation
#conda install -c conda-forge pandas matplotlib tqdm

MARKER='all_bacteria'

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
			MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775
	do
	echo ${sample}
	/nesi/nobackup/uoo02328/meriam/conda_environments/pysamstats_env/bin/python BAMPlotter_pydamage_1000readsonly.py \
		-b bam_files/${sample}_maponly.bam \
    	            -d pydamage_results/${sample}_${MARKER}_pydamage_results.csv \
		-o BAMPlotter/${sample}_${MARKER}_BAMPlotter_pydamage_1000+reads.pdf
done 

conda deactivate	
