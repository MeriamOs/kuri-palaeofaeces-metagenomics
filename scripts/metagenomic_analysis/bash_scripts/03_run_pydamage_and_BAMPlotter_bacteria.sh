#!/bin/bash -e
#SBATCH --account=uoo02328
#SBATCH --job-name=bamplotter
#SBATCH --time=2:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64G
#SBATCH --mail-type=ALL
#SBATCH --output bamplotter.%j.out # CHANGE map1 part each run
#SBATCH --error bamplotter.%j.err # CHANGE map1 part each runmodule purge 

##############
## BACTERIA ##
##############

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/bacterial_mapping/all_bacteria
mkdir mapping_stats
mkdir BAMPlotter
mkdir pydamage_stats

#####################
## EXTRA FILTERING ##
#####################

module load SAMtools/1.23.1-GCC-12.3.0

mkdir bam_files
cd bam_files

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775
	do
    samtools view -b -S -q 30 ../${sample}_maponly.bam > ${sample}_maponly_q30.bam
done

for bam in *_q30.bam
    do
    samtools index ${bam}
done

cd ..

###################
## MAPPING STATS ##
###################

module load SAMtools/1.12-GCC-9.2.0

samtools coverage MS11669_maponly.bam | head -n 1 > header.txt

for ref in $(grep '>' ../../references/bacteria_all.fasta | awk '{print substr($1,2)}'); 
    do
        echo "$ref"
        cat header.txt | sed "s/#rname/$ref/" > ${ref}_stats.txt
        for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
                        MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
            do 
            samtools coverage ${sample}_maponly.bam | grep $ref | sed "s/${ref}/$sample/" >> mapping_stats/${ref}_stats.txt
        done
done

for ref in $(grep '>' ../../references/bacteria_all.fasta | awk '{print substr($1,2)}'); 
    do
        echo "$ref"
        cat header.txt | sed "s/#rname/$ref/" > ${ref}_stats.txt
        for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
                        MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
            do 
            samtools coverage bam_files/${sample}_maponly_q30.bam | grep $ref | sed "s/${ref}/$sample/" >> mapping_stats/${ref}_q30_stats.txt
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
    pydamage analyze -f --no_ga ${sample}_maponly.bam
    mv pydamage_results/pydamage_results.csv pydamage_stats/${sample}_${MARKER}_pydamage_results.csv
done

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
            MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775
    do 
    pydamage analyze -f --no_ga bam_files/${sample}_maponly_q30.bam
    mv pydamage_results/pydamage_results.csv pydamage_stats/${sample}_q30_${MARKER}_pydamage_results.csv
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
			MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
	do
	echo ${sample}
	/nesi/nobackup/uoo02328/meriam/conda_environments/pysamstats_env/bin/python BAMPlotter_pydamage_1000readsonly.py \
		-b ${sample}_maponly.bam \
    	-d pydamage_stats/${sample}_${MARKER}_pydamage_results.csv \
		-o BAMPlotter/${sample}_${MARKER}_BAMPlotter_pydamage_1000+reads.pdf
done 

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
			MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
	do
	echo ${sample}
	/nesi/nobackup/uoo02328/meriam/conda_environments/pysamstats_env/bin/python BAMPlotter_pydamage_1000readsonly.py \
		-b bam_files/${sample}_maponly_q30.bam \
    	-d pydamage_stats/${sample}_q30_${MARKER}_pydamage_results.csv \
		-o BAMPlotter/${sample}_${MARKER}_q30_BAMPlotter_pydamage_1000+reads.pdf
done 

conda deactivate	
