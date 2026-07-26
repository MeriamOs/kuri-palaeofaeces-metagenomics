#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J coproid
#SBATCH --time 0:10:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3G
#SBATCH --mail-type=ALL
#SBATCH --output sourcepredict.%j.out # CHANGE map1 part each run
#SBATCH --error sourcepredict.%j.err # CHANGE map1 part each runmodule purge

module purge 
module load Miniconda3

source $(conda info --base)/etc/profile.d/conda.sh
export PYTHONNOUSERSITE=1

#conda create -p /nesi/nobackup/uoo02328/meriam/conda_environments/sourcepredict python=3.10 -c conda-forge -c bioconda sourcepredict

conda activate /nesi/nobackup/uoo02328/meriam/conda_environments/sourcepredict

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/sourcepredict

sourcepredict -s sp_sources_september2025.csv -l sp_labels_september2025.csv \
    -o coproid.kuri.standard2024.june2026.1000cutoff.sourcepredict.csv \
    -e coproid.kuri.standard2024.june2026.1000cutoff.embedding.csv -t 8 \
    kraken_kuri_june2026_1000count_minimizer.csv 

conda deactivate
