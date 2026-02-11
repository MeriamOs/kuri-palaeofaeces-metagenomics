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

conda activate sourcepredict

sourcepredict -s sp_sources_september2025.csv -l sp_labels_september2025.csv \
    -o coproid.kuri.standard2024.september2025.50cutoff.sourcepredict.csv \
    -e coproid.kuri.standard2024.sepetmber2025.50cutoff.embedding.csv -t 8 \
    kraken_reports/kraken_kuri_50cutoff.csv 