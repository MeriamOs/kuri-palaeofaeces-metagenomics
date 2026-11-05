module load SAMtools/1.12-GCC-9.2.0

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/bacterial_mapping/all_bacteria

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
    MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
    do 
    samtools coverage ${sample}_maponly.bam >> map_stats_per_sample/${sample}_mapping_stats.txt
done