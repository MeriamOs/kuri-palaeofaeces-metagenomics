### Script to remove the duplicate sequences of the unmapped reads

module purge
module load fastp/1.1.0-GCC-12.3.0

for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 \
			MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775
	do
    fastp -i ${sample}.unmapped.fastq.gz_lowcomplexityremoved.fq.gz --dedup --dup_calc_accuracy 3 \
        -o fastp/${sample}.unmapped.lowcomplexityremoved.dedup.fq.gz
done

cd fastp

module load FastQC/0.11.9

fastqc *.dedup.fq.gz 

module load MultiQC/1.24.1-foss-2023a-Python-3.11.6

multiqc .
