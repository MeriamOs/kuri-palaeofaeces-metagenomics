module load SAMtools

for ref in $(grep '>' ../references/bacteria_all.fasta | awk '{print substr($1,2)}'); 
    do
        echo "$ref"
        cat header.txt | sed "s/#rname/$ref/" > ${ref}_stats.txt
        for sample in MS11669 MS11670 MS11673 MS11674 MS11675 MS11676 MS11677 MS11678 MS11679 MS11683 MS11684 MS11686 MS11770 MS11771 MS11774 MS11775;
            do 
            samtools coverage ${sample}_maponly.bam | grep $ref | sed "s/${ref}/$sample/" >> ${ref}_stats.txt
        done
done
