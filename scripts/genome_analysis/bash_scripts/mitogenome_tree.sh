# Script to download coyote and dog references and create a ML tree
# Created by Meriam van Os

module load entrez-direct/13.3

for acc in KM061574	KM061588 EU789740 DQ480496 EU789756 KM061501 EU789721 EU408280 \
            EU789759 KM061506 KM061529 DQ480497 KU290442 KM061567 EU789686 EU789717 \
            KU290608 KM061542 EU789677 KU290617 EU789673 KT168369 EU789776 EU789681 \
            KM113774 EU789692 EU789669 EU789666 EU408300 EU789651 EU789757 EU408268 \
            EU789646 KU290510 KJ637137 KU290450 KX379529 EU789655 EU789773 EU789662 \
            AB499817 NC_008093.1; #42
        do
        esearch -db nucleotide -query "$acc" | \
        efetch -format fasta > "${acc}.fasta"
done

for acc in KT168372 KT168373 KT168374 KU215693 KU215698 KU215702 KY798508 KY798509 KY798510 KY798516; #10
        do
        esearch -db nucleotide -query "$acc" | \
        efetch -format fasta > "${acc}.fasta"
done

cat ind_fastas/*.fasta >> coyote_dogs_palaeofaeces.fasta
cat kuri_all_NZ_June_2024_aligned.fasta >> coyote_dogs_palaeofaeces.fasta

module purge
module load MAFFT/7.505-gimkl-2022a-with-extensions

mafft coyote_dogs_palaeofaeces.fasta > coyote_dogs_palaeofaeces_aligned.fasta

module purge
module load seqmagick/0.8.4-gimkl-2020a-Python-3.8.2

seqmagick convert coyote_dogs_palaeofaeces_aligned.fasta coyote_dogs_palaeofaeces_aligned.phy

module purge 
module load PhyML/3.3.20211231-gimpi-2022a

phyml -i coyote_dogs_palaeofaeces_aligned.phy
