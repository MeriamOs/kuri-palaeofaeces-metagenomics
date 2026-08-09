# Script to download bacterial references found with kraken2 and prepare for mapping

module load datasets/18.23.0

for genome in GCF_000005845.2 GCF_000006765.1 GCF_000012825.1 GCF_000013425.1 \
        GCF_000026745.1 GCF_000253315.1 GCF_000317975.2 GCF_000412675.1 GCF_000468955.1 \
        GCF_000756775.1 GCF_001278715.1 GCF_002024345.1 GCF_002393505.1 GCF_002838165.1 \
        GCF_003544815.1 GCF_004328515.1 GCF_008151785.1 GCF_009911755.1 GCF_010537335.1 \
        GCF_011040435.1 GCF_015679285.1 GCF_016406325.1 GCF_017161405.1 GCF_017309675.2 \
        GCF_017915135.1 GCF_018437745.1 GCF_019355055.1 GCF_019880365.1 GCF_023499275.1 \
        GCF_030848765.1 GCF_036884295.1 GCF_049959445.1 GCF_055413515.1 GCF_900002575.1 \
        GCF_900637935.1 GCF_011064845.1 GCF_025147325.1 GCF_029024925.1 GCF_002443155.1 \
        GCF_002949495.1 GCF_048319845.1 GCF_000013285.1 GCF_025723065.1 #43
do
    datasets download genome accession "$genome" \
        --include genome \
        --filename "${genome}.zip"
    unzip -o "${genome}.zip"
done

module load entrez-direct/13.3

for acc in NC_004337.2 NC_003197.2 NC_016901.1 NZ_LN679998.1 NZ_CM003594.1 \
        NZ_CP014946.1 NZ_CP012749.1 NZ_CP011036.1 NZ_OX460909.1 NZ_CP038996.1 \
        NZ_CP059408.1 NZ_CP064354.1 NZ_CP065637.1 NZ_CP066013.1 NZ_CP049902.1 \
        NZ_CP082939.1 NZ_CP090968.1 NZ_CP055055.1 NZ_CP089203.1 NZ_AP024854.1 \
        NZ_CP117568.1 NZ_CP118963.1 NZ_CP146077.1 NZ_CP166085.1 #24
do
  esearch -db nucleotide -query "$acc" | \
  efetch -format fasta > "${acc}.fasta"
done

## improve header for bamplotter readability ##

for fna in GCF*/*.fna; do
  [ -e "$fna" ] || continue
  awk '
  /^>/ {
    n = 0
    for (i = 1; i <= length($0); i++) {
      c = substr($0, i, 1)
      if (c == " " && n < 2) {
        printf "_"
        n++
      } else {
        printf "%s", c
      }
    }
    printf "\n"
    next
  }
  { print }
  ' "$fna" > "${fna%.fna}_clean.fna"
done

#### For efetch data
for fasta in *.fasta; do
  [ -e "$fasta" ] || continue
  awk '
  /^>/ {
    n = 0
    for (i = 1; i <= length($0); i++) {
      c = substr($0, i, 1)
      if (c == " " && n < 2) {
        printf "_"
        n++
      } else {
        printf "%s", c
      }
    }
    printf "\n"
    next
  }
  { print }
  ' "$fasta" > "${fasta%.fasta}_clean.fasta"
done

for fna in GCF*/*_clean.fna; do
  if [ "$(grep -c '^>' "$fna")" -gt 1 ]; then
    echo "$fna"
  fi
done

cat individual_species/ncbi_dataset/data/GCF*/*_clean.fna individual_species/efetch_data/*_clean.fasta >> bacteria_all.fasta
