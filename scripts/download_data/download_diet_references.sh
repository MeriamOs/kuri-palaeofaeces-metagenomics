# Script to download diet screening references and prepare for mapping
# Created by Meriam van Os

mkdir mitochondrion
mkdir COI
mkdir plants

module load entrez-direct/13.3

#####################
### MITOCHONDRION ###
#####################

cd mitochondrion

for acc in NC_002008.4	NC_004023.1	NC_008418.1	NC_002503.2	NC_005272.1	NC_005270.1	NC_023830.1	KC312626.1 \
            NC_006930.1	NC_012389.1	CM019797.1	NC_026899.1	NC_088427.1	NC_002813.1	NC_066972.1	NC_036931.1 \
            NC_018045.1	NC_020046.1	NC_024236.1	NC_008449.1	NC_022508.1	NC_085253.1	NC_057117.1	NC_035059.1 \
            NC_017879.1	JQ639063.1	CM100769.1	NC_015787.1	NC_023527.1	MF278960.1	NC_002672.1	MK294166.1 \
            NC_013806.1	NC_027267.1	NC_025917.1	NC_057528.1	NC_052809.1	NC_004538.1	MK290247.1	KC875856.1 \
            NC_029144.1	MK761003.1	MK290258.1	NC_023952.1	NC_003190.1	GU071054.1	NC_027845.1	PP929789.1 \
            NC_013244.1	NC_029141.1	NC_007006.1	MK760984.1	MK760983.1	NC_043899.1	GU810158.1	NC_059921.1 \
            NC_037185.1	NC_036483.1	MK577964.1	NC_068097.1	NC_013725.1	NC_003182.1	NC_023448.1	GU071056.1 \
            PP929796.1	NC_005932.1	NC_035423.1	NC_031361.1	ON990049.1	NC_033393.1	NC_088542.1	NC_002012.1\
            EU398541.1	PP533446.1; #74
        do
        esearch -db nucleotide -query "$acc" | \
        efetch -format fasta > "${acc}.fasta"
done

## improve header for later ##

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

##################
### COI REGION ###
##################

cd ../COI

for acc in HM422251.1	GU806141.1	MN123447.1	MH707750.1	EF609337.1	EF609336.1	DQ107893.1	DQ108054.1 \
          GU806164.1	HQ956258.1	KX781853.1	KP267655.1	KJ012398.1	JF494694.1	EF609428.1	EF609444.1 \
          MK261845.1	JN801370.1	MK261804.1	MK261953.1	MK261898.1	MK262152.1	AB588860.1	HM887971.1 \
          JN200827.1	AY858100.1	AF546039.1	HM180704.1	GU324171.1	HM888033.1	AB433649.1	GQ387270.1 \
          MK261882.1	MK261892.1	EF101681.2	MK261782.1	DQ108051.1	AY858094.1	AM049374.1	OR910872.1 \
          MH485162.1	EU391573.1	DQ916503.1	KC669662.1	HQ956117.1	DQ107900.1	MK101227.1	MK262546.1 \
          JQ175767.1	JN801453.1	MK262500.1	DQ837002.1	AM403873.1	KP980809.1	MK016469.1	MW039509.1 \
          GQ249699.1	DQ298507.1	EU265768.1	DQ917614.1	KP844998.1	KT883635.1	AF156244.1	GU679400.1; #64
        do
        esearch -db nucleotide -query "$acc" | \
        efetch -format fasta > "${acc}.fasta"
done

## improve header for later ##

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

##############
### PLANTS ###
##############

cd ../plants

for acc in NC_036808.1	OM808940.1	NC_035998.1	OP168865.1	NC_070422.1	NC_000932.1	PP100725.1	NC_045119.1 \
          MZ662865.1	NC_016753.1	NC_014807.1	NC_039707.1	MW470977.1	NC_062424.1	MT621036.1	NC_050985.1 \
          MW191857.1	MW191854.1	MZ144895.1	MZ144929.1	MW218462.1	MW218464.1	MW218459.1	MT106774.1 \
          NC_050678.1	AP018905.1	NC_034942.1	NC_020361.1	OR074083.1	NC_063825.1	OQ848587.1	NC_036991.1 \
          MW214669.1	MW214671.1	PV409670.1	OL938727.1	MW528270.1	MW491882.1	MN239902.1	NC_050955.1 \
          MT385084.1	NC_051904.1	NC_036681.1	NC_047178.1	MT385079.1	MW191877.1	MW221973.1; #47
        do
        esearch -db nucleotide -query "$acc" | \
        efetch -format fasta > "${acc}.fasta"
done

## improve header for later ##

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

cd ..

cat mitochondrion/*_clean.fasta >> diet_screen.fasta
cat COI/*_clean.fasta >> COI_diet_screen.fasta
cat plants/*_clean.fasta >> plants.fasta
