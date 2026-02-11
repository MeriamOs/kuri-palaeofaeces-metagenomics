#!/bin/bash -e
#SBATCH -A uoo02328
#SBATCH -J krk2
#SBATCH --time 6:00:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -n 1
#SBATCH --mem=90G
#SBATCH --mail-type=ALL
#SBATCH --output krk2_modern.%j.out # CHANGE map1 part each run
#SBATCH --error krk2_modern.%j.err # CHANGE map1 part each run

module purge

module load Kraken2
# note: adding blast module to get dustmasker program to finish
module load BLAST/2.13.0-GCC-11.3.0
## go to directory

cd /nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis

#DATADIR='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/laos_modern/bam_nuclear'
DATADIR2='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/yarlagadda_modern/mapping'
DATADIR3='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/human_modern/trimmed_data'
DATADIR4='/nesi/nobackup/uoo02328/meriam/coprolite_analysis/04-analysis/soil'

## reference
ref='/nesi/nobackup/uoo02328/meriam/kraken/databases/2024-standard'

## run kraken

#### dogs ###
for sample in SRR14842324 SRR14842325 SRR14842327 SRR14842328 SRR14842334 SRR14842335 SRR14842336 SRR14842337 SRR14842338 SRR14842340 \
		SRR14842341 SRR14842342 SRR14842344 SRR14842345 SRR14842346 SRR14842347 SRR14842348 SRR14842350 SRR14842351 SRR14842352 SRR14842363 \
		SRR14842374 SRR14842385 SRR14842396 SRR14842407 SRR14842415 SRR14842424 SRR14842426 SRR14842400 SRR14842401 SRR14842402 SRR14842403 \
		SRR14842404 SRR14842405 SRR14842406 SRR14842408 SRR14842409 SRR14842410 SRR14842411 SRR14842412 SRR14842413 SRR14842414 SRR14842416 \
		SRR14842417 SRR14842418 SRR14842419 SRR14842320 SRR14842420 SRR14842321 SRR14842421 SRR14842322 SRR14842422 SRR14842323 SRR14842423 \
		SRR14842425 SRR14842326 SRR14842329 SRR14842330 SRR14842331 SRR14842332 SRR14842333 SRR14842339 SRR14842343 SRR14842349 SRR14842353 \
		SRR14842354 SRR14842355 SRR14842356 SRR14842357 SRR14842358 SRR14842359 SRR14842360 SRR14842361 SRR14842362 SRR14842364 SRR14842365 \
		SRR14842366 SRR14842367 SRR14842368 SRR14842369 SRR14842370 SRR14842371 SRR14842372 SRR14842373 SRR14842375 SRR14842376 SRR14842377 \
		SRR14842378 SRR14842379 SRR14842380 SRR14842381 SRR14842382 SRR14842383 SRR14842384 SRR14842386 SRR14842387 SRR14842388 SRR14842389 \
		SRR14842390 SRR14842391 SRR14842392 SRR14842393 SRR14842394 SRR14842395 SRR14842397 SRR14842398 SRR14842399;
	do 
	kraken2 --db $ref \
		--threads 32 \
		--classified-out kraken/kraken-2024-standard/reports_yarlagadda/${sample}_kraken_standard2024_conf0.50.fasta \
		--report kraken/kraken-2024-standard/reports_yarlagadda/${sample}_kraken_standard2024_conf0.50.txt \
		--output kraken/kraken-2024-standard/reports_yarlagadda/${sample}_kraken_standard2024_conf0.50.50 \
		--use-names \
		--report-minimizer-data \
		--minimum-base-quality 30 \
		--confidence 0.50 \
		${DATADIR2}/${sample}_unmapped.fq
done

### humans

#SRR1952344 SRR1952455 SRR1952457 SRR1952459 SRR1952521 SRR1952522 SRR1952575 \
#            SRR1952613 SRR1952617 SRR1952618 SRR1952620 SRR1952621 SRR2175654 SRR2175670 \
#            SRR2175725 SRR2175726 SRR2175727 SRR2175750

for sample in SRR2175767 SRR2175768 SRR2175776 SRR2175773 \
        SRR2175777 SRR2175779 SRR2175780 SRR2175804 SRR2240287 SRR2240290 SRR2240728 \
        SRR2240729 SRR2240837 SRR2240921 SRR2241023 SRR2241109 SRR2241305 SRR2241306 \
        SRR2241409 SRR2241505 SRR2175774 SRR2175647 SRR2175649 SRR2175658 SRR2175792 \
		SRR1944873 SRR1952058 SRR1952418 SRR1952501 SRR1952518 SRR1952531 SRR1952562 \
		SRR2175723 SRR2175724 SRR2175749 SRR2175755 SRR2175756 SRR2175764 SRR2175766 ; 
	do 
	kraken2 --db $ref \
		--threads 32 \
		--classified-out kraken/kraken-2024-standard/reports_human/${sample}_kraken_standard2024_conf0.50.fasta \
		--report kraken/kraken-2024-standard/reports_human/${sample}_kraken_standard2024_conf0.50.txt \
		--output kraken/kraken-2024-standard/reports_human/${sample}_kraken_standard2024_conf0.50.50 \
		--use-names \
		--report-minimizer-data \
		--gzip-compressed \
		--minimum-base-quality 30 \
		--confidence 0.50 \
		${DATADIR3}/${sample}_trimmed/${sample}.collapsed.gz
done

### environmental

for sample in ERR12814476 ERR12814474 ERR12814471 ERR12814470 ERR12814467 \
		ERR12074386 ERR12074385 ERR12074384 ERR1883478	ERR1883438	ERR1883421 \
		ERR3761403	ERR3761401	ERR3761400	ERR3761405	SRR7774471	SRR7774470 \
        SRR7774477	SRR7774473	SRR7774472	SRR7774474	SRR7774476	SRR7774469 \
        SRR25620681	SRR25620676	SRR25620756	SRR25620753	SRR25620795	SRR25620790	\
        SRR25620786	ERR3761404	ERR3761402	ERR3761406	ERR10114881	ERR10878162	\
        ERR1883436	ERR1883422	ERR1883479	ERR1883424	ERR1883420	ERR1883477	\
        ERR1883486	
	do
	kraken2 --db $ref \
		--threads 32 \
		--classified-out kraken/kraken-2024-standard/reports_soil/${sample}_kraken_standard2024_conf0.50.fasta \
		--report kraken/kraken-2024-standard/reports_soil/${sample}_kraken_standard2024_conf0.50.txt \
		--output kraken/kraken-2024-standard/reports_soil/${sample}_kraken_standard2024_conf0.50.50 \
		--use-names \
		--report-minimizer-data \
		--gzip-compressed \
		--minimum-base-quality 30 \
		--confidence 0.50 \
		${DATADIR4}/${sample}_trimmed/${sample}.collapsed.gz
done

for sample in SRR24222629 SRR24222627 SRR24222626 SRR24222625 SRR24222623 ERR1560093 ERR1560094 ERR1560095
	do
	kraken2 --db $ref \
		--threads 32 \
		--classified-out kraken/kraken-2024-standard/reports_soil/${sample}_kraken_standard2024_conf0.50.fasta \
		--report kraken/kraken-2024-standard/reports_soil/${sample}_kraken_standard2024_conf0.50.txt \
		--output kraken/kraken-2024-standard/reports_soil/${sample}_kraken_standard2024_conf0.50.50 \
		--use-names \
		--report-minimizer-data \
		--gzip-compressed \
		--minimum-base-quality 30 \
		--confidence 0.50 \
		${DATADIR4}/${sample}_trimmed/${sample}.truncated.gz
done
