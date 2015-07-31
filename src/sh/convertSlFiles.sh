#!/bin/bash

#SBATCH --job-name conv
#SBATCH --ntasks 6
#SBATCH --mail-type=ALL

outdirTom="data/reads/tomato"

#first decompress
tgzFiles=("$outdirTom/*.tgz")
for tgz in $tgzFiles; do
	bn=$(basename $tgz .tgz)
	cmd="tar --gunzip --get --verbose --to-stdout --file $tgz > $outdirTom/$bn.fastq"
	echo "srun --ntasks=1 --cpus-per-task=1 --exclusive $cmd &"
done

#tar --gunzip --get --verbose --to-stdout --file data/reads/tomato/wt_fm_r1_R.tgz | gzip --best > data/reads/tomato/wt_fm_r1_R.fastq.gz
