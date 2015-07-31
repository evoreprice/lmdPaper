#!/bin/bash

#SBATCH --job-name conv
#SBATCH --ntasks=6
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --output /tmp/dl.%N.%j.out


outdirTom="data/reads/tomato"

#first decompress
tgzFiles=("$outdirTom/*.tgz")
for tgz in $tgzFiles; do
	bn=$(basename $tgz .tgz)
	cmd="tar --gunzip --get --verbose --to-stdout --file $tgz > $outdirTom/$bn.fastq"
	srun --ntasks=1 --cpus-per-task=1 --exclusive $cmd &
done

wait

# then recompress
find $outdirTom -name *fastq -type f -exec bash -c 'srun --exclusive --ntasks=1 --cpus-per-task=1 gzip --best {} &' \;

wait

# tidy up
find $outdirTom -name "*.fastq" -delete
find $outdirTom -name "*.tgz" -delete

exit 0
