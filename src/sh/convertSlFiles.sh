#!/bin/bash

#SBATCH --job-name conv
#SBATCH --ntasks=6
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=1
#SBATCH --output /tmp/dl.%N.%j.out
#SBATCH --error /tmp/dl.%N.%j.out
#SBATCH --open-mode append

outdirTom="data/reads/tomato"

#first decompress
tgzFiles=("$outdirTom/*.tgz")
for tgz in $tgzFiles; do
	bn=$(basename $tgz .tgz)
	cmd="tar --gunzip --get --verbose --to-stdout --file $tgz"
	srun --ntasks=1 --cpus-per-task=1 --output=$outdirTom/$bn.fastq --error=/dev/null --exclusive $cmd &
done

wait

# then recompress
find $outdirTom -name "*fastq" -type f -exec bash -c 'srun --exclusive --ntasks=1 --cpus-per-task=1 gzip --best {} &' \;

wait

# tidy up
find $outdirTom -name "*.tgz" -delete

exit 0
