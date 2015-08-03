#!/bin/bash

#SBATCH --job-name shfcount
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output /tmp/shfcount.%N.%j.out
#SBATCH --mail-type=ALL

# try to find the shuffled GTF
shopt -s nullglob
folders=(output/shuffle*/)
shopt -u nullglob
if (( ${#folders[@]} == 0 )); then
	echo -e "[ "$(date)": Can't find shuffled gff3 ]"
	exit 1
fi
shuffDir="${folders[-1]}"
GTF="$shuffDir"/Osativa_204_v7.shuffled.gff3
if [[ ! -e "$GTF" ]]; then
	echo -e "[ "$(date)": Can't find shuffled gff3 ]"
	exit 1
fi

echo -e "[ "$(date)": Found shuffled gff3 "$GTF" ]"

# make today's output directory
outdir="$shuffDir/htseq-count-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
	Script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	htseq-count version, $(htseq-count --help 2>&1 | tail -n 1)
	gff3,$GTF
	output,$outdir
_EOF_


# find STAR files
# choose the most recent cutadapt output
shopt -s nullglob
folders=(output/cutadapt*)
shopt -u nullglob
# stop if there is no cutadapt folder
if (( ${#folders[@]} == 0 )); then
	echo -e "[ "$(date)": No cutadapt folder found ]"
	exit 1
fi
cutadapt_dir="${folders[-1]}"

shopt -s nullglob
folders=("$cutadapt_dir"/STAR*)
shopt -u nullglob
# stop if there is no star folder
if (( ${#folders[@]} == 0 )); then
	echo -e "[ "$(date)": No star folder found ]"
	exit 1
fi
star_dir="${folders[-1]}"

echo -e "[ "$(date)": Using star folder $star_dir ]"

shopt -s nullglob
bamfiles=(""$star_dir"/*.Aligned.out.bam")
shopt -u nullglob

# run htseq-count

for bam in $bamfiles; do
	n="$(basename $bam)"
	lib_name=${n:0:4}
	cat <<- _EOF_
	[ $(date): Submitting htseq-count job ]
	 bamfile: $bam
	lib_name: $lib_name
	  output: $outdir/$lib_name.htseq-count
_EOF_
	cmd="htseq-count -t CDS -f bam -s no -i ID $bam $GTF"
	srun --ntasks=1 --cpus-per-task=1 --exclusive --output=$outdir/$lib_name.htseq-count --error=$outdir/$lib_name.htseq-count.error $cmd &
done

echo -e "[ "$(date)": Waiting for htseq-count jobs to finish ]"
wait

echo -e "[ "$(date)": Finished, tidying up ]"

# email output
cat <<- _EOF_ | mail -s "[Tom@SLURM] Job "$SLURM_JOBID" finished" tom
	Job $SLURM_JOBID started at $THEN has finished.
	Output:

	$(cat /tmp/shfcount."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)
_EOF_

mv /tmp/shfcount."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out $outdir/shfcount."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out

exit 0
