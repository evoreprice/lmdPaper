#!/bin/bash

#SBATCH --job-name cutadapt
#SBATCH --ntasks 6
#SBATCH --output /tmp/cutadapt.%N.%j.out
#SBATCH --open-mode=append

# prep
THEN="$(date)"

# make today's output directory
outdir="output/cutadapt-"$(date +%F)""
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
version="$(cutadapt --version)"
cat -t <<- _EOF_ > $outdir/METADATA.csv
	script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	cutadapt version,$version
_EOF_

# catch SIGTERM etc.
clean_up() {
	# remove temp files before exit
	echo -e "[ "$(date)" : Script aborted ]"
	# email output
	NOW="$(date)"
	MESSAGE="$(cat /tmp/cutadapt."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)"
	cat <<- _EOF_ | mail -s "[Tom@SLURM] Job "$SLURM_JOBID" aborted" tom
	Job $SLURM_JOBID submitted at $THEN was aborted.
	Concatenated stdout files:
	$MESSAGE
_EOF_
	mv /tmp/cutadapt."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out "$outdir"/
	exit 1
}
trap clean_up SIGHUP SIGINT SIGTERM

# parameters
adaptor='Illumina_Universal_Adapter=AGATCGGAAGAG'
trim_qualities=20
minimum_length=25

echo -e "[ "$(date)": Submitting cutadapt jobs ]"
readFiles=("data/reads/*.fastq.gz")
for readFile in $readFiles; do
	fFile="$(basename $readFile)"
	lib_name="${fFile:0:4}"
	output="$outdir/$lib_name.fastq.gz"
	# print some info
	cat <<- _EOF_
	Running cutadapt:
	[ lib_name ]:	$lib_name
	[ readFile ]:	$readFile
	[ output ]:		$output
_EOF_
	# run cutadapt
	cmd="cutadapt -a $adaptor --quality-cutoff $trim_qualities --minimum-length $minimum_length --output=$output $readFile"
	srun --output $outdir/$lib_name.out --exclusive --ntasks=1 --cpus-per-task=1 $cmd &
done

echo -e "[ "$(date)": Waiting for jobs to finish ]"

wait

echo -e "[ "$(date)": Jobs finished, tidying up ]"

# email output
NOW="$(date)"
MESSAGE="$(cat /tmp/cutadapt."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out)"
cat <<- _EOF_ | mail -s "[Tom@SLURM] Job "$SLURM_JOBID" finished" tom
	Job $SLURM_JOBID submitted at $THEN is finished.
	Concatenated stdout files:
	$MESSAGE
_EOF_

mv /tmp/cutadapt."$SLURM_JOB_NODELIST"."$SLURM_JOBID".out "$outdir"/

echo -e "[ "$(date)": Exiting ]"
exit 0
