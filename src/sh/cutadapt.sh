#!/bin/bash

set -e

# how many CPUs we got?
if [[ $SLURM_JOB_CPUS_PER_NODE ]]; then
	maxCpus="$SLURM_JOB_CPUS_PER_NODE"
	echo -e "[ "$(date)": Running with "$maxCpus" CPUs ]"
else
	maxCpus=1
fi

# cleanup functions
exit_error() {
	echo -e "[ "$(date)": Script aborted ]"
	exit 1
}

# catch exit codes
trap_exit() {
	exitCode=$?
	if (( "exitCode" == 0 )) ; then
		exit 0
	else
		exit_error
	fi
}

# traps
trap exit_error SIGHUP SIGINT SIGTERM
trap trap_exit EXIT

# handle waiting
FAIL=0
fail_wait() {
for job in $(jobs -p); do
  wait $job || let "FAIL+=1"
done
if [[ ! "$FAIL" == 0 ]]; then
  exit 1
fi
}

### CODE STARTS HERE ------------------------------------------------------------------

# make today's output directory
outdir="output/cutadapt"
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
cat -t <<- _EOF_ > $outdir/METADATA.csv
	script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	cutadapt version,$(cutadapt --version)
_EOF_

echo -e "[ "$(date)": Adaptor trimming with cutadapt ]"

# parameters
adaptor='Illumina_Universal_Adapter=AGATCGGAAGAG'
trim_qualities=20
minimum_length=25

echo -e "[ "$(date)": Submitting cutadapt jobs ]"
readFiles=("data/reads/os/*.fastq.gz")
FAIL=0
for readFile in $readFiles; do
	fFile="$(basename $readFile)"
	lib_name="${fFile%".fastq.gz"}"
	output="$outdir/$lib_name.fastq.gz"
	# print some info
	cat <<- _EOF_
	[ $(date): Submitting cutadapt job ]
	lib_name:  $lib_name
	readFile:  $readFile
	  output:  $output
_EOF_
	# run cutadapt
	cmd="cutadapt -a $adaptor --quality-cutoff $trim_qualities --minimum-length $minimum_length --output=$output $readFile"
	srun --output $outdir/$lib_name.out --exclusive --ntasks=1 --cpus-per-task=1 $cmd &
done

echo -e "[ "$(date)": Waiting for jobs to finish ]"
fail_wait

echo -e "[ "$(date)": Jobs finished, exiting ]"

exit 0
