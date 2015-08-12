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

echo -e "[ "$(date)": Two-step mapping with STAR ]"

# stop if there is no STAR index
star_index_dir="output/star-index"
if [[ ! -d "$star_index_dir"]]; then
	echo -e "[ "$(date)": No STAR index found ]"
	exit 1
fi
echo -e "[ "$(date)": Using STAR index $star_index_dir ]"

# stop if there is no cutadapt folder
cutadapt_dir="output/cutadapt"
if [[ ! -d "$cutadapt_dir" ]]; then
	echo -e "[ "$(date)": No cutadapt folder found ]"
	exit 1
fi
echo -e "[ "$(date)": Using cutadapt folder $cutadapt_dir ]"

# make today's output directory
outdir="output/STAR"
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
	Script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	STAR version,$(STAR --version)
	STAR index,$star_index_dir
	cutadapt folder,$cutadapt_dir
_EOF_

# load genome
echo -e "[ "$(date)": Loading genome into shared memory ]"
cmd="STAR --runThreadN "$maxCpus" --genomeDir $star_index_dir --genomeLoad LoadAndExit --outFileNamePrefix $outdir/gLoad."
srun --ntasks=1 --exclusive --cpus-per-task="$maxCpus" $cmd

# STAR options
OPTIONS="--runThreadN "$maxCpus" --genomeDir "$star_index_dir" --outSAMtype BAM Unsorted --alignIntronMax 5000 --outSJfilterReads Unique --outSJfilterCountUniqueMin 5 5 5 5 --outSJfilterCountTotalMin 5 5 5 5 --outSJfilterIntronMaxVsReadN 5000 --readFilesCommand zcat"

### STEP 1

# make directory for step 1
if [[ ! -d $outdir/step1 ]]; then
	mkdir -p $outdir/step1
fi

echo -e "[ "$(date)": Submitting step 1 mapping jobs ]"

# find the read files
shopt -s nullglob
fastq_files=("$cutadapt_dir/*.fastq.gz")
shopt -u nullglob

FAIL=0
for read_file in $fastq_files; do
	n=$(basename $read_file)
	library_name="${n%".fastq.gz"}"
	cat <<- _EOF_
	[ $(date): Submitting STAR run ]
	library_name:   $library_name
	read_file:      $read_file	
_EOF_
	cmd="STAR $OPTIONS --genomeLoad LoadAndKeep --readFilesIn $read_file --outFileNamePrefix $outdir/step1/$library_name."
	srun --output $outdir/step1/$library_name.out --exclusive --ntasks=1 --cpus-per-task="$maxCpus" $cmd &	
done

echo -e "[ "$(date)": Waiting for step 1 jobs to finish ]"
fail_wait

echo -e "[ "$(date)": Step 1 finished. Removing index from memory ]"
srun --exclusive --ntasks=1 --cpus-per-task="$maxCpus" \
	STAR --runThreadN "$maxCpus" --genomeDir $star_index_dir --genomeLoad Remove --outFileNamePrefix $outdir/gRem.

### STEP 2

echo -e "[ "$(date)": Submitting step 2 mapping jobs ]"

# remove step 1 bamfile (waste of space)
rm "$outdir"/step1/*.out.bam

# get the SJ.out.tab files from step 1
sjTabs=("$outdir/step1/*.SJ.out.tab")
cat <<- _EOF_
	[ SJ.out.tab files from step 1 ]
	$(for tab in $sjTabs; do echo $tab; done)
_EOF_

FAIL=0
for read_file in $fastq_files
do
	n=$(basename $read_file)
	library_name="${n%".fastq.gz"}"
	cat <<- _EOF_
	[ $(date): Submitting step 2 STAR run ]
	library_name:   $library_name
	read_file:      $read_file	
_EOF_
	cmd="STAR $OPTIONS --genomeLoad NoSharedMemory --sjdbFileChrStartEnd $sjTabs --quantMode GeneCounts --readFilesIn $read_file --outFileNamePrefix $outdir/$library_name."
	srun --output $outdir/$library_name.out --exclusive --ntasks=1 --cpus-per-task="$maxCpus" $cmd &	
done

echo -e "[ "$(date)": Waiting for step 2 jobs to finish ]"
fail_wait
echo -e "[ "$(date)": Jobs finished. Exiting ]"
exit 0
