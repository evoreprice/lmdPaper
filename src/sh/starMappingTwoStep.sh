#!/bin/bash

set -eu

# bash traceback code from https://docwhat.org/tracebacks-in-bash/

_showed_traceback=f

_exit_trap () {
  local _ec="$?"
  if [[ $_ec != 0 && "${_showed_traceback}" != t ]]; then
    traceback 1
  fi
}

_err_trap() {
  local _ec="$?"
  local _cmd="${BASH_COMMAND:-unknown}"
  traceback 1
  _showed_traceback=t
  echo "The command ${_cmd} exited with exit code ${_ec}." 1>&2
}

traceback() {
  # Hide the traceback() call.
  local -i start=$(( ${1:-0} + 1 ))
  local -i end=${#BASH_SOURCE[@]}
  local -i i=0
  local -i j=0

  echo "Traceback (last called is first):" 1>&2
  for ((i=${start}; i < ${end}; i++)); do
    j=$(( $i - 1 ))
    local function="${FUNCNAME[$i]}"
    local file="${BASH_SOURCE[$i]}"
    local line="${BASH_LINENO[$j]}"
    echo "     ${function}() in ${file}:${line}" 1>&2
  done
}

# traps
trap _err_trap SIGHUP SIGINT SIGTERM
trap _exit_trap EXIT
trap _err_trap ERR

# how many CPUs we got?
if [[ $SLURM_JOB_CPUS_PER_NODE ]]; then
	maxCpus="$SLURM_JOB_CPUS_PER_NODE"
	echo -e "[ "$(date)": Running with "$maxCpus" CPUs ]"
else
	maxCpus=1
fi

# handle waiting
FAIL=0
fail_wait() {
for job in $(jobs -p); do
	wait $job || let "FAIL+=1"
done
if [[ ! "$FAIL" == 0 ]]; then
	echo -e "[ "$(date)": Detected fail in background job ]"
	exit 1
fi
}

### CODE STARTS HERE ------------------------------------------------------------------

echo -e "[ "$(date)": Two-step mapping with STAR ]"

# stop if there is no STAR index
star_index_dir="output/star-index"
if [[ ! -d "$star_index_dir" ]]; then
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

# STAR options
OPTIONS="--runThreadN "$maxCpus" --genomeDir "$star_index_dir" --outSAMtype BAM Unsorted --alignIntronMax 5000 --outSJfilterReads Unique --outSJfilterCountUniqueMin 5 5 5 5 --outSJfilterCountTotalMin 5 5 5 5 --outSJfilterIntronMaxVsReadN 5000 --readFilesCommand zcat"

# step 1 function for easy resume
star_step1() {
# load genome
echo -e "[ "$(date)": Loading genome into shared memory ]"
cmd="STAR --runThreadN "$maxCpus" --genomeDir $star_index_dir --genomeLoad LoadAndExit --outFileNamePrefix $outdir/gLoad."
srun --ntasks=1 --exclusive --cpus-per-task="$maxCpus" $cmd

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

# remove step 1 bamfiles (waste of space)
rm "$outdir"/step1/*.out.bam
}

# check if step 1 is already complete
if [[ ! -d "$outdir"/step1 ]]; then
	star_step1
fi

# locate step 1 splice juncions
sjTabs=("$outdir/step1/*.SJ.out.tab")
cat <<- _EOF_
	[ $(date): Found SJ.out.tab files from step 1 ]
	$(for tab in $sjTabs; do echo $tab; done)
_EOF_

### STEP 2

# find the read files
shopt -s nullglob
fastq_files=("$cutadapt_dir/*.fastq.gz")
shopt -u nullglob

echo -e "[ "$(date)": Submitting step 2 mapping jobs ]"

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
