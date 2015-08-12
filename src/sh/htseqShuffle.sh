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

# check for shuffled GTF

shuffDir="output/shuffle"
if [[ ! -d "$shuffDir" ]]; then
	echo -e "[ "$(date)": Can't find shuffled gff3 ]"
	exit 1
fi
GTF="$shuffDir"/Osativa_204_v7.shuffled.gff3
if [[ ! -e "$GTF" ]]; then
	echo -e "[ "$(date)": Can't find shuffled gff3 ]"
	exit 1
fi

echo -e "[ "$(date)": Found shuffled gff3 "$GTF" ]"

# make today's output directory
outdir="output/dnaTpm/shuffledCounts"
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

# stop if there is no STAR output
star_dir="output/STAR"
if [[ ! -d "$star_dir" ]]; then
	echo -e "[ "$(date)": No STAR outputfound ]"
	exit 1
fi
echo -e "[ "$(date)": Using star folder $star_dir ]"

shopt -s nullglob
bamfiles=(""$star_dir"/*.Aligned.out.bam")
shopt -u nullglob

# run htseq-count
FAIL=0
for bam in $bamfiles; do
	n="$(basename $bam)"
	lib_name="${n%".Aligned.out.bam"}"
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
fail_wait

echo -e "[ "$(date)": Finished, exiting ]"
exit 0
