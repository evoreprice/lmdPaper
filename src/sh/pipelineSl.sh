#!/bin/bash

set -eu

_showed_traceback=f

function _exit_trap
{
  local _ec="$?"
  if [[ $_ec != 0 && "${_showed_traceback}" != t ]]; then
    traceback 1
  fi
}

function _err_trap
{
  local _ec="$?"
  local _cmd="${BASH_COMMAND:-unknown}"
  traceback 1
  _showed_traceback=t
  echo "The command ${_cmd} exited with exit code ${_ec}." 1>&2
}

function traceback
{
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

echo -e "[ "$(date)": Starting pipeline for Solanum lycopersicum ]"

# 1. Generate genome for STAR

# define a function for easy testing
genome_generate() {

echo -e "[ "$(date)": Genome generation with STAR ]"

# make output directory
outdir="output/madsComp/sl/star-index"
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# set parameters
genomeFastaFiles="data/genome/sl/Slycopersicum_225_iTAGv2.40.fa"
sjdbGTFfile="data/genome/sl/Slycopersicum_225_iTAGv2.3.gene_exons.cuffcomp.gtf"
sjdbGTFtagExonParentTranscript="oId"
sjdbGTFtagExonParentGene="gene_name"
sjdbOverhang=49
outFileNamePrefix="$outdir/"

# log metadata
cat -t <<- _EOF_ > $outdir/METADATA.csv
	Script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	STAR version,$(STAR --version)
	genomeFastaFiles,$genomeFastaFiles
	sjdbGTFfile,$sjdbGTFfile
	sjdbOverhang,$sjdbOverhang
_EOF_

cat <<- _EOF_
[ $(date): Generating genome ]
              genomeFastaFiles:    $genomeFastaFiles
                   sjdbGTFfile:    $sjdbGTFfile
sjdbGTFtagExonParentTranscript:    $sjdbGTFtagExonParentTranscript
                  sjdbOverhang:    $sjdbOverhang
_EOF_

FAIL=0
cmd="STAR --runThreadN "$maxCpus" --runMode genomeGenerate --genomeDir $outdir --genomeFastaFiles $genomeFastaFiles --sjdbGTFfile $sjdbGTFfile --sjdbGTFtagExonParentTranscript $sjdbGTFtagExonParentTranscript --sjdbGTFtagExonParentGene $sjdbGTFtagExonParentGene --sjdbOverhang $sjdbOverhang --outFileNamePrefix $outFileNamePrefix"
srun --output $outdir/stargg.out --exclusive --ntasks=1 --cpus-per-task="$maxCpus" $cmd &

echo -e "[ "$(date)": Waiting... ]"
fail_wait
echo -e "[ "$(date)": Genome generation finished ]"
}

# check if genome exists already
outdir="output/madsComp/sl/star-index"
if [[ -d $outdir ]]; then
	echo -e "[ "$(date)": Found STAR index ]\n"$outdir""
else
	genome_generate
fi

# 2. Trim adaptors

cutadapt_trim() {

echo -e "[ "$(date)": Trimming reads with cutadapt ]"

# make today's output directory
if [[ ! -d $outdir ]]; then
	mkdir -p $outdir
fi

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
	script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	cutadapt version,$(cutadapt --version)
_EOF_

# parameters
adaptorFwd='TruSeq_indexed_adaptor=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
adaptorRev='TruSeq_universal_adapter=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
trim_qualities=20
minimum_length=50

echo -e "[ "$(date)": Submitting cutadapt jobs ]"

shopt -s nullglob
readFiles=("data/reads/sl/*L.fastq.gz")
shopt -u nullglob
FAIL=0
for fwd_reads in $readFiles; do
	fFile="$(basename $fwd_reads)"
	lib_name="${fFile%"_L.fastq.gz"}"
	rev_reads="data/reads/sl/$(basename $fwd_reads L.fastq.gz)R.fastq.gz"
	if [[ ! -e "$rev_reads" ]]; then echo -e "Error: rev_reads not found\n[ lib_name ]:\t$lib_name\n[ rev_reads ]:\t$rev_reads"; exit 1; fi
	output="$outdir/$lib_name.R1.fastq.gz"
	paired_output="$outdir/$lib_name.R2.fastq.gz"
	# print some info
	cat <<- _EOF_
	[ $(date): Running cutadapt ]
	     lib_name:    $lib_name
	    fwd_reads:    $fwd_reads
	    rev_reads:    $rev_reads
	       output:    $output
	paired_output:    $paired_output
_EOF_
	# run cutadapt
	cmd="cutadapt -a $adaptorFwd -A $adaptorRev --quality-cutoff $trim_qualities --minimum-length $minimum_length --output=$output --paired-output=$paired_output $fwd_reads $rev_reads"
	srun --output $outdir/$lib_name.out --exclusive --ntasks=1 --cpus-per-task=1 $cmd &
done

echo -e "[ "$(date)": Waiting for adaptor trimming to finish ]"
fail_wait
echo -e "[ "$(date)": Adaptor trimming finished ]"
}

# check if trimming has been done
outdir="output/madsComp/sl/cutadapt"
if [[ -d $outdir ]]; then
	trimmedReads=(""$outdir"/*.fastq.gz")
	cat <<- _EOF_
	[ $(date): Found cutadapt output ]
	$(for rf in $trimmedReads; do echo $rf; done)
_EOF_
else
	cutadapt_trim
fi

# 3. Run 2-step mapping

echo -e "[ "$(date)": Two-step mapping with STAR ]"

# stop if there is no STAR index
star_index_dir="output/madsComp/sl/star-index"
if [[ ! -d "$star_index_dir"]]; then
	echo -e "[ "$(date)": No STAR index found ]"
	exit 1
fi
echo -e "[ "$(date)": Using STAR index $star_index_dir ]"

# stop if there is no cutadapt folder
cutadapt_dir="output/madsComp/sl/cutadapt"
if [[ ! -d "$cutadapt_dir" ]]; then
	echo -e "[ "$(date)": No cutadapt folder found ]"
	exit 1
fi
echo -e "[ "$(date)": Using cutadapt folder $cutadapt_dir ]"

# make today's output directory
outdir="output/madsComp/sl/STAR"
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
fastq_files=("$cutadapt_dir/*R1.fastq.gz")
shopt -u nullglob

FAIL=0
for fwd_read_file in $fastq_files; do
	n=$(basename $fwd_read_file)
	library_name="${n%".R1.fastq.gz"}"
	rev_read_file=${fwd_read_file/$library_name.R1/$library_name.R2}
	# double check rev_read_file exists
	if [[ ! -e $rev_read_file ]]; then
		echo -e "[ "$(date)" : Error! Couldn't find reverse read file $rev_read_file for library $library_name ]"
		exit 1
	fi
	cat <<- _EOF_
	[ $(date): Submitting STAR run ]
	 library_name:    $library_name
	fwd_read_file:    $fwd_read_file
	rev_read_file:    $rev_read_file
_EOF_
	cmd="STAR $OPTIONS --genomeLoad LoadAndKeep --readFilesIn $fwd_read_file $rev_read_file --outFileNamePrefix $outdir/step1/$library_name."
	srun --output $outdir/step1/$library_name.out --exclusive --ntasks=1 --cpus-per-task="$maxCpus" $cmd &	
done

echo -e "[ "$(date)": Waiting for step 1 jobs to finish ]"
fail_wait

echo -e "[ "$(date)": Step 1 finished. Removing index from memory ]"
srun --exclusive --ntasks=1 --cpus-per-task="$maxCpus" \
	STAR --runThreadN "$maxCpus" --genomeDir $star_index_dir --genomeLoad Remove --outFileNamePrefix $outdir/gRem.

### STEP 2

echo -e "[ "$(date)": Submitting step 2 mapping jobs ]"

# remove step 1 bamfiles (waste of space)
rm $outdir/step1/*.out.bam

# get the SJ.out.tab files from step 1
sjTabs=("$outdir/step1/*.SJ.out.tab")
cat <<- _EOF_
	[ SJ.out.tab files from step 1 ]
	$(for tab in $sjTabs; do echo $tab; done)
_EOF_

FAIL=0
for fwd_read_file in $fastq_files; do
	n=$(basename $fwd_read_file)
	library_name="${n%".R1.fastq.gz"}"
	rev_read_file=${fwd_read_file/$library_name.R1/$library_name.R2}
	# double check rev_read_file exists
	if [[ ! -e $rev_read_file ]]; then
		echo -e "[ "$(date)" : Error! Couldn't find reverse read file $rev_read_file for library $library_name ]"
		exit 1
	fi
	cat <<- _EOF_
	[ $(date): Submitting step 2 STAR run ]
	 library_name:    $library_name
	fwd_read_file:    $fwd_read_file
	rev_read_file:    $rev_read_file
_EOF_
	cmd="STAR $OPTIONS --genomeLoad NoSharedMemory --sjdbFileChrStartEnd $sjTabs --quantMode GeneCounts --readFilesIn $fwd_read_file $rev_read_file --outFileNamePrefix $outdir/$library_name."
	srun --output $outdir/$library_name.out --exclusive --ntasks=1 --cpus-per-task="$maxCpus" $cmd &	
done

echo -e "[ "$(date)": Waiting for step 2 jobs to finish ]"
fail_wait

echo -e "[ "$(date)": Step 2 jobs finished ]"
echo -e "[ "$(date)": Done, exiting ]"

exit 0
