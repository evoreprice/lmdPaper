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
set +u
if [[ $SLURM_JOB_CPUS_PER_NODE ]]; then
	maxCpus="$SLURM_JOB_CPUS_PER_NODE"
	echo -e "[ "$(date)": Running with "$maxCpus" CPUs ]"
else
	maxCpus=1
fi
set -u

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

outdir="output/quantStats"
if [[ ! -d "$outdir" ]]; then
	mkdir -p "$outdir"
fi

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
	Script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	faidx version,$(faidx --version)
	wgsim version,$(wgsim 2>&1 | sed '3q;d')
	bedtools version,$(bedtools --version)
	tophat2 version,$(tophat2 --version)
	bowtie2-build version,$(bowtie2-build --version | head -n 1)
	samtools version,$(samtools --version | head -n 1)
_EOF_

# download oryza repeat database
echo -e "[ $(date): Downloading plant repeat database ]"
wget --output-document $outdir/TIGR_Oryza_Repeats.v3.3_0_0.fsa ftp://ftp.plantbiology.msu.edu/pub/data/TIGR_Plant_Repeats/TIGR_Oryza_Repeats.v3.3

# find the rRNA sequences headers in the repeat database
echo -e "[ $(date): Finding rRNA genes in plant repeat database ]"
IFS="
"
headers=( $(grep "rRNA" "$outdir"/TIGR_Oryza_Repeats.v3.3_0_0.fsa) )
for header in "${headers[@]}"; do b=$(echo ${header#>}); echo "${b%%\ *}"; done > "$outdir"/rrnaHeaders.tmp

# extract rRNA sequences from repeat database
faidx "$outdir"/TIGR_Oryza_Repeats.v3.3_0_0.fsa $(cat "$outdir"/rrnaHeaders.tmp) > "$outdir"/rRna.fasta

# simulate rRNA reads
echo -e "[ $(date): Simulating rRNA reads from plant repeat database ]"
wgsim -e 0 -1 50 -2 1 -r 0 -R 0 -X 0 -d 0 -s 0 "$outdir"/rRna.fasta $outdir/rRna.wgsim.fq /dev/null

# map rRNA reads
if [[ ! -e "data/genome/os/Osativa_204_v7.0.1.bt2" ]]; then
	echo -e "[ $(date): Need to build a bowtie2 index ]"
	bowtie2-build data/genome/os/Osativa_204_v7.0.fa data/genome/os/Osativa_204_v7.0
fi
echo -e "[ $(date): Mapping simulated reads ]"
srun --ntasks=1 --cpus-per-task="$maxCpus" --exclusive tophat2 --num-threads $maxCpus --output-dir $outdir/tophat data/genome/os/Osativa_204_v7.0 $outdir/rRna.wgsim.fq

# convert tophat accepted hits to bed format
bedtools bamtobed -i $outdir/tophat/accepted_hits.bam > $outdir/accepted_hits.bed6

# download tRNA/rRNA annotation from Rap-DB
echo -e "[ $(date): Downloading Rap-DB annotation ]"
# Rap-DB is down, for now use a symlink in the data/ folder to a version downloaded earlier (22/5/15)
#wget --output-document $outdir/irgsp1_rRNA_tRNA.gff.gz http://rapdb.dna.affrc.go.jp/download/archive/irgsp1_rRNA_tRNA.gff.gz &
cp data/irgsp1_rRNA_tRNA.gff $outdir/irgsp1_rRNA_tRNA.gff

# remove genes with wonky coordinates from the tRNA/rRNA database 
sed '282d' $outdir/irgsp1_rRNA_tRNA.gff | sed '537d' | sed '913d' | gff2bed > $outdir/irgsp1_rRNA_tRNA.bed

# extract rRNA regions from rap-db annotation and fix chromosome names
echo -e "[ $(date): Extracting rRNA regions and converting to BED6 ]"
grep "rRNA" $outdir/irgsp1_rRNA_tRNA.bed | sed -e 's/chr/Chr/' | sed -e 's/Chr0/Chr/' > "$outdir"/rap_rRNA.bed9
cat "$outdir"/rap_rRNA.bed9 | cut -f1-6 > "$outdir"/rap_rRNA.bed6

# extract tRNA regions from rap-db annotation and fix chromosome names
grep "tRNA" $outdir/irgsp1_rRNA_tRNA.bed | sed -e 's/chr/Chr/' | sed -e 's/Chr0/Chr/' > "$outdir"/rap_tRNA.bed9

# combine mapped rRNA reads with rap-db regions
sort-bed "$outdir"/rap_rRNA.bed6 "$outdir"/accepted_hits.bed6 | bedtools merge -i stdin > $outdir/rRna.combined.bed

# count rRNA and tRNA, write to temporary files
bamfiles=("output/STAR/*.Aligned.out.bam")
for bam in $bamfiles; do
	library=$(basename "$bam" ".Aligned.out.bam")
	echo -e "[ $(date): Counting rRNA and tRNA for "$library" ]"
	cmdR="samtools view -c -F0x100 -L "$outdir"/rRna.combined.bed "$bam""
	srun --ntasks=1 --exclusive --output="$outdir"/"$library".rrna.tmp $cmdR &
	cmdT="samtools view -c -F0x100 -L "$outdir"/rap_tRNA.bed9 "$bam""
	srun --ntasks=1 --exclusive --output="$outdir"/"$library".trna.tmp $cmdT &
done
echo -e "[ $(date): Waiting for jobs to finish ]"

fail_wait

# read results from temporary fines into array
lines=()
for bam in $bamfiles; do
	library=$(basename "$bam" ".Aligned.out.bam")
	rrna=$(cat $outdir/$library.rrna.tmp)
	trna=$(cat $outdir/$library.trna.tmp)
	lines+=("$library","$rrna","$trna")
	rm "$outdir"/"$library".rrna.tmp
	rm "$outdir"/"$library".trna.tmp
done

# output results
echo -e "[ $(date): Printing results to "$outdir"/quantStats.csv ]"
cat <<- _EOF_ > $outdir/quantStats.csv
	library,rRNA,tRNA
	$(for line in "${lines[@]}"; do echo "$line"; done)
_EOF_

echo -e "[ $(date): Done ]"
exit 0
