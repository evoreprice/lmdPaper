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

### CODE STARTS HERE ------------------------------------------------------------------

# stop if there is no madsPeptides
outdir="output/madsComp/clustal"
if [[ ! -d "$outdir" ]]; then
  echo -e "[ "$(date)": outdir $outdir not found ]"
  exit 1
fi
if [[ ! -e "$outdir"/madsPeptides.fasta ]]; then
  echo -e "[ "$(date)": madsPeptides not found ]"
  exit 1
fi

# log metadata
cat <<- _EOF_ > $outdir/METADATA.csv
  Script,${0}
  branch,$(git rev-parse --abbrev-ref HEAD)
  hash,$(git rev-parse HEAD)
  date,$(date +%F)
  clustal-omega version,$(clustalo --version)
_EOF_

# align

clustalo -i "$outdir"/madsPeptides.fasta --full --force --outfmt=clustal --outfile="$outdir"/madsPeptides.aln

exit 0
