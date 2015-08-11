#!/bin/bash

set -e

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

# arabidopsis
outdirAt="data/reads/at"
if [[ ! -d $outdirAt ]]; then
	mkdir -p $outdirAt
fi

echo -e "[ "$(date)": Downloading arabidopsis reads ]"
FAIL=0

# IM_R1
wget --output-document "$outdirAt"/IM_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR351/ERR351334/ERR351334.fastq.gz &

# IM_R2
wget --output-document "$outdirAt"/IM_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR351/ERR351333/ERR351333.fastq.gz &

# FM_R1
wget --output-document "$outdirAt"/FM_R1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR351/ERR351336/ERR351336.fastq.gz &

# FM_R2
wget --output-document "$outdirAt"/FM_R2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR351/ERR351331/ERR351331.fastq.gz &

fail_wait

exit 0
