#!/bin/bash

# catch email address (-e) and password (-p) for jgi passed to bash script
# if the jgi password has special characters this script won't work

while [ "$1" != "" ]; do
	case $1 in
		-e )	shift
				EMAIL=$1
				;;
		-p )	shift
				PASSWORD=$1
				;;
		* )		echo "Bad input"
				exit 1
	esac
	shift
done

genome_dir='data/genome'
genome_url='http://genome.jgi.doe.gov/Osativa/download/_JAMO/53112abc49607a1be00559bc/Osativa_204_v7.0.fa.gz'
annot_url='http://genome.jgi.doe.gov/Osativa/download/_JAMO/53112ab649607a1be00559b0/Osativa_204_v7.0.gene_exons.gff3.gz'

genome_file="$(basename $genome_url .fa.gz)"
annotation_file="$(basename $annot_url .gff3.gz)"

# download genomes
if [ ! -d $genome_dir ]; then
	mkdir -p $genome_dir
fi

curl https://signon.jgi.doe.gov/signon/create --data-ascii \
	login="$EMAIL"\&password="$PASSWORD" \
	-b cookies -c cookies > /dev/null

curl $genome_url -b cookies -c cookies > $genome_dir/$genome_file.fa.gz
curl $annot_url -b cookies -c cookies > $genome_dir/$annotation_file.gff3.gz

gunzip $genome_dir/$genome_file.fa.gz
gunzip $genome_dir/$annotation_file.gff3.gz

rm cookies

# make cuffcomp gtf
cuffcompare -s $genome_dir/$genome_file.fa -CG -r $genome_dir/$annotation_file.gff3 \
	-o $genome_dir/$annotation_file.cuffcomp $genome_dir/$annotation_file.gff3

# remove Chr9 rRNA 'genes'
sed '/LOC_Os09g01000/d' $genome_dir/$annotation_file.cuffcomp.combined.gtf \
	| sed '/LOC_Os09g00999/d' \
	> $genome_dir/gtf_final.tmp

# remove cuffcomp intermediates
rm $genome_dir/*cuffcomp*
mv $genome_dir/gtf_final.tmp $genome_dir/$annotation_file.cuffcomp.rRNAremoved.gtf

# log metadata
cat -t <<- _EOF_ > $genome_dir/METADATA.csv
	script,${0}
	branch,$(git rev-parse --abbrev-ref HEAD)
	hash,$(git rev-parse HEAD)
	date,$(date +%F)
	cuffcomp version,$(cuffcompare -h 2>&1 | head -n 1)
_EOF_

exit 0
