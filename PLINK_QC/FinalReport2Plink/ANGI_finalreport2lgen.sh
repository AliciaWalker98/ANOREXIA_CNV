#!/bin/bash

set -eu

# convert individual finalreports to lgen, then plink binary

REPORTDIR=$1
OUTPREFIX=$2
PLINKFILEDIR=/cifs/P/Pgced_cnvs/PLINK_GENO

find $REPORTDIR -name '*FinalReport[0-9]*.txt.bz2' -print0 | xargs -0 bzcat |
    # drop Illumina multi-line header
    sed '/\[Header\]/,/\[Data\]/d' |
    awk 'BEGIN {FS="\t";RS="\r\n"}
    	# column headers
    	/Sample ID/ {
	    if (! snpcol) {firstfile = 1} else {firstfile = 0}
	    for (i = 1; i <= NF; i++) {
	        if ($(i) == "SNP Name") snpcol = i
		if ($(i) == "Sample ID") samplecol = i
		if ($(i) == "Chr") chrcol = i
		if ($(i) == "Position") poscol = i
		if ($(i) == "Allele1 - Forward") a1col = i
		if ($(i) == "Allele2 - Forward") a2col = i
	    }
	    print "Processing file number " ++fno > "/dev/stderr"
	    if (! (snpcol && samplecol && chrcol && poscol && a1col && a2col)) {print "All expected columns not found"; exit}
	}
	# data lines
	!/Sample ID/ {
	     # map file record
	     if (firstfile) {
		 print $(chrcol), $(snpcol), 0, $(poscol) > "'$TMPDIR/$OUTPREFIX'" ".map"
	     }
	     # fam file record
	     if (lastsample != $(samplecol)) {
	         print $(samplecol), $(samplecol), 0, 0, 0, -9 > "'$TMPDIR/$OUTPREFIX'" ".fam"
	     }
	     lastsample = $(samplecol)
	     # genotypes (lgen)
	     print  $(samplecol), $(samplecol), $(snpcol), ($(a1col) == "-")? 0 : $(a1col), ($(a2col) == "-") ? 0 : $(a2col) > "'$TMPDIR/$OUTPREFIX'" ".lgen"
	}
	# wrap up
	END {
	}'

mkdir -p $PLINKFILEDIR/$OUTPREFIX
plink --threads 10 --lfile $TMPDIR/$OUTPREFIX --reference $PLINKFILEDIR/GSA-24v1-0_C1.fwd.ref --make-bed --out $PLINKFILEDIR/$OUTPREFIX/$OUTPREFIX
rm $TMPDIR/$OUTPREFIX.{fam,map,lgen}
