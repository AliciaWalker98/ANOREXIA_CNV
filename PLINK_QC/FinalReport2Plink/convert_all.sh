#!/bin/bash

set -eu

# convert all batches of FinalReports to Plink binary format
# NOTE: meant to Slurm this, but cannot access /cifs from compute nodes

for REPORTDIR in \
    /cifs/P/Pgced_cnvs/ANGI_Bulik_1v3_190430/PGCED_SE_ANGI_Bulik_20190506/Reports_ANGI_Bulik_20190506 \
	/cifs/P/Pgced_cnvs/ANGI_Kaye_190430/PGCED_SE_ANGI_Kaye_190430/Reports_PGCED_SE_ANGI_Kaye_190430 \
	/cifs/P/Pgced_cnvs/ANGI_Landen-1/PGCED_SE_Landen_AN_Sample1_20190312/PGCED_SE_Landen_AN_Sample1/PGCED_SE_Landen_AN_Sample1/Reports_PGCED_SE_Landen_AN_Sample1 \
	/cifs/P/Pgced_cnvs/ANGI_Landen-2/PGCED_SE_Landen_AN_Sample2_20190318/PGCED_SE_Landen_AN_Sample2_Reports \
	/cifs/P/Pgced_cnvs/ANGI_Martin_v1_3_190429/PGCED_SE_ANGI_Martin_v1_3_190515/Reports_ANGI_Martin_v1_3_190515 \
	/cifs/P/Pgced_cnvs/ANGI_PhaseII_Pedersen_Anorexia_GSA-MD_wave1/PGCED_SE_Pedersen_AN_20190118/FinalPGCED_SE_Pedersen_AN_20190118/Reports_SE_Pedersen_AN_20190118 ; do
    OUTPREFIX=$(dirname $REPORTDIR | xargs basename | sed 's/^Final//')

    #sbatch -p core -n 10 -J $OUTPREFIX -t 24:00:00 ./ANGI_finalreport2lgen.sh $REPORTDIR $OUTPREFIX
    echo Starting $OUTPREFIX
    ./ANGI_finalreport2lgen.sh $REPORTDIR $OUTPREFIX
done
