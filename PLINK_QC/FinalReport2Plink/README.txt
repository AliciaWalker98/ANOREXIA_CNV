README for folder P:\Pgced_cnvs\PLINK_GENO

created by Robert Karlsson <robert.karlsson@ki.se> on 2022-06-01

This folder contains plink binary format genotypes extracted and
converted from Illumina FinalReport files for each batch of ANGI
samples.

The following scripts were used to generate the files:

ANGI_manifest2ref.R - prepare .ref file with allele codes for each SNP
ANGI_finalreport2lgen.sh - convert a single folder of .txt.bz2 FinalReports
			 to plink binary format
convert_all.sh - run the above conversion for all batches

Note that the generated genotype files have only one ID (labelled
"chip_well_barcode" in metadata files), and that sex and case/control
status fields of genotypes are not filled in.

See metadata files for each batch (in its subfolder in the main
project directory) for additional sample IDs and reported sex.
