#!/usr/bin/env Rscript

# Create ".ref" file from Illumina manifest to get all allele codes with lgen format
# (forward strand)

library(data.table)
library(stringr)

flip <- function(x) c("A"="T", "C"="G", "G"="C", "T"="A", "D"="D", "I"="I")[x]

m <- fread("sed '0,/\\[Assay\\]/d; /\\[Controls\\]/,$d' GSA-24v1-0_C1.csv")
# (T/B/M/P)_(F/R) in IlmnID indicate whether the SNP field alleles are forward 
# or reverse coded - if reverse, flip to get forward
m[, tfbr := str_match(IlmnID, "_([TBMP]_[FR])_")[,-1]]
m[, c("IlmnA1", "IlmnA2") := as.data.frame(str_match(SNP, "\\[(.)/(.)\\]")[, -1])]
# write .ref format
fwrite(m[, list(Name,
                ifelse(grepl("R$", tfbr), flip(IlmnA1), IlmnA1),
                ifelse(grepl("R$", tfbr), flip(IlmnA2), IlmnA2))],
       "GSA-24v1-0_C1.fwd.ref",
       quote = FALSE, sep = " ", col.names = FALSE)


# code below was used to check that flipping made sense
if (FALSE) { 
    m <- fread('~/ANGI_CNV/GSA-24v1-0_C1.csv', skip=7, nrows=618540)
    bim <- fread('/scratch/tmp/robkar/ANGI_Bulik_20190506_PLINKB.bim', select=c(1,2,4:6), header=FALSE, col.names = c("CHR", "SNP", "POS", "A1", "A2"))
    mrg <- merge(bim, m[, list(IlmnID, Name, IlmnSNP=SNP, IlmnStrand, SourceStrand, RefStrand)], by.x="SNP", by.y="Name")
    
    mrg[, c("IlmnA1", "IlmnA2") := as.data.frame(str_match(IlmnSNP, "\\[(.)/(.)\\]")[, -1])]
    mrg[, sameA := (A1 == "0" | A1 == IlmnA1 | A1 == IlmnA2) & (A2 == "0" | A2 == IlmnA1 | A2 == IlmnA2)]
    
    mrg[, tfbr := str_match(IlmnID, "([TBMP]_[FR])")[,-1]]
    
    mrg[, IlmnA1F := if_else(grepl("R$", tfbr), flip(IlmnA1), IlmnA1)]
    mrg[, IlmnA2F := if_else(grepl("R$", tfbr), flip(IlmnA2), IlmnA2)]
    mrg[, sameAF := (A1 == "0" | A1 == IlmnA1F | A1 == IlmnA2F) & (A2 == "0" | A2 == IlmnA1F | A2 == IlmnA2F)]
}    
