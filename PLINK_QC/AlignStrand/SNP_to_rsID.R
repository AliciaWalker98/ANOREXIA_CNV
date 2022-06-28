# Julia Sidorenko
# 04/07/19
# This script uses dbSNP build 144 (dbSNP144.GRCh37) to rename genotyped SNPs on Illumina array to rsIDs

##PLEASE NOTE 1: the small insertions and deletions (I/D) are left unchanged in this pipeline 
##PLEASE NOTE 2: If the alleles for insertions and deletions (>1bp change) are other than "I" / "D" ,
# they won't be filtered out and can be replaced by SNPs (1bp change) at that specific position 

library(plyr)
#BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

ref <- SNPlocs.Hsapiens.dbSNP144.GRCh37
snpcount(ref)
seqlevels(ref)


# Read in bim file
bimname = commandArgs(trailingOnly = TRUE)
bimfile <- read.table(bimname, colClasses=c("numeric", "character", "numeric", "numeric", "character", "character"))

#separte the insertions/deletions
d<-subset(bimfile,V5=="D" | V5=="I" | V6=="D" | V6=="I" )
bim<-subset(bimfile,V5!="D" & V5!="I" & V6!="D" & V6!="I" )
bim$index <- 1:nrow(bim)
head(bim)
tail(bim)

## change 23 and 25 to X
## CAUTION: if there are other chr than 1~23,25, remove them use plink. Don't remove them only in bim file.
bim$V1<-ifelse(bim$V1=="23", "X",bim$V1)
bim$V1<-ifelse(bim$V1=="25", "X",bim$V1)
bim$V1 <- ifelse(bim$V1=="24", "Y",bim$V1)
bim$V1 <- ifelse(bim$V1=="26", "MT",bim$V1)

# For each chromosome download the SNP locations and match to bim file
a <- ddply(bim, .(V1), .progress="text", function(x)
{
  x <- mutate(x)
  chr <- as.character(x$V1[1])
  snps <- data.frame(snpsBySeqname(ref, chr))
  
  snps <- subset(snps, pos %in% x$V4, select=-c(alleles_as_ambig, strand))
  snps <- subset(snps, !duplicated(pos))
  snps <- subset(snps, !duplicated(RefSNP_id))
  x <- merge(x, snps, by.x="V4", by.y="pos", all.x=T)
  x <- x[order(x$index), ]
  index <- !is.na(x$RefSNP_id)
  x$V2[index] <- x$RefSNP_id[index]
  x <- subset(x, select=c(V1, V2, V3, V4, V5, V6))
  return(x)
})

# If there are duplicated SNPs for any reason then label them as _2, _3 etc
temp <- rle(a$V2)
temp2 <- paste0(rep(temp$values, times = temp$lengths), "_", unlist(lapply(temp$lengths, seq_len)))
temp2 <- gsub("_1", "", temp2)
a$V2 <- temp2


# Rename X and XY to 23
a[a[,1]=="X",1]=23
a[a[,1]=="Y",1]=24
a[a[,1]=="MT",1]=26


# bind the insertions/deletions back
c<-rbind (a,d)

# add PAR boundaries boundaries 'b37'/'hg19': GRCh37/UCSC human genome 19, boundaries 2699520 and 154931044 
c$V1<-ifelse(c$V1==23 & c$V4<2699520, "25", c$V1)
c$V1<-ifelse(c$V1==23 & c$V4>154931044, "25", c$V1)

# reorder the SNPs
sorted.2=c[with(c, order(as.numeric(V1),V4)),]

#check the order
length(which(bimfile$V4!=sorted.2$V4))

# Save file
write.table(sorted.2, file=paste(bimname, "_rsed_dbSNP144.GRCh37", sep="" ), row=F, col=F, qu=F)
