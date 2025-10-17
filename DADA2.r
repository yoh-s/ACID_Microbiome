#DADA 2 workflow
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # set current dir as WD
# setwd("~/_PROJETS/Documents/projectethiopia")


# -------------------------------------------------------------
# 1. Packages Imports
# -------------------------------------------------------------
library(dada2); packageVersion("dada2")
library(tictoc)

# -------------------------------------------------------------
# 2. Config
# -------------------------------------------------------------
gene="16S"  # ITS, 16S ou EF1
source("__config.R")

#follow the following step if cutadapt is not necessary
##We split the data between forward (fnFs) and reverse (fnRs) sequences
fnFs <- sort(list.files(rawdata, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(rawdata, pattern="_R2.fastq.gz", full.names = TRUE))

##extract the sample names from the forward reads
sample.names5 <- sapply(strsplit(basename(fnFs), "-"), `[`, 5)
sample.names6 <- sapply(strsplit(basename(fnFs), "-"), `[`, 6)
sample.names6_n <- sapply(strsplit(basename(sample.names6), "_"), `[`, 1)

#sample names
sample.names = paste(sample.names5, sample.names6_n, sep = "-")

##Display of graphs representing the quality profiles 
plotQualityProfile(fnFs[1:3]) 
plotQualityProfile(fnRs[1:3])

##Place the filtered sequences in a 'Filtered' sub-directory
filtFs <- file.path(data_folder, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(data_folder, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

##Filter and trim
if (file.exists(file.path(output_results,"out_maxee2.rds"))) {
  out_maxee2=readRDS(file.path(output_results,"out_maxee2.rds"))
} else {
  out_maxee2 <- filterAndTrim(fnFs, filtFs, 
                              fnRs, filtRs, 
                              truncLen=c(266,261),
                              maxN=0, maxEE=c(2,2), 
                              truncQ=2, rm.phix=TRUE, 
                              compress=TRUE, 
                              multithread=TRUE, verbose=TRUE) # On Windows set multithread=FALSE
  
  saveRDS(out_maxee2, file.path(output_results,"out_maxee2.rds"))
}

head(out_maxee2)

filt_trim_maxee2 <- tibble(SampleIDs = rownames(out_maxee2), Input = out_maxee2[,1], Filtered = out_maxee2[,2])
#display the # of reads pre- and post-filtering
filt_trim_maxee2

write.csv(filt_trim_maxee2,  "filt_trim_maxee2.csv")

------------------------------------------------------------------------------
# Learn the error rates from maxee2 (MAXEE 2,2)
filtFs 
filtRs
# forward reads error

set.seed(1988)
if (file.exists(file.path(output_results,"errF_maxee2.rds"))) {
  errF_maxee2=readRDS(file.path(output_results,"errF_maxee2.rds"))
} else {
  errF_maxee2 <- learnErrors(filtFs, randomize = TRUE, multithread=TRUE)
  saveRDS(errF_maxee2,file.path(output_results,"errF_maxee2.rds"))
}

# reverse reads error

if (file.exists(file.path(output_results,"errR_maxee2.rds"))) {
  errR_maxee2=readRDS(file.path(output_results,"errR_maxee2.rds"))
} else {
  errR_maxee2 <- learnErrors(filtRs, randomize = TRUE, multithread=TRUE)
  saveRDS(errR_maxee2,file.path(output_results,"errR_maxee2.rds"))
}
------------------------------------------------------------------------------
# Dereplication of the reads (MAXEE 2,2)
#filtFs 
#filtRs

#FORWARD
if (file.exists(file.path(output_results, "derepFs_maxee2.rds"))) {
  derepFs_maxee2=readRDS(file.path(output_results,"derepFs_maxee2.rds"))
} else {
  derepFs_maxee2 <- derepFastq(filtFs, verbose = TRUE)
  saveRDS(derepFs_maxee2,file.path(output_results,"derepFs_maxee2.rds"))
}

#REVERSE
if (file.exists(file.path(output_results, "derepRs_maxee2.rds"))) {
  derepRs_maxee2=readRDS(file.path(output_results,"derepRs_maxee2.rds"))
} else {
  derepRs_maxee2 <- derepFastq(filtRs, verbose = TRUE)
  saveRDS(derepRs_maxee2,file.path(output_results,"derepRs_maxee2.rds"))
}


# Name the derep-class objects by the sample names ##should run
names(derepFs_maxee2) <- sample.names
names(derepRs_maxee2) <- sample.names

--------------------------------------------------------------------------------
## PSEUDO POOLING
## Inference of Amplicon Sequence Variants - truncLen = c(266, 261)
# maxEE(2, 2)
derepFs_maxee2
derepRs_maxee2

if (file.exists(file.path(output_results, "dadaF_maxee2.rds"))) {
  dadaF_maxee2=readRDS(file.path(output_results,"dadaF_maxee2.rds"))
} else {
  dadaF_maxee2 <- dada(derepFs_maxee2, err = errF_maxee2, pool = "pseudo", multithread = TRUE)
  saveRDS(dadaF_maxee2,file.path(output_results,"dadaF_maxee2.rds"))
}

#REVERSE

if (file.exists(file.path(output_results, "dadaR_maxee2.rds"))) {
  dadaR_maxee2=readRDS(file.path(output_results,"dadaR_maxee2.rds"))
} else {
  dadaR_maxee2 <- dada(derepRs_maxee2, err = errR_maxee2, pool = "pseudo", multithread = TRUE)
  saveRDS(dadaR_maxee2,file.path(output_results,"dadaR_maxee2.rds"))
}
## 

--------------------------------------------------------------------------------
#Merging of the forward and reverse reads
if (file.exists(file.path(output_results,"maxeeon_mergers_pseudo.rds"))) {
  maxeeon_mergers_pseudo=readRDS(file.path(output_results,"maxeeon_mergers_pseudo.rds"))
} else {
  maxeeon_mergers_pseudo <- mergePairs(dadaF_maxee2, derepFs_maxee2, dadaR_maxee2, derepRs_maxee2, verbose=TRUE)
  saveRDS(maxeeon_mergers_pseudo,file.path(output_results,"maxeeon_mergers_pseudo.rds"))
}

#Inspect the merger data.frame from the first sample
head(maxeeon_mergers_pseudo[[1]])

--------------------------------------------------------------------------------
###Construct a sequence table
maxeeon_mergers_pseudo
if (file.exists(file.path(output_results,"maxeeon_seqtab_pseudo.rds"))) {
  maxeeon_seqtab_pseudo=readRDS(file.path(output_results,"maxeeon_seqtab_pseudo.rds"))
} else {
  maxeeon_seqtab_pseudo <- makeSequenceTable(maxeeon_mergers_pseudo)
  saveRDS(maxeeon_seqtab_pseudo, file.path(output_results,"maxeeon_seqtab_pseudo.rds"))
}

# Inspect distribution of sequence lengths
table(nchar(getSequences(maxeeon_seqtab_pseudo)))

--------------------------------------------------------------------------------
#Remove chimeras
#PSUEDO POOLING
maxeeon_seqtab_pseudo
if (file.exists(file.path(output_results,"maxeeon_seqtabpseudo_nochim.rds"))) {
  maxeeon_seqtabpseudo_nochim=readRDS(file.path(file.path(output_results,"maxeeon_seqtabpseudo_nochim.rds")))
} else {
  maxeeon_seqtabpseudo_nochim <- removeBimeraDenovo(maxeeon_seqtab_pseudo, method="consensus", multithread=TRUE, verbose=FALSE)
  saveRDS(maxeeon_seqtabpseudo_nochim,file.path(output_results,"maxeeon_seqtabpseudo_nochim.rds"))
}
dim(maxeeon_seqtabpseudo_nochim)

--------------------------------------------------------------------------------
#
##Track reads through the pipeline 
##PSEUDO POOLING

getN <- function(x) sum(getUniques(x))
track_maxeeon_pseudo <- cbind(out_maxee2, sapply(dadaF_maxee2, getN), sapply(dadaR_maxee2, getN), sapply(maxeeon_mergers_pseudo, getN),
                              rowSums(maxeeon_seqtabpseudo_nochim))
colnames(track_maxeeon_pseudo) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track_maxeeon_pseudo) <- sample.names
head(track_maxeeon_pseudo)
track_maxeeon_pseudo = as.data.frame(track_maxeeon_pseudo)
write.csv2(track_maxeeon_pseudo, file.path(output_results,"track_maxeeon_pseudo.csv"))

--------------------------------------------------------------------------------
#Export dada2 results to QIIME2 for taxonomy assignment
#write ASVs to FASTA file
#create a sub-directory to store the dada2 outputs
dir.create(file.path(output_results, "maxeeon_dada2"))

#PSEUDO POOLING
asv_ids_pseudo <- paste0("ASV13_", seq(length(getSequences(maxeeon_seqtabpseudo_nochim))))
uniquesToFasta(maxeeon_seqtabpseudo_nochim, file.path(output_results, "maxeeon_dada2/maxeeon_pseudo_asvs.fasta"), ids = asv_ids_pseudo)

#Export counts table to biom
library(biomformat)
seqtabpseudo.asvids <- t(maxeeon_seqtabpseudo_nochim)
rownames(seqtabpseudo.asvids) <- asv_ids_pseudo
maxeeon_stpseudo.biom <- make_biom(seqtabpseudo.asvids)
write_biom(maxeeon_stpseudo.biom, file.path(output_results, "maxeeon_dada2/maxeeon_table_pseudo.biom"))

##new export using ASV sequence instead of ASV number
dir.create(file.path(output_results, "mock_community"))
uniquesToFasta(maxeeon_seqtabpseudo_nochim, file.path(output_results, "mock_community/mockcommunity_pseudo_asvs.fasta"), 
              ids = colnames(maxeeon_seqtabpseudo_nochim))
