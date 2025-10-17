## Load required library
library(dada2)
set.seed(531)

## Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
path <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), 'fastq')
fns <- list.files(path, pattern = ".fastq.gz$", full.names = TRUE)
fns <- sort(fns)

## Separate forward and reverse reads
fnFs <- fns[grepl("_1.fastq.gz", fns)]
fnRs <- fns[grepl("_2.fastq.gz", fns)]

## Plot quality profiles of the first two samples
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

## Filter and trim reads (matches your methods)
filt_path <- file.path(path, "filtered")
dir.create(filt_path, showWarnings = FALSE)
filtFs <- file.path(filt_path, basename(fnFs))
filtRs <- file.path(filt_path, basename(fnRs))

out <- filterAndTrim(
  fnFs, filtFs,
  fnRs, filtRs,
  truncLen = c(266, 261),
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE,
  pool = "pseudo"
)
write.csv(out, file.path(filt_path, "filtering_summary.csv"))

## Learn error rates
errF <- learnErrors(filtFs, nbases = 1e5, multithread = TRUE)
errR <- learnErrors(filtRs, nbases = 1e5, multithread = TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

## Dereplicate reads
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
names(derepFs) <- sapply(strsplit(basename(filtFs), "_1.fastq.gz"), `[`, 1)
names(derepRs) <- sapply(strsplit(basename(filtRs), "_2.fastq.gz"), `[`, 1)

## Denoise
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = "pseudo")
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = "pseudo")

## Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

## Make sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE)

## Save RDS for downstream use and export for QIIME2
dir.create("final16S", showWarnings = FALSE)
saveRDS(seqtab_nochim, "final16S/seqtab_all_no_chimeras.rds")
write.table(t(seqtab_nochim), "final16S/seqtab_all_no_chimeras.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
uniquesToFasta(seqtab_nochim, fout='final16S/repseqs_perrunB.fasta', ids=colnames(seqtab_nochim))