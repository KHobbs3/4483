# Control dada2 workflow
# set up ----
source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
a
no
library(dada2)

path <- "/Volumes/data/khobbs3/control/even_mock+spikes"
taxpath<-"/Volumes/data/annotationDB/dada2/silva_nr_v123_train_set.fa.gz" #agrajag
reads<-"demultiplex_reads"

# list files
list.files(path)

# sort fwd and rvs reads
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names=TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names=TRUE))

# get file names only (remove path)
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# check for duplications
any(duplicated(sample.names))

# check read quality of a random subset of 4 samples ----
ids<-round(runif(4,1,length(sample.names)))

setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/control/dada2-output")

pdf("qualprofiles.pdf")
plotQualityProfile(fnFs[ids])
plotQualityProfile(fnRs[ids])
dev.off()

# make file names for filtered reads ----
filtFs <- paste0(reads, "/", sample.names, "-F-filt.fastq")
filtRs <- paste0(reads, "/", sample.names, "-R-filt.fastq")

# filter based on QC ----
out<-filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                   truncLen=c(200,150),		
                   truncQ=2,
                   maxN=0,
                   maxEE=c(2,2),
                   compress=TRUE, verbose=TRUE, multithread=TRUE)

write.table(out, file="after_filter.txt", sep="\t", col.names=NA, quote=F)

# learn error rates ----
errF <- learnErrors(filtFs, multithread=TRUE, randomize=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE, randomize=TRUE)

pdf("err.pdf")
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
dev.off()

# dereplication (combine identicals to make unique sequences) ----
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# infer SVs ----
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# make sequence table
seqtab <- makeSequenceTable(mergers)
table(nchar(getSequences(seqtab)))

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", verbose=TRUE, multithread=TRUE)
dim(seqtab.nochim)

write.table(seqtab.nochim, file="temp_dada2_nochim.txt", sep="\t", col.names=NA, quote=F)

# How many reads made it through pipeline ----
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
write.table(track, file="control-readsout.txt", sep="\t", col.names=NA, quote=F)

# Assign taxonomy ---- 
taxa <- assignTaxonomy(seqtab.nochim, taxpath, multithread=TRUE)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

#merge columns 1 to 6 to get the full taxonomy (to genus)
un.tax <- unname(taxa)
tax.vector <- apply(un.tax, 1, function(x){paste(x[1:6], collapse=":")})

# add new row of taxonomy names then transpose so that its a col
seqtab.nochim.tax<-rbind(seqtab.nochim, tax.vector)
t.seqtab.nochim.tax<-t(seqtab.nochim.tax)

# replace SV rownames with arbitrary numbers
sv.seqs<-rownames(t.seqtab.nochim.tax)
sv.num<-paste("SV", seq(from = 0, to = nrow(t.seqtab.nochim.tax)-1), sep="_")

rownames(t.seqtab.nochim.tax)<-sv.num

write.table(t.seqtab.nochim.tax, file="dada2_nochim_tax.txt", sep="\t", col.names=NA, quote=F)
write.table(sv.seqs, file="sv_seqs.txt", sep="\t", row.names=sv.num, col.names=F,  quote=F)
