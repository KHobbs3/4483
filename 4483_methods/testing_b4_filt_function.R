# testing b4 filt function
library(dada2)
setwd("/Volumes/data/khobbs3/control/even_mock+spikes")
path <- "/Volumes/data/khobbs3/control/even_mock+spikes"
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

# dereplicate reads
derepF <- derepFastq(fnFs, verbose = TRUE)
derepR <- derepFastq(fnRs, verbose = TRUE)

# name the derep-class objects by the sample names 
names(derepF) <- sample.names
names(derepR) <- sample.names

# make sequence table
# forward reads
seqtabF <- makeSequenceTable(derepF)
table(nchar(getSequences(seqtabF)))

write.table(seqtabF, file="dada2_before_filt_fwd.txt", sep="\t", col.names=NA, quote=F)

# reverse reads
seqtabR <- makeSequenceTable(derepR)
table(nchar(getSequences(seqtabR)))

write.table(seqtabR, file="dada2_before_filt_rvs.txt", sep="\t", col.names=NA, quote=F)


# import fasta spike-in sequence files ----
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/control/spike-in_seqs/")

spike <- read.table("LC140942.txt")
x <- intersect(seqtabF, spike)
y <- intersect(seqtabR, spike)
View(x)
library(dada2)
setwd("/Volumes/data/khobbs3/control/even_mock+spikes")
path <- "/Volumes/data/khobbs3/control/even_mock+spikes"
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

# dereplicate reads
derepF <- derepFastq(fnFs, verbose = TRUE)
derepR <- derepFastq(fnRs, verbose = TRUE)

# name the derep-class objects by the sample names 
names(derepF) <- sample.names
names(derepR) <- sample.names

# make sequence table
# forward reads
seqtabF <- makeSequenceTable(derepF)
table(nchar(getSequences(seqtabF)))

write.table(seqtabF, file="dada2_before_filt_fwd.txt", sep="\t", col.names=NA, quote=F)

# reverse reads
seqtabR <- makeSequenceTable(derepR)
table(nchar(getSequences(seqtabR)))

write.table(seqtabR, file="dada2_before_filt_rvs.txt", sep="\t", col.names=NA, quote=F)


# import fasta spike-in sequence files ----
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/control/spike-in_seqs/merged/")

spike <- read.table("LC140942_merged.txt")
x <- intersect(seqtabF, spike)
y <- intersect(seqtabR, spike)
View(x)
