# negative control: finding the spike-ins before filtering
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


## make list of files
spikeins <- list.files(pattern= "LC1409")

## make empty directory 
data <- list()

## read all files in directory and add names of the files
for (i in 1:length(spikeins))
{
  data[[i]] <- read.table(spikeins[i], sep ="\t")
}

names(data)<-spikeins

# paste sequences together so they exist in one col
data.1 <- list()
data.2 <- list()


for (i in data)
{
  for (k in 2:nrow(i))
  {
    data.1 <- paste(data.1, i[k,], sep ="")
  }
  data.2 <- c(data.2, data.1)
  data.1 <- list()
}

names(data.2)<-spikeins


# merging datasets to identify spikes in dada2 output----
# for fwd reads (seqtabF):
data.3 <- list()

for (j in data)
{
  for (h in j)
  {
    data.3 <- intersect(seqtabF, data.2)
  }
}

View(data.3)

# for rvs reads (seqtabR):
data.4 <- list()

for (a in data)
{
  for (t in a)
  {
    data.4 <- intersect(seqtabR, data.2)
  }
}

View(data.4)
