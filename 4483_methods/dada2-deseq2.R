# DESeq2 Workflow

#source("http://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
library( "DESeq2" )

# 2.4.3 build DESeqDataSet from a counts table ----
## import dada2 counts table:
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/Dada2-Output/")
counts <- read.table("dada2_nochim_tax.txt", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, row.names=1, comment.char="")
counts.d <- counts[,1:ncol(counts)-1]
dim(counts.d)

## subsetting data:
controls <- c("SRR3501908", "SRR3501909", "SRR3501910", "SRR3501911", "SRR3501912", "SRR3501913", "SRR3501914", "SRR3501915", "SRR3501916", "SRR3501925", "SRR3501936", "SRR3501945", "SRR3501946", "SRR3501947", "SRR3501948", "SRR3501949", "SRR3501950", "SRR3501951", "SRR3501952", "SRR3501953", "SRR3501954", "SRR3501955", "SRR3501956", "SRR3501957", "SRR3501958", "SRR3501959", "SRR3501969", "SRR3501973", "SRR3501974", "SRR3501975", "SRR3501976", "SRR3501977", "SRR3501978", "SRR3501979", "SRR3501980", "SRR3501981", "SRR3501982", "SRR3501983", "SRR3501984", "SRR3501985", "SRR3501986", "SRR3501987", "SRR3501991", "SRR3502002")
untreated <- c("SRR3501917", "SRR3501918", "SRR3501921", "SRR3501923", "SRR3501924", "SRR3501932", "SRR3501933", "SRR3501944", "SRR3501960", "SRR3501961", "SRR3501962", "SRR3501963", "SRR3501964", "SRR3501965", "SRR3501966", "SRR3501967", "SRR3501988", "SRR3501989", "SRR3501990", "SRR3501992", "SRR3501993", "SRR3501994", "SRR3501995", "SRR3501996", "SRR3501997", "SRR3501998", "SRR3501999", "SRR3502000", "SRR3502010")
copaxone <- c("SRR3501919", "SRR3501927", "SRR3501930", "SRR3501931", "SRR3501935", "SRR3501940", "SRR3501942", "SRR3501943", "SRR3501970", "SRR3501971", "SRR3502001", "SRR3502003", "SRR3502005", "SRR3502007")
interferon <- c("SRR3501920", "SRR3501922", "SRR3501926", "SRR3501928", "SRR3501929", "SRR3501934", "SRR3501937", "SRR3501938", "SRR3501939", "SRR3501941", "SRR3501968", "SRR3501972", "SRR3502004", "SRR3502006", "SRR3502008", "SRR3502009", "SRR3502011", "SRR3502012")

## import metadata table:
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/")
sampleInfo <- read.csv("sampleInfo.csv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, row.names=1, comment.char="")

library(stringr)
sampleInfo <- str_split_fixed(rownames(sampleInfo), ",", n = 2)
sampleInfo <- as.data.frame(sampleInfo)
colnames(sampleInfo) <- c("run", "condition")

sampleInfo2 <- DataFrame(sampleInfo[,2])
sampleInfo2

## change rownames so that they are the same as colnames of counts (check to make sure they are the same before running function)
rownames(sampleInfo2) <- NULL
rownames(sampleInfo2)  <- colnames(counts.d)
colnames(sampleInfo2) <- "condition"
head(sampleInfo2)

## combine matrix and metatable (counts and sampleInfo2)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = counts.d,
  colData = sampleInfo2,
  design= ~ condition
)
head(ddsFullCountTable)



# 3.1 Preparing data for DESeq2 pipeline ----
## did not subset by col since only one col represents col metadata (colData)
ddsFullCountTable$condition <- droplevels(ddsFullCountTable$condition)


# 3.2 Running DESeq2 pipeline ----
deseq.full <- DESeq(ddsFullCountTable)
?DESeq



# 3.3 Obtaining results ----
deseq.tax <- deseq.full
rownames(deseq.tax)<-paste(counts$tax.vector, rownames(deseq.tax), sep=":")

res <- results(deseq.tax, contrast=c("condition", "untreated", "controls"))
res

setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/deseq2-output/")
write.table(res.tax, file="deseqcvu.txt", sep="\t", quote=F, col.names=NA)

## obtain other "contrasts" (comparisons between conditions):
### copaxone v controls (control must be named second to appear on denominator)
res.cvcop <- results(deseq.tax, contrast=c("condition", "copaxone", "controls"))
res.cvcop

write.table(res.cvcop, file="deseqcvcop.txt", sep="\t", quote=F, col.names=NA)

### interferon v controls
res.cvint <- results(deseq.tax, contrast=c("condition", "interferon", "controls"))
res.cvint

write.table(res.cvint, file="deseqcvint.txt", sep="\t", quote=F, col.names=NA)



# 3.5 Inspecting results ----
# find number of significant results BEFORE adjusted p
sum( res$pvalue < 0.01, na.rm=TRUE )
sum( res.cvcop$pvalue < 0.01, na.rm=TRUE )
sum( res.cvint$pvalue < 0.01, na.rm=TRUE )

# find number of significant results AFTER adjusted p
sum( res$padj < 0.1, na.rm=TRUE )
sum( res.cvcop$padj < 0.1, na.rm=TRUE )
sum( res.cvint$padj < 0.1, na.rm=TRUE )

# subset the results to sort by log2 fold change est to get sig genes with strongest DECREASE/DOWN REG
resSig <- res[ which(res$padj < 0.1 ), ]
head( resSig[ order( resSig$log2FoldChange ), ] )
write.table(resSig, file="cvu-sig-down.txt", sep="\t", quote=F, col.names=NA)

resSig.cvcop <- res.cvcop[ which(res.cvcop$padj < 0.1 ), ]
head( resSig.cvcop[ order( resSig.cvcop$log2FoldChange ), ] )
write.table(resSig.cvcop, file="cvcop-sig-down.txt", sep="\t", quote=F, col.names=NA)

resSig.cvint <- res.cvint[ which(res.cvint$padj < 0.1 ), ]
head( resSig.cvint[ order( resSig.cvint$log2FoldChange ), ] )
write.table(resSig.cvint, file="cvint-sig-down.txt", sep="\t", quote=F, col.names=NA)

# subset to get strongest INCREASE/UP-REG
resSigUp <- tail( resSig[ order( resSig$log2FoldChange ), ] )
write.table(resSigUp, file="cvu-sig-up.txt", sep="\t", quote=F, col.names=NA)

resSigUp.cvcop <- tail( resSig.cvcop[ order( resSig.cvcop$log2FoldChange ), ] )
write.table(resSigUp.cvcop, file="cvcop-sig-up.txt", sep="\t", quote=F, col.names=NA)

resSigUp.cvint <- tail( resSig.cvint[ order( resSig.cvint$log2FoldChange ), ] )
write.table(resSigUp.cvint, file="cvint-sig-up.txt", sep="\t", quote=F, col.names=NA)



# 3.6 plots ----
## MA plots:
pdf("deseq2 cvu MA plot.pdf")
plotMA( res, alpha=0.1, main="MA plot for controls and untreated individuals", ylim = c(-1,1))
dev.off()

pdf("deseq2 cvcop MA plot.pdf")
plotMA( res.cvcop, alpha=0.1, main="MA plot for controls and copaxone-treated individuals", ylim = c(-1,1))
dev.off()

pdf("deseq2 cvint MA plot.pdf")
plotMA( res.cvint, alpha=0.1, main="MA plot for controls and interferon-treated individuals", ylim = c(-1,1))
dev.off()


## Dispersion plots:
dds <- estimateSizeFactors(ddsFullCountTable)
dds <- estimateDispersions(dds)

pdf("deseq2 dispersion plot.pdf")
plotDispEsts( dds, ylim = c(1e-3, 1e5))
dev.off()


## Volcano plots:
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))

pdf("cvcop volcano plot.pdf")
plot(res.cvcop$log2FoldChange, -log10(res.cvcop$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res$log2FoldChange) > 2.5 & res$padj < alpha
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.3)
dev.off()


