#Part Two: ALDEx2 for Differential Expression Analysis ----
library("ALDEx2")
setwd('/Users/kt/Documents/Documents/Undergrad/4/4483E/test/mothuranalysis')
mothur <- read.table("SRP075039_taxonomy_abundances_v3.0.tsv", header=T, sep="\t", stringsAsFactors=F, quote = "", check.names=F, row.names=1, comment.char="")

### Subsetting data:
controls <- c("SRR3501908", "SRR3501909", "SRR3501910", "SRR3501911", "SRR3501912", "SRR3501913", "SRR3501914", "SRR3501915", "SRR3501916", "SRR3501925", "SRR3501936", "SRR3501945", "SRR3501946", "SRR3501947", "SRR3501948", "SRR3501949", "SRR3501950", "SRR3501951", "SRR3501952", "SRR3501953", "SRR3501954", "SRR3501955", "SRR3501956", "SRR3501957", "SRR3501958", "SRR3501959", "SRR3501969", "SRR3501973", "SRR3501974", "SRR3501975", "SRR3501976", "SRR3501977", "SRR3501978", "SRR3501979", "SRR3501980", "SRR3501981", "SRR3501982", "SRR3501983", "SRR3501984", "SRR3501985", "SRR3501986", "SRR3501987", "SRR3501991", "SRR3502002")
untreated <- c("SRR3501917", "SRR3501918", "SRR3501921", "SRR3501923", "SRR3501924", "SRR3501932", "SRR3501933", "SRR3501944", "SRR3501960", "SRR3501961", "SRR3501962", "SRR3501963", "SRR3501964", "SRR3501965", "SRR3501966", "SRR3501967", "SRR3501988", "SRR3501989", "SRR3501990", "SRR3501992", "SRR3501993", "SRR3501994", "SRR3501995", "SRR3501996", "SRR3501997", "SRR3501998", "SRR3501999", "SRR3502000", "SRR3502010")
copaxone <- c("SRR3501919", "SRR3501927", "SRR3501930", "SRR3501931", "SRR3501935", "SRR3501940", "SRR3501942", "SRR3501943", "SRR3501970", "SRR3501971", "SRR3502001", "SRR3502003", "SRR3502005", "SRR3502007")
interferon <- c("SRR3501920", "SRR3501922", "SRR3501926", "SRR3501928", "SRR3501929", "SRR3501934", "SRR3501937", "SRR3501938", "SRR3501939", "SRR3501941", "SRR3501968", "SRR3501972", "SRR3502004", "SRR3502006", "SRR3502008", "SRR3502009", "SRR3502011", "SRR3502012")

mothurd<-mothur[,1:105]
rownames(mothurd)<-paste(mothur$SampleID, rownames(mothur))
dim(mothurd)

# controls versus untreated ====
# Step 1: create different pairwise condition vectors since ttests can only run for 2 populations
cvu <- c(rep("controls", length(controls)), 
         rep("untreated", length(untreated))
)

# Step 2: create new aldex df containing only controls and untreated columns
aldex.cvu <- mothurd[,c(controls, untreated)] 
colnames(aldex.cvu)

#aldex: 
# removed rows with sums == 0
# computed centre with all features - converts instances using centred log-ratio transform
# generates Monte Carlo samples of Dirichlet distribution for each sample

# Step 3: get the clr values
cvu.x <- aldex.clr(aldex.cvu, cvu, mc.samples=128, verbose=TRUE) 

# Step 4: perform t-tests: Welches, Wilcoxon, and Benjamini-Hochberg multiple test correction
cvu.x.tt <- aldex.ttest(cvu.x, cvu, paired.test=FALSE)

# Step 5: estimate effect size
cvu.x.effect <- aldex.effect(cvu.x, cvu, include.sample.summary=TRUE, verbose=TRUE)

# Step 6: merge data
cvu.x.all <- data.frame(cvu.x.tt, cvu.x.effect)

# Step 7: create table of merged data
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/repeat")
write.table(cvu.x.all, file="cont_untreat aldex_ttest.txt", sep = "\t")

# Step 8: see MA plots, both with test ="welch" and "wilcoxon" 
pdf("controls_untreated_wilcoxon.pdf")
aldex.plot(cvu.x.all, type="MA", test="wilcoxon")
dev.off()

# Step 9: get features passing a specified significance level 
sig <- which(cvu.x.all$we.eBH < 0.05)

# Step 10: get significant points that only point in the positive direction
psig <- which(cvu.x.all$we.eBH < 0.05 & cvu.x.all$diff.btw > 0)

# Step 11:
#plot diff btwn vs diff within
#plot significant points in a different color
#add the effect=1 and -1 lines
pdf("with text-control_untreated_MA_wilcoxon btw-win.pdf")
aldex.plot(cvu.x.all, type="MA", test="wilcoxon")
dev.off()

pdf("with sig points - control_untreated btw-win.pdf")
aldex.plot(cvu.x.all, type="MW", test="welch")
dev.off()

# Step 12: make a table of significant taxa

cvu.table <- cvu.x.all[sig, c(4:7, 10,11)] 
caption = "Table of significant taxa" 
digits=3
label="sig.table"
align=c("1",rep("r",6) )

write.table(cvu.table, file="cvu.table_sig_taxa.txt")

# repeat for other pairwise comparisons (annotated as subset versus subset - svs):

# Step 13: note how many SVs with zeros were removed
cvu.rmvd <- rownames(aldex.cvu)[rowSums(aldex.cvu) == 0]
cvu.rmvd # = 78

# cvcop = controls v copaxone ====
cvcop <- c(rep("controls", length(controls)),
           rep("copaxone", length(copaxone))
)
aldex.cvcop <- mothurd[,c(controls, copaxone)] 
cvcop.x <- aldex.clr(aldex.cvcop, cvcop, mc.samples=128, verbose=TRUE)
cvcop.x.tt <- aldex.ttest(cvcop.x, cvcop, paired.test=FALSE)
cvcop.x.effect <- aldex.effect(cvcop.x, cvcop, include.sample.summary=TRUE, verbose=TRUE)
cvcop.x.all <- data.frame(cvcop.x.tt, cvcop.x.effect)
cvcop.table <- write.table(cvcop.x.all, file="cont_copax aldex_ttest.txt", sep="\t", quote=F, col.names=NA)

pdf("controls_copaxone_MW_wilcoxon.pdf")
aldex.plot(cvcop.x.all, type="MW", test="wilcoxon")
dev.off()


sig <- which(cvcop.x.all$we.eBH < 0.05)
psig <- which(cvcop.x.all$we.eBH < 0.05 & cvcop.x.all$diff.btw > 0)

pdf("control_copaxone MA_wilcoxon.pdf")
aldex.plot(cvcop.x.all, type="MA", test="wilcoxon")
dev.off()

cvcop.rmvd <- rownames(aldex.cvcop)[rowSums(aldex.cvcop) == 0]
dim(cvcop.rmvd) # = 0
length(cvcop.rmvd)  




#cvint = controls v interferon ====
cvint <- c(rep("controls", length(controls)),
           rep("interferon", length(interferon))
)
aldex.cvint <- mothurd[,c(controls, interferon)] 
cvint.x <- aldex.clr(aldex.cvint, cvint, mc.samples=128, verbose=TRUE)
cvint.x.tt <- aldex.ttest(cvint.x, cvint, paired.test=FALSE)
cvint.x.effect <- aldex.effect(cvint.x, cvint, include.sample.summary=TRUE, verbose=TRUE)
cvint.x.all <- data.frame(cvint.x.tt, cvint.x.effect)
write.table(cvint.x.all, file="cont_interf aldex_ttest.txt", sep="\t", quote=F, col.names=NA)

pdf("controls_interferon MA_wilcoxon.pdf")
aldex.plot(cvint.x.all, type="MA", test="wilcoxon")
dev.off()

sig <- which(cvint.x.all$we.eBH < 0.05) 
psig <- which(cvint.x.all$we.eBH < 0.05 & cvint.x.all$diff.btw > 0)

pdf("control_interferon btw-win.pdf")
aldex.plot(cvint.x.all, type="MW", test="welch")
dev.off()

cvint.rmvd <- rownames(aldex.cvint)[rowSums(aldex.cvint) == 0]
length(cvint.rmvd) # = 106 rmvd