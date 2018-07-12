#Part Two: ALDEx2 for Differential Expression Analysis ----
library("ALDEx2")
setwd('/Users/kt/Documents/Documents/Undergrad/4/4483E/test/CoDa-Output/Part 2 - ALDEx/')

# cvu = controls v untreated ====
# Step 1: create different pairwise condition vectors since ttests can only run for 2 populations
#(must have subsetted data first):
cvu <- c(rep("controls", length(controls)), 
         rep("untreated", length(untreated))
)

# Step 2: create new aldex df containing only controls and untreated columns
# run codapt1-biplots.R first to get ddf 
aldex.cvu <- ddf[,c(controls, untreated)] 
colnames(aldex.cvu)
dim(aldex.cvu)

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
dim(cvu.x.effect)

# Step 6: merge data
cvu.x.all <- data.frame(cvu.x.tt, cvu.x.effect)

# Step 7: create table of merged data
write.table(cvu.x.all, file="cont_untreat aldex_ttest.txt")

# Step 8: see plots 
pdf("controls_untreated.pdf")
aldex.plot(cvu.x.all, type="MA", test="welch")
dev.off()

# Step 9: get features passing a specified significance level 
sig <- which(cvu.x.all$we.eBH < 0.05)

# Step 10: get significant points that only point in the positive direction
psig <- which(cvu.x.all$we.eBH < 0.05 & cvu.x.all$diff.btw > 0)

# Step 11:
#plot diff btwn vs diff within
#plot significant points in a different color
#add the effect=1 and -1 lines
pdf("with text-control_untreated btw-win.pdf")
aldex.plot(cvu.x.all, type="MW", test="welch")
text(cvu.x.all$diff.win[psig], cvu.x.all$diff.btw[psig], labels=row.names(common3), cex= 0.5, pos=3, col="blue") #only positive sig points
text(cvu.x.all$diff.win[psig], cvu.x.all$diff.btw[psig], labels=row.names(common4), cex= 0.5, pos=3, col="red") #meaning only points where values >mean of sample
dev.off()

pdf("with sig points - control_untreated btw-win.pdf")
aldex.plot(cvu.x.all, type="MW", test="welch")
text(cvu.x.all$diff.win, cvu.x.all$diff.btw, labels=row.names(common3), cex= 0.5, pos=3, col="blue") #all points
text(cvu.x.all$diff.win, cvu.x.all$diff.btw, labels=row.names(common4), cex= 0.5, pos=3, col="red") 
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
cvu.rmvd #98

# cvcop = controls v copaxone ====
cvcop <- c(rep("controls", length(controls)),
           rep("copaxone", length(copaxone))
)
aldex.cvcop <- ddf[,c(controls, copaxone)] 
cvcop.x <- aldex.clr(aldex.cvcop, cvcop, mc.samples=128, verbose=TRUE)
cvcop.x.tt <- aldex.ttest(cvcop.x, cvcop, paired.test=FALSE)
cvcop.x.effect <- aldex.effect(cvcop.x, cvcop, include.sample.summary=TRUE, verbose=TRUE)
cvcop.x.all <- data.frame(cvcop.x.tt, cvcop.x.effect)
cvcop.table <- write.table(cvcop.x.all, file="cont_copax aldex_ttest.txt", sep="\t", quote=F, col.names=NA)

pdf("controls_copaxone.pdf")
aldex.plot(cvcop.x.all, type="MA", test="welch")
dev.off()

sig <- which(cvcop.x.all$we.eBH < 0.05)
psig <- which(cvcop.x.all$we.eBH < 0.05 & cvcop.x.all$diff.btw > 0)

pdf("control_copaxone btw-win.pdf")
aldex.plot(cvcop.x.all, type="MW", test="welch")
dev.off()

cvcop.rmvd <- rownames(aldex.cvcop)[rowSums(aldex.cvcop) == 0]
dim(cvcop.rmvd)
length(cvcop.rmvd)  



#cvint = controls v interferon ====
cvint <- c(rep("controls", length(controls)),
           rep("interferon", length(interferon))
)
aldex.cvint <- ddf[,c(controls, interferon)] 
cvint.x <- aldex.clr(aldex.cvint, cvint, mc.samples=128, verbose=TRUE)
cvint.x.tt <- aldex.ttest(cvint.x, cvint, paired.test=FALSE)
cvint.x.effect <- aldex.effect(cvint.x, cvint, include.sample.summary=TRUE, verbose=TRUE)
cvint.x.all <- data.frame(cvint.x.tt, cvint.x.effect)
write.table(cvint.x.all, file="cont_interf aldex_ttest.txt", sep="\t", quote=F, col.names=NA)

pdf("controls_interferon.pdf")
aldex.plot(cvint.x.all, type="MA", test="welch")
dev.off()

sig <- which(cvint.x.all$we.eBH < 0.05) 
psig <- which(cvint.x.all$we.eBH < 0.05 & cvint.x.all$diff.btw > 0)

pdf("control_interferon btw-win.pdf")
aldex.plot(cvint.x.all, type="MW", test="welch")
dev.off()

cvint.rmvd <- rownames(aldex.cvint)[rowSums(aldex.cvint) == 0]
length(cvint.rmvd)
