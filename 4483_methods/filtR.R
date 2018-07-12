devtools::install_github('bjoris33/filtR')
library(ALDEx2)
library(propr)

# TEST ----
## dada2
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/Dada2-Output/")
mytable <- filtR::filtR(count_file = "dada2_nochim_tax.txt", rho_CO = 0.5, clr_CO = 3)

setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/filtR")
write.table(mytable, file="dada2-rho0.5-clr5.txt", sep="\t")

setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/mothuranalysis/")
mthrtable <- filtR::filtR(count_file = 'SRP076838_taxonomy_abundances_v3.0.tsv', rho_CO = 0.5, clr_CO = 3)

setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/filtR")
write.table(mthrtable, file="mothur-rho0.5-clr3.txt", sep="\t")


## mothur

setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/mothuranalysis/")
mothurf <- mothur[,1:105]
write.table(mothurf, "mothur_illumina_only.txt", sep="\t")
mytable <- filtR::filtR(count_file = "SRP075039_taxonomy_abundances_v3.0.tsv", rho_CO = 0.7, clr_CO = 5)
View(mytable)

setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/test/filtR")
write.table(mytable, file="mothur-rho0.7-clr5.txt", sep="\t")


# CONTROL ----
## dada2
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/control/dada2-output")
dada2.control <- filtR::filtR(count_file = "dada2_nochim_tax.txt", rho_CO = 0.3, clr_CO = 1)



## QIIME
setwd("/Users/kt/Documents/Documents/Undergrad/4/4483E/control")
qiime <- filtR::filtR("control_24_samples_tax.tsv", rho_CO = 0.3, clr_CO = 1)
