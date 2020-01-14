#################################################################################
#################################################################################
#################################################################################
################### ELIFE malaria paper #########################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

#samples in order for mara: a, c, naive

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65928
# paper: https://elifesciences.org/articles/07218

### FACS

#A   atypical B-cell  (CD10- CD19+ CD20+ CD21- CD27-)
#C   classical B-cell (CD10- CD19+ CD20+ CD21+ CD27+)
#N   naive B-cell     (CD10- CD19+ CD20+ CD21+ CD27-)

#### classical, atypical and naive






# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Aug 9 04:55:33 EDT 2019

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE65928", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL16686", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
#gsms <- "22222222222222222222000000000000000000001111111111111111111"
gsms <- "AAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCNNNNNNNNNNNNNNNNNNN"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

############## Patients
sml_patient <- str_split_fixed(pData(gset)$title," ",4)[,4]  #c(1:20,1:20,1:20)
sml_patient <- paste("P", sml_patient, sep="")    
sml_patient <- as.factor(sml_patient)





# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
fl <- as.factor(sml)
gset$description <- fl
gset$patient <- sml_patient
design <- model.matrix(~ description + patient + 0, gset)
#design <- model.matrix(~ description + 0, gset)
colnames(design)[1:3] <- levels(fl)
fit <- lmFit(gset, design)





#############################
### Comparison: Classival over atypical
#############################
#A   atypical B-cell  (CD10- CD19+ CD20+ CD21- CD27-)
#C   classical B-cell (CD10- CD19+ CD20+ CD21+ CD27+)
tT <- compareGEOcontrast(makeContrasts(C-A, levels=design), xlab="elife atyp vs classic, high is atypical", fname="out/elife_atyp_vs_classic.pdf")
comparisonElifeCA <- tT

sum(tT$P.Value<0.001)
#Largest difference of 3 comparisons


tT[tT$SYMBOL=="SELL",]

#CD27 and CD21 up as expected (CD21 = CR2)
tT[tT$SYMBOL=="CR2",]
tT[tT$SYMBOL=="CD27",]

tT[tT$SYMBOL=="SELL",]  #Unchanged

######Atypical genes up:
#PLXNC1 - receptor
#NRCAM - cell adhesion. see https://www.jneurosci.org/content/34/34/11274    connects to Sema 3f, and a plexin
#SIGLEC6 - sialyl binding
#FOXP4 - TF
#NKG7 - Natural killer cell granule protein 
#CST7 - cysteine protease inhibitor
#MAP4K3
#NEB - nebulin. complicated
#ENC1 - oxidative stress, toward Nrf2
#SYT1 - endocytotis, Ca2+ sensor
#FGR - downregulates migration and adhesion

####### Classical genes up
#CXCR5
#CCR7
#IL4R
#ALOX5, arachidonic acid
#CR1 and CR2
#FLI1 - TF
#TRAF5 - 
#EML6 - May modify the assembly dynamics of microtubules, such that microtubules are slightly longer, but more dynamic - for migration?
#CD44
#MYC      *** confirmed by mara
#CASK - poorly studied outside neurons. calcium regulation
#CD38 - well studied.  cyclic ADP ribose hydrolase. calcium???? NAD+ and NAADP

#############################
### Comparison: Classical over Naive
#############################
#C   classical B-cell (CD10- CD19+ CD20+ CD21+ CD27+)
#N   naive B-cell     (CD10- CD19+ CD20+ CD21+ CD27-)
tT <- compareGEOcontrast(makeContrasts(C-N, levels=design))
comparisonElifeCN <- tT

sum(tT$P.Value<0.001)

###### Up in classic
#CD27 goes up as expected
#LYZ    The lysozome
#CXCR3
#SOCS3 
#AIM2   DNA-binder, part of inflammosome
#VSIG1
#CD226   binds  CD112 and CD155   https://www.cell.com/cancer-cell/pdf/S1535-6108(14)00467-X.pdf   for antigen presentation

###### Up in naive
#very little!!!
#
#IL4R (so naive > classic > atypical)
#CD72
#YBX3   Binds also to full-length mRNA and to short RNA sequences containing the consensus site 5'-UCCAUCA-3'   Binds to the GM-CSF promoter




#############################
### Comparison: Atypical over Naive
#############################
#A   atypical B-cell  (CD10- CD19+ CD20+ CD21- CD27-)
#N   naive B-cell     (CD10- CD19+ CD20+ CD21+ CD27-)
tT <- compareGEOcontrast(makeContrasts(A-N, levels=design))
comparisonElifeAN <- tT

sum(tT$P.Value<0.001)
#More balanced than classic vs naive. More genes different. Atypical is further away!

###### Up in atypical
#CST7
#NRCAM
#MAP4K3
#NEB


###### Up in naive
#ALOX5
#CR1, CR2
#IL4R
#CCR7, 5


