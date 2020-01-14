
#################################################################################
#################################################################################
#################################################################################
################### PLOS malaria paper ##########################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################

#Paper: https://www.ncbi.nlm.nih.gov/pubmed/25993340
#Probes not expressed above background (normalized intensity of 128) in either sample group were removed from the data set


#typical MBCs   (CD19+CD20+CD21-CD27-IgG+IgD-)
#classical MBCs (CD19+CD20+CD21+CD27+IgG+IgD-)
#transitional B cell (CD19+CD20+CD10+), 

# CD10 is present on pre-B, disappears on mature


# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Aug 9 04:57:52 EDT 2019

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE64493", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL16699", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
#gsms <- "000000222222"
gsms <- "CCCCCCAAAAAA"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# patient names for all samples
gsms <- "abcdefabcdef"
sml_patient <- c()
for (i in 1:nchar(gsms)) { sml_patient[i] <- substr(gsms,i,i) }
sml_patient <- paste("P", sml_patient, sep="")    # set group names
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
#sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
gset$patient     <- sml_patient

design <- model.matrix(~ description + patient + 0, gset)
colnames(design)[1:2] <- levels(fl)
fit <- lmFit(gset, design)




#############################
### Comparison: Classic vs atypical
#############################
tT <- compareGEOcontrast(makeContrasts(C-A, levels=design),shownum = 1000)
comparisonPlosCA <- tT

#### in the other paper
#A   atypical B-cell  (CD10- CD19+ CD20+ CD21- CD27-)
#C   classical B-cell (CD10- CD19+ CD20+ CD21+ CD27+)

#### Here
#typical MBCs   (CD19+CD20+CD21-CD27-IgG+IgD-)
#classical MBCs (CD19+CD20+CD21+CD27+IgG+IgD-)
### what about CD10??


###### Up in classical
#SELL
#ALOX5
#IL11


###### Up in atypical
#IL10 , really up!
#HSPB1
#FOXF1
#BCL6
#GSTM3
#FCRL6
#IKZF2
#ITGB2
#ZAP70
#FGR
#





#############################
### Comparison: Classic vs atypical    merged Elife & Plos
#############################

comparisonPlosCA


comparisonPlosCA <- mergeprobe(comparisonPlosCA)
comparisonElifeCA <- mergeprobe(comparisonElifeCA)
colnames(comparisonPlosCA)  <- c("p_plos","fc_plos","SYMBOL")
colnames(comparisonElifeCA) <- c("p_elife","fc_elife","SYMBOL")

compTotalCA <- merge(comparisonPlosCA, comparisonElifeCA)

#They correlate!
plot(compTotalCA$fc_plos, compTotalCA$fc_elife)

compTotalCA <- compTotalCA[compTotalCA$fc_plos*compTotalCA$fc_elife>0,]
compTotalCA$logFC <- (compTotalCA$fc_plos+compTotalCA$fc_elife)/2
compTotalCA$P.Value <- sqrt(compTotalCA$p_plos*compTotalCA$p_elife)

plot_volcano(compTotalCA,pcut = 1e-8,xlab = "foldFC", fname="out/plos_alife_atyp_vs_classical.pdf")

#compTotalCA[compTotalCA$SYMBOL=="SELL",]


###### Up in classical
#CD27, good!
#CR1, CR2
#ALOX5
#CD44
#CD24
#CASK
#MYC
#MAPKAPK2
#CCR7
#CD96
#S1PR1    https://en.wikipedia.org/wiki/S1PR1
      #S1PR1-deficient thymocytes do not emigrate from the thymus
      #ligand is S1P. GPCR, RAS, ERK


###### Up in atypical
#TOX2
#NRCAM
#BCL11B
#FCRL3
#FGR
#FGFR1
#NKG7
#CST7
#ENC1
#MAP4K3
#ZAP70
#CD68
#SYT1
#IL10
#NEB
#RBPMS2 - binds NOG, connected to BMP
#https://en.wikipedia.org/wiki/ZBTB32    known in immune system. t1d. 
