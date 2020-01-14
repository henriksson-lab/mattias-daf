#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE38697
#paper  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3447782/

#Here we show that expression of CXCR4 and CD83 also distinguishes LZ and DZ populations in humans. 
#Claim LZ B cells more like lymphomas
#Talks about cell cycle phase LZ/DZ

#cell cycle markers
#   cyclins D2/D3/A2/B1/B2, and AURKA
#surface markers
#   CXCR5, GPR183, and SLAMF1
#DNA edinig
#   AICDA, POLH, and LIG4
#apoptosis
#   BCL2LA1, BIM, and NFKBIA
#TFs
#   MYC, EGR1-3, and BATF


#what about these markers vs the inner center of B GC, DAF+?




# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Fri Aug 9 09:29:54 EDT 2019

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE38697", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "LLLLDDDD"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

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
# design <- model.matrix(~ description + 0, gset)
# colnames(design) <- levels(fl)
# fit <- lmFit(gset, design)


design <- model.matrix(~ description + 0, gset)  ######  + patient
colnames(design)[1:2] <- levels(fl)
fit <- lmFit(gset, design)



#############################
### Comparison: Light over dark
#############################
tT <- compareGEOcontrast(makeContrasts(L-D, levels=design),shownum = 1000,xlab = "Light vs dark, high is LZ", fname="out/compare_lz_dz.pdf")
comparisonLD <- tT

#According to paper: CXCR4hi CD83lo  should be DZ. fits

tT[tT$SYMBOL=="CXCR4",]  #more in dark
tT[tT$SYMBOL=="CD83",]   #more in light

#Here we show that expression of CXCR4 and CD83 also distinguishes LZ and DZ populations in humans. 

#https://www.nature.com/articles/cdd2011158  BCL2A2, important!


###### Up in light (right side)
#CD83
#BCL2A1
#EBI3
#SEMA7A
#TRAF4
#MYC
#DUSP2
#EGR3


###### Up in dark (left side)
#CCNB1, 2 -- proliferation
#PRR11
#




