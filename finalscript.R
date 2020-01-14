if(FALSE){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("pd.clariom.d.human")
  BiocManager::install("affy")
  BiocManager::install("affycoretools")
  BiocManager::install("GEOquery")
  
  BiocManager::install("org.Bt.eg.db")
  BiocManager::install("AnnotationFuncs")
  
}

library(pd.clariom.d.human)
library(Biobase)
library(sqldf)
library(oligo)
library(arrayQualityMetrics)
library(limma)
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(matrixStats)
library(genefilter)
library(affycoretools)
library(reshape2)
library(Biobase)
library(GEOquery)
library(limma)


###########################################################
## mapping gene names to symbols
map_entrez_genesym <- read.csv("map_entrez_genename.csv")
map_entrez_genesym
colnames(map_entrez_genesym) <- c("ensid","GB_ACC","SYMBOL")


###########################################################
## In: microarray limma DE
## out: limited column set
format_de <- function(x){
  x <- x[!is.na(x$SYMBOL),]
  x <- x[str_sub(x$SYMBOL,0,1) == str_to_upper(str_sub(x$SYMBOL,0,1)),]
  subset(x,select=c("SYMBOL","logFC","P.Value"))
}

###########################################################
## Merge the probes for each gene
mergeprobe <- function(x){
  merge(x,
        sqldf("select SYMBOL,logFC,min(`P.Value`) as `P.Value` from x group by SYMBOL"))
}

###########################################################
## Produce a DE volcano plot
plot_volcano <- function(x,cex=0.8,pcut=0.001,fname="",xlab="logFC"){
  x <- x[x$P.Value<pcut,]
  plot(x$logFC, log10(x$P.Value),cex=0)
  text(x$logFC, log10(x$P.Value),labels=x$SYMBOL,cex=cex, xlab=xlab)
  
  if(fname!=""){
    pdf(fname)
    plot(x$logFC, log10(x$P.Value),cex=0)
    text(x$logFC, log10(x$P.Value),labels=x$SYMBOL,cex=cex*0.5, xlab=xlab)
    dev.off()
  }
}

###########################################################
## Compare GEO data, given contrast
compareGEOcontrast <- function(cont.matrix,shownum=2000,xlab="logFC", fname=""){  ########## note the fit!
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000)
  
  if("Gene.symbol" %in% colnames(tT)){
    #used in light dark
    colnames(tT)[colnames(tT)=="Gene.symbol"] <- "SYMBOL"
    tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","logFC","SYMBOL"))
  } else {
    tT <- tT[tT$GB_ACC!="",]
    tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","logFC","GB_ACC"))
    tT <- subset(merge(tT, map_entrez_genesym),select=c("P.Value","logFC","SYMBOL"))
  }
  
  
  tT <- subset(tT, select=c("P.Value","logFC","SYMBOL"))
  tTs <- mergeprobe(tT)
  tTs <- tTs[order(tTs$P.Value),]
  tTs <- tTs[1:shownum,]
  
  plot(tTs$logFC,log10(tTs$P.Value),cex=0,xlab=xlab)
  text(tTs$logFC,log10(tTs$P.Value),labels=tTs$SYMBOL)
  
  if(fname!=""){
    pdf(fname)
    plot(tTs$logFC,log10(tTs$P.Value),cex=0,xlab=xlab)
    text(tTs$logFC,log10(tTs$P.Value),labels=tTs$SYMBOL,cex=0.5)
    dev.off()
  }
  
  unique(tT[order(tT$P.Value),])
}

###############################################################################################
########### Mattias DE genes (new for this study) #############################################
###############################################################################################

# 
inp_pdata <- read.csv("BEA19P113_JM.pdata.csv", row.names = 1)
inp_pdata$comb_type <- str_c(inp_pdata$DAF_fraction,":",inp_pdata$Celltype)
phenoData <- AnnotatedDataFrame(inp_pdata)
celfiles <- list.files(pattern = "*.CEL","BEA19P113_JM", full.names = TRUE)
raw_data <- oligo::read.celfiles(celfiles,verbose =FALSE, phenoData = phenoData)
normData <- rma(raw_data)
normData <- annotateEset(normData, pd.clariom.d.human)


###############################################################################################
########### DAF diff, inside and outside GC ###################################################
###############################################################################################

#DAF+ is 1. Such logFC>0 means higher in DAF+ fraction

d <- normData[, pData(normData)$Celltype=="GC_B"]
design <- model.matrix(~DAF_fraction+Patient,pData(d))
fit <- lmFit(d, design)  
efit <- eBayes(fit)        
deg_vsDAFinside   <- topTable(efit, coef=2, number = 100000,confint = TRUE)  ## vs DAF_fraction
#deg_vsBtype <- topTable(efit, coef=3, number = 100000)  ## vs DAF_fraction
deg_vsDAFinside$SYMBOL[is.na(deg_vsDAFinside$SYMBOL)] <- ""

table_daf_inside <- format_de(deg_vsDAFinside)
plot_volcano(format_de(deg_vsDAFinside),pcut = 0.01,
             fname="out/forGC_DAFplus_vs_minus.pdf", 
             xlab="in GC, compare DAF, High means more in DAF+")

###############
###############  Plot it for a lot of IGH genes
###############
sum(str_starts(deg_vsDAFinside$SYMBOL,"IGH"),na.rm = TRUE)
sum(str_starts(deg_vsDAFinside$SYMBOL,"IGH"),na.rm = TRUE)

subtab <- deg_vsDAFinside
subtab$SYMBOL <- as.character(subtab$SYMBOL)
subtab$SYMBOL[is.na(subtab$SYMBOL)] <- ""
subtab <- subtab[str_starts(subtab$SYMBOL,"IGH") | str_starts(subtab$SYMBOL,"CD55"),]

subtab <- subtab[order(subtab$logFC),]
subtab$SYMBOL <- factor(subtab$SYMBOL,levels = subtab$SYMBOL)
plot(sort(subtab$logFC))



p<- ggplot(subtab, aes(x=SYMBOL, y=logFC)) + theme(axis.text.x = element_text(angle = 90)) + 
  geom_point()+
  geom_errorbar(aes(ymin=CI.L, ymax=CI.R), width=.2,
                position=position_dodge(.9)) 
print(p)
pdf("new_igh.pdf",w=15)
print(p)
dev.off()



d <- normData[, pData(normData)$Celltype=="Bmem"]
design <- model.matrix(~DAF_fraction+Patient,pData(d))
fit <- lmFit(d, design)  
efit <- eBayes(fit)        
deg_vsDAFoutside   <- topTable(efit, coef=2, number = 100000)  ## vs DAF_fraction
#deg_vsBtype <- topTable(efit, coef=3, number = 100000)  ## vs DAF_fraction

table_daf_outside <- format_de(deg_vsDAFoutside)
plot_volcano(format_de(deg_vsDAFoutside),pcut=0.01)


ti <- table_daf_inside
to <- table_daf_outside
ti <- mergeprobe(ti)
to <- mergeprobe(to)
colnames(ti) <- c("SYMBOL","fc_inside","p_inside")
colnames(to) <- c("SYMBOL","fc_outside","p_outside")

tio <- merge(ti,to)
tio <- tio[tio$p_inside<0.1 | tio$p_outside<0.1,]

plot(tio$fc_inside, tio$fc_outside, cex=0,
     xlab="DAF+ vs -, inside GC (high is plus)", 
     ylab="DAF+ vs -, outside GC (high is +)")
text(tio$fc_inside, tio$fc_outside, labels = tio$SYMBOL, cex=0.8)
pdf("out/compare_daf_2way.pdf")
plot(tio$fc_inside, tio$fc_outside, cex=0,
     xlab="DAF+ vs -, inside GC (high is plus)", 
     ylab="DAF+ vs -, outside GC (high is +)")
text(tio$fc_inside, tio$fc_outside, labels = tio$SYMBOL, cex=0.8)
dev.off()

###The least in DAF+, from the least to more
tio <- tio[order(tio$fc_inside + tio$fc_outside),]
as.data.frame(tio$SYMBOL[tio$fc_inside + tio$fc_outside < -0.8*2])

###The most in DAF+, from the most to less
tio <- tio[order(tio$fc_inside + tio$fc_outside,decreasing = TRUE),]
as.data.frame(tio$SYMBOL[tio$fc_inside + tio$fc_outside > 0.8*2])





###############################################################################################
########### DAF diff, inside GC, focus on IGH per patient  ####################################
###############################################################################################


get_for_p <- function(pid, celltype="GC_B"){
  d1 <- normData[, pData(normData)$Celltype==celltype & pData(normData)$Patient==pid & pData(normData)$DAF_fraction=="Daf+"]
  t1 <- unname(exprs(d1))[,1]
  d2 <- normData[, pData(normData)$Celltype==celltype & pData(normData)$Patient==pid & pData(normData)$DAF_fraction=="Daf-"]
  t2 <- unname(exprs(d2))[,1]
  t1/t2  
}


temp <- featureData(normData)@data
igh_perp <- data.frame(
  p1=get_for_p("p1"), 
  p3=get_for_p("p3"),
  meanprobe=apply(exprs(normData),1,mean),
  SYMBOL=temp$SYMBOL)
#hist(igh_perp$meanprobe)   #above 7?
igh_perp <- igh_perp[!is.na(igh_perp$SYMBOL) & str_starts(igh_perp$SYMBOL,"IGH") & igh_perp$meanprobe>0,]
igh_perp$p1 <- igh_perp$p1 - mean(igh_perp$p1)  #normalization of some sort
igh_perp$p3 <- igh_perp$p3 - mean(igh_perp$p3)

pdf("out/igh_gc.pdf")
plot(
  igh_perp$p1,   #not normalized enough
  igh_perp$p3,
  pch=19,
  xlab="Patient 1 GC_B Daf+/Daf- IGH*",
  ylab="Patient 3 GC_B Daf+/Daf- IGH*"
)
dev.off()
# text(
#   igh_perp$p1,   #not normalized enough
#   igh_perp$p3,
#   labels = igh_perp$SYMBOL,
#   cex=0.5,col="red"   #IGHV3-73  stands out
# )






temp <- featureData(normData)@data
igh_perp <- data.frame(
  p1=get_for_p("p1","Bmem"), 
  p3=get_for_p("p3","Bmem"),
  meanprobe=apply(exprs(normData),1,mean),
  SYMBOL=temp$SYMBOL)
igh_perp <- igh_perp[!is.na(igh_perp$SYMBOL) & str_starts(igh_perp$SYMBOL,"IGH") & igh_perp$meanprobe>0,]
igh_perp$p1 <- igh_perp$p1 - mean(igh_perp$p1)  #normalization of some sort
igh_perp$p3 <- igh_perp$p3 - mean(igh_perp$p3)

pdf("out/igh_bmem.pdf")
plot(
  igh_perp$p1,   #not normalized enough
  igh_perp$p3,
  pch=19,
  xlab="Patient 1 Bmem Daf+/Daf- IGH*",
  ylab="Patient 3 Bmem Daf+/Daf- IGH*"
)
dev.off()


###############################################################################################
########### inside and outside GC, vs DAF ###################################################
###############################################################################################

#GC_B is 1. Such logFC>0 means higher in the GC fraction

d <- normData[, pData(normData)$DAF_fraction=="Daf+"]
design <- model.matrix(~Celltype+Patient,pData(d))
fit <- lmFit(d, design)  
efit <- eBayes(fit)        
table_dafplus   <- topTable(efit, coef=2, number = 100000)

table_dafplus <- format_de(table_dafplus)
plot_volcano(format_de(table_dafplus),pcut=0.01,
             fname="out/forDAFplus_inside_vs_outside.pdf", 
             xlab="High means more GC fraction")


d <- normData[, pData(normData)$DAF_fraction=="Daf-"]
design <- model.matrix(~Celltype+Patient,pData(d))
fit <- lmFit(d, design)  
efit <- eBayes(fit)        
table_dafminus   <- topTable(efit, coef=2, number = 100000)

table_dafminus <- format_de(table_dafminus)
plot_volcano(format_de(table_dafminus),pcut=0.001,
             fname="out/forDAFminus_inside_vs_outside.pdf", 
             xlab="High means more GC fraction")

tp <- table_dafplus
tm <- table_dafminus
tp <- mergeprobe(tp)
tm <- mergeprobe(tm)
colnames(tp) <- c("SYMBOL","fc_plus","p_plus")
colnames(tm) <- c("SYMBOL","fc_minus","p_minus")

tpm <- merge(tp,tm)
tpm <- tpm[tpm$p_plus<0.1 | tpm$p_minus<0.1,]

plot(tpm$fc_plus, tpm$fc_minus, cex=0)
text(tpm$fc_plus, tpm$fc_minus, labels = tpm$SYMBOL, cex=0.8)

###The least in DAF+, from the least to more
tpm <- tpm[order(tpm$fc_plus + tpm$fc_minus),]
as.data.frame(tpm$SYMBOL[tpm$fc_plus + tpm$fc_minus < -0.8*2])

###The most in DAF+, from the most to less
tpm <- tio[order(tpm$fc_plus + tpm$fc_minus,decreasing = TRUE),]
as.data.frame(tpm$SYMBOL[tpm$fc_plus + tpm$fc_minus > 0.8*2])







##############################################
########## mara for mattias, #################
##############################################


mara.mattias <- read.table("mara/active_matrices.mattias")
colnames(mara.mattias) <- c("MOTIF","m")

mara.mattias <- read.table("mara/activity_table.mattias")
colnames(mara.mattias) <- c("MOTIF","m")

motifs_interesting <- read.csv("interesting_motif.csv",sep="\t", stringsAsFactors = FALSE)

mara.all.normf <- rowMeans(mara.all[,-1])
mara.all.norm <- mara.all
for(i in 2:ncol(mara.all.norm)){
  mara.all.norm[,i] <- mara.all.norm[,i]/mara.all.normf
}

mara.all.m <- melt(mara.all.norm)

p <- ggplot(mara.all.m[mara.all.m$MOTIF %in% motifs_interesting$tf,], aes(variable, MOTIF)) +
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient(low = "white",high = "steelblue")
p



##############################################
########## mara - mattias - dzlz #############
##############################################


mara.dzlz <- t(read.table("mara/activity_table.dzlz",stringsAsFactors = FALSE))
mara.dzlz <- data.frame(lz=rowMeans(mara.dzlz[,1:4]), dz=rowMeans(mara.dzlz[,5:8]))
mara.dzlz$tf <- rownames(mara.dzlz)
colnames(mara.dzlz) <- c("LZ","DZ","tf")

mara.dzlz.red <- melt(mara.dzlz, id.vars = c("tf"))
mara.dzlz.red <- mara.dzlz.red[mara.dzlz.red$tf %in% motifs_interesting$tf,]

p <- ggplot(mara.dzlz.red, aes(variable, tf)) +
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient(low = "white",high = "steelblue")
p
pdf("out/maraheatmap_dzdlz.pdf",height = 5,width=5)
p
dev.off()

















##############################################
##############################################
####### Light vs dark. DE analysis ###########
##############################################
##############################################

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
fl <- as.factor(sml)
gset$description <- fl

design <- model.matrix(~ description + 0, gset)  ######  + patient
colnames(design)[1:2] <- levels(fl)
fit <- lmFit(gset, design)



#############################
### Comparison: Light over dark
#############################
tT <- compareGEOcontrast(makeContrasts(L-D, levels=design),shownum = 1000,xlab = "Light vs dark, high is LZ", fname="out/compare_lz_dz.pdf")
comparisonLD <- tT

tT[tT$SYMBOL=="CXCR4",]  #more in dark
tT[tT$SYMBOL=="CD83",]   #more in light



##############################################
##############################################
########## Meta analysis #####################
##############################################
##############################################



normorder <- function(x) x[,c("P.Value","logFC","SYMBOL")]


plot_2de <- function(listao,listbo,pcut=0.001,dotext=TRUE,
                     fname="",
                     xlab=deparse(substitute(listao)),
                     ylab=deparse(substitute(listbo))){
  lista <- listao
  listb <- listbo
  lista <- mergeprobe(lista)
  listb <- mergeprobe(listb)
  colnames(lista) <- c("p_a","fc_a","SYMBOL")
  colnames(listb) <- c("p_b","fc_b","SYMBOL")
  listab <- merge(lista,listb,all=TRUE)
  listab$p_a[is.na(listab$p_a)] <- 1
  listab$p_b[is.na(listab$p_b)] <- 1
  listab$fc_a[is.na(listab$fc_a)] <- 0
  listab$fc_b[is.na(listab$fc_b)] <- 0
  listab <- listab[listab$p_a < pcut | listab$p_b < pcut,]
  
  if(dotext){
    cex<-0
  }
  plot(listab$fc_a, listab$fc_b,cex=cex,xlab=xlab,ylab=ylab)
  if(dotext){
    text(listab$fc_a, listab$fc_b,labels = listab$SYMBOL,cex=0.7)
  }
  
  if(fname!=""){
    cex<-0
    pdf(fname)
    plot(listab$fc_a, listab$fc_b,cex=cex,xlab=xlab,ylab=ylab)
    text(listab$fc_a, listab$fc_b,labels = listab$SYMBOL,cex=0.3)
    dev.off()
  }
  
  
  print(cor.test(listab$fc_a, listab$fc_b))
}



#within Bmem vs Bgc, DAF+ or DAF-. logFC>0 means higher in DAF+ fraction
table_daf_inside  <- normorder(table_daf_inside)
table_daf_outside <- normorder(table_daf_outside)

#within DAF+ or DAF-, compare Bmem vs Bgc. logFC>0 means higher in the GC fraction
table_dafplus     <- normorder(table_dafplus)
table_dafminus    <- normorder(table_dafminus)

plot_2de(table_daf_inside, table_daf_outside, pcut=0.01)  #not correlated. and weaker
plot_2de(table_dafplus, table_dafminus, pcut=0.01)        #correlated


plot_2de(comparisonLD, table_daf_inside,  #correlated. ~0.2     ### we really care about this one
         fname="out/comp2_LD_inside_dafPLUS_vs_MINUS.pdf",
         xlab="Light vs dark, high is LZ",
         ylab="in GC, compare DAF, High means more in DAF+")
plot_2de(comparisonLD, table_daf_outside, #not correlated
         fname="out/comp2_LD_outside_dafPLUS_vs_MINUS.pdf",
         xlab="Light vs dark, high is LZ",
         ylab="outside GC, compare DAF, High means more in DAF+")
plot_2de(comparisonLD, table_dafplus,     #negatively correlated. ~0.3
         fname="out/comp2_LD_dafplus_inside_vs_outside.pdf",
         xlab="Light vs dark, high is LZ",
         ylab="for daf+, inside vs outside, High means more in GC fraction")
plot_2de(comparisonLD, table_dafminus,    #negatively correlated. ~0.3
         fname="out/comp2_LD_dafminus_inside_vs_outside.pdf",
         xlab="Light vs dark, high is LZ",
         ylab="for daf-, inside vs outside, High means more in GC fraction")





########################################
########################################
########################################
########################################


plot_2de_new <- function(listao,listbo,pcut=0.0001,dotext=TRUE,
                         fname="",
                         xlab=deparse(substitute(listao)),
                         ylab=deparse(substitute(listbo))){
  lista <- listao
  listb <- listbo
  lista <- mergeprobe(lista)
  listb <- mergeprobe(listb)
  colnames(lista) <- c("p_a","fc_a","SYMBOL")
  colnames(listb) <- c("p_b","fc_b","SYMBOL")
  listab <- merge(lista,listb,all=TRUE)
  listab$p_a[is.na(listab$p_a)] <- 1
  listab$p_b[is.na(listab$p_b)] <- 1
  listab$fc_a[is.na(listab$fc_a)] <- 0
  listab$fc_b[is.na(listab$fc_b)] <- 0
  listab <- listab[listab$p_a < pcut | listab$p_b < pcut,]
  
  listab <- listab[listab$fc_a !=0 & listab$fc_b != 0,]
  
  
  cex<-0
  pdf(fname)
  plot(listab$fc_a, listab$fc_b,xlab=xlab,ylab=ylab,cex=0.5,pch=19,col="red")
  
  listab <- listab[listab$SYMBOL %in% interesting_gene,]
  text(listab$fc_a, listab$fc_b,labels = listab$SYMBOL,cex=1)
  
  dev.off()
}


interesting_gene <- read.csv("interesting_genes.csv",sep="\t", stringsAsFactors = FALSE)$gene
interesting_gene <- c(interesting_gene,as.character(table_daf_inside$SYMBOL[grep("IGH.*",table_daf_inside$SYMBOL)]))


plot_2de_new(
  comparisonLD, table_daf_inside,  #correlated. ~0.2     ### we really care about this one
  fname="out/comp2new_LD_inside_dafPLUS_vs_MINUS.pdf",
  xlab="Light vs dark, high is LZ",
  ylab="in GC, compare DAF, High means more in DAF+")










########################################
####### Final plots wanted #############
########################################






