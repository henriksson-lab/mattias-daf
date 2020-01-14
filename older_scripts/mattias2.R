map_entrez_genesym <- read.csv("map_entrez_genename.csv")
map_entrez_genesym
colnames(map_entrez_genesym) <- c("ensid","GB_ACC","SYMBOL")



###### FACS:
#First (cd3/cd14)-, cd19+
#
#Then:
#  cd38+igd- (GC)
#  cd38-igd- (BMem)
#
# Finally: daf+ vs daf-     daf = CD55

# Review of CD55:
#   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5833118/
#   https://www.uniprot.org/uniprot/P08174  CD55, GPI anchor 
#   Binds to CD97

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
library(topGO)
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


#clariom example ish
#https://support.bioconductor.org/p/93272/

format_de <- function(x){
  x <- x[!is.na(x$SYMBOL),]
  x <- x[str_sub(x$SYMBOL,0,1) == str_to_upper(str_sub(x$SYMBOL,0,1)),]
  subset(x,select=c("SYMBOL","logFC","P.Value"))
}

mergeprobe <- function(x){
  merge(x,
        sqldf("select SYMBOL,logFC,min(`P.Value`) as `P.Value` from x group by SYMBOL"))
}

plot_volcano <- function(x,cex=0.8,pcut=0.001,fname="",xlab="logFC"){
  #x <- format_de(deg_vs2)
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
  #tT
}

###############################################################################################
########### Juice mattias DE genes #############################################################
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
deg_vsDAFinside   <- topTable(efit, coef=2, number = 100000)  ## vs DAF_fraction
#deg_vsBtype <- topTable(efit, coef=3, number = 100000)  ## vs DAF_fraction

table_daf_inside <- format_de(deg_vsDAFinside)
plot_volcano(format_de(deg_vsDAFinside),pcut = 0.01,
             fname="out/forGC_DAFplus_vs_minus.pdf", 
             xlab="in GC, compare DAF, High means more in DAF+")

#pData(normData)

#### Higher inside
#CASP10


#### Higher outside


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
#tio$SYMBOL[tio$p_inside + tio$p_outside < 0.005*2]

###The most in DAF+, from the most to less
tio <- tio[order(tio$fc_inside + tio$fc_outside,decreasing = TRUE),]
as.data.frame(tio$SYMBOL[tio$fc_inside + tio$fc_outside > 0.8*2])
#tio$SYMBOL[tio$p_inside + tio$p_outside < 0.005*2]



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

# pData(normData)
# pData(d)
# design

### Higher outside


#### Higher inside


d <- normData[, pData(normData)$DAF_fraction=="Daf-"]
design <- model.matrix(~Celltype+Patient,pData(d))
fit <- lmFit(d, design)  
efit <- eBayes(fit)        
table_dafminus   <- topTable(efit, coef=2, number = 100000)

table_dafminus <- format_de(table_dafminus)
plot_volcano(format_de(table_dafminus),pcut=0.001,
  fname="out/forDAFminus_inside_vs_outside.pdf", 
  xlab="High means more GC fraction")



#### cd96!!! fcrl4, cd38, 
# notch2 unclear, but a known gene

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
#tio$SYMBOL[tio$p_inside + tio$p_outside < 0.005*2]

###The most in DAF+, from the most to less
tpm <- tio[order(tpm$fc_plus + tpm$fc_minus,decreasing = TRUE),]
as.data.frame(tpm$SYMBOL[tpm$fc_plus + tpm$fc_minus > 0.8*2])
#tio$SYMBOL[tio$p_inside + tio$p_outside < 0.005*2]





##############################################
########## mara pairwise - mattias ###########
##############################################

mara.mattias <- (read.table("mara/activity_table",stringsAsFactors = FALSE))

pData(normData)$Celltype
pData(normData)$DAF_fraction
##ignore patient for now


meanmara <- NULL
for(celltypei in 1:2){
  for(daftypei in 1:2){

    celltype <- c("GC_B","Bmem")[celltypei]
    daftype <- c("Daf+","Daf-")[daftypei]

    groupname <- paste(c("GCB","Bmem")[celltypei],c("Dafplus","Dafminus")[daftypei])
  
    takei <- pData(normData)$Celltype==celltype & pData(normData)$DAF_fraction==daftype
  
    meanmara <- cbind(meanmara,colMeans(mara.mattias[takei,]))
    colnames(meanmara)[ncol(meanmara)]<-groupname
      
  }
}

meanmara <- as.data.frame(meanmara)

plot(meanmara$`GCB Dafplus`,meanmara$`GCB Dafminus`,cex=0)
text(meanmara$`GCB Dafplus`,meanmara$`GCB Dafminus`,labels = rownames(meanmara))

plot_mara_comp <- function(dirx,diry){
  plot(meanmara[,dirx],meanmara[,diry],cex=0)
  text(meanmara[,dirx],meanmara[,diry],labels = rownames(meanmara))
  
  pdf(sprintf("out/maracomp_mattias/%s_%s.pdf",dirx,diry))
  plot(meanmara[,dirx],meanmara[,diry],cex=0,xlab=dirx,ylab=diry)
  text(meanmara[,dirx],meanmara[,diry],labels = rownames(meanmara),cex=0.5)
  dev.off()
}

#correlated, because fairly similar
plot_mara_comp("GCB Dafplus","GCB Dafminus")
plot_mara_comp("Bmem Dafplus","Bmem Dafminus")

#anticorrelated - fairly different programs
plot_mara_comp("Bmem Dafplus","GCB Dafplus")
plot_mara_comp("Bmem Dafminus","GCB Dafminus")





##############################################
########## mara - mattias - heatmap ##########
##############################################

motifs_interesting <- read.csv("interesting_motif.csv",sep="\t", stringsAsFactors = FALSE)

mara.all.normf <- rowMeans(meanmara)
mara.all.norm <- meanmara
# for(i in 1:ncol(mara.all.norm)){
#   #mara.all.norm[,i] <- mara.all.norm[,i]/mara.all.normf
# }
mara.all.norm$tf <- rownames(meanmara)
colnames(mara.all.norm) <- c("GC-B DAF+","GC-B DAF-","Bmem DAF+","Bmem DAF-","tf")
#melt(mara.all)

mara.all.m <- melt(mara.all.norm, id.vars = c("tf"))
#nba.m <- ddply(nba.m, .(variable), transform, rescale = rescale(value))
#mara.all.m

mara.red <- mara.all.m[mara.all.m$tf %in% motifs_interesting$tf,]

#data <-
#ord <- hclust( dist(scale(meanmara), method = "euclidean"), method = "ward.D" )$order

p <- ggplot(mara.red, aes(variable, tf)) +
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient(low = "white",high = "steelblue")
p
pdf("out/maraheatmap_mattias.pdf",height = 5,width=6)
p
dev.off()














##############################################
########## mara pairwise - mattias ###########
##############################################
# 
# mara.mattias <- (read.table("mara/activity_table",stringsAsFactors = FALSE))
# 
# pData(normData)$Celltype
# pData(normData)$DAF_fraction
# ##ignore patient for now
# 
# meanmara <- NULL
# for(celltypei in 1:2){
#   for(daftypei in 1:2){
#     
#     celltype <- c("GC_B","Bmem")[celltypei]
#     daftype <- c("Daf+","Daf-")[daftypei]
#     
#     groupname <- paste(c("GCB","Bmem")[celltypei],c("Dafplus","Dafminus")[daftypei])
#     
#     takei <- pData(normData)$Celltype==celltype & pData(normData)$DAF_fraction==daftype
#     
#     meanmara <- cbind(meanmara,colMeans(mara.mattias[takei,]))
#     colnames(meanmara)[ncol(meanmara)]<-groupname
#     
#   }
# }
# 
# meanmara <- as.data.frame(meanmara)
# 
# plot(meanmara$`GCB Dafplus`,meanmara$`GCB Dafminus`,cex=0)
# text(meanmara$`GCB Dafplus`,meanmara$`GCB Dafminus`,labels = rownames(meanmara))
# 
# plot_mara_comp <- function(dirx,diry){
#   plot(meanmara[,dirx],meanmara[,diry],cex=0)
#   text(meanmara[,dirx],meanmara[,diry],labels = rownames(meanmara))
#   
#   pdf(sprintf("out/maracomp_mattias/%s_%s.pdf",dirx,diry))
#   plot(meanmara[,dirx],meanmara[,diry],cex=0,xlab=dirx,ylab=diry)
#   text(meanmara[,dirx],meanmara[,diry],labels = rownames(meanmara),cex=0.5)
#   dev.off()
# }
# 
# #correlated, because fairly similar
# plot_mara_comp("GCB Dafplus","GCB Dafminus")
# plot_mara_comp("Bmem Dafplus","Bmem Dafminus")
# 
# #anticorrelated - fairly different programs
# plot_mara_comp("Bmem Dafplus","GCB Dafplus")
# plot_mara_comp("Bmem Dafminus","GCB Dafminus")
# 
# 
# 




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

#Up in DZ: E2F7_E2F1 & E2F2_E2F5 -- also in GC-B DAF-
#But POU* inverted
#MYC very LZ, but not affected by DAF
#The IRF* up in LZ. these are particularly low in GC-B DAF-

#XBP1 moderately up in LZ. mainly in GC-B DAF+



