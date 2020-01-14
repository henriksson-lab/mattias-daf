########################################
########################################
########################################
########################################

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

plot_2de(comparisonLD, comparisonElifeAN)
plot_2de(comparisonLD, comparisonElifeCA)
plot_2de(comparisonLD, comparisonElifeCN)

plot_2de(comparisonLD, compTotalCA[,c("P.Value","logFC","SYMBOL")])


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

### it would make sense to merge dafplus & dafminus



comparisonTotalCA <- normorder(compTotalCA)

plot_2de(comparisonTotalCA, table_daf_inside)  #not correlated
plot_2de(comparisonTotalCA, table_daf_outside) #very slightly positively correlated?
plot_2de(comparisonTotalCA, table_dafplus)     #not correlated
plot_2de(comparisonTotalCA, table_dafminus)    #not correlated

plot_2de(comparisonTotalCA, comparisonLD)      #not correlated


plot_2de(comparisonElifeAN, table_daf_outside)      #not correlated

#### So:
## * Gene is up in the light zone <=> Gene is up in DAF+ as compared to DAF-, for B_GC (20% correlation). No such correlation B_mem 
##    => DAF+ cells are more LZ-like, but only if within the GC
## * Gene is up in light zone <=> Gene is down in B GC fraction  (30% correlation, independent of DAF status)
##    => This is weird. It means light zone more similar to Bmem?
## * DAF+/- are markers for the same subset of genes, independent of if the cell is in the GC or outside
##    => these genes only have an effect within the GC, or they are set by the activation process but do not influence it
##
##  Atypical B cells do not correlate with any other axis defined. We could at best speculate on individual DE genes but
##  I would not try to present any numbers
##
##  The two datasets on atypical B cells correlate well and appear of high quality
##
##  The CD55-gated microarrays appear of sufficient quality. The limitations are not in sample noise but in the
##  difficulity of interpreting data, especially when different arrays have been used. Even if this was not the case,
##  sample prep introduces batch effects that would be hard or impossible to regress out, unless all measurements are
##  done anew in one single go












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

