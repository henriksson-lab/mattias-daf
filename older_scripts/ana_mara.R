mara.a_vs_n <- read.table("mara/active_matrices.a_vs_n")


mara.a_vs_n <- read.table("mara/active_matrices.a_vs_n")
mara.c_vs_a <- read.table("mara/active_matrices.c_vs_a")
mara.c_vs_n <- read.table("mara/active_matrices.c_vs_n")

colnames(mara.a_vs_n) <- c("MOTIF","an")
colnames(mara.c_vs_a) <- c("MOTIF","ca")
colnames(mara.c_vs_n) <- c("MOTIF","cn")


mara.all <- merge(
  merge(mara.a_vs_n,mara.c_vs_a,all=TRUE),
  mara.c_vs_n,all=TRUE)

plot(mara.all$an, mara.all$ca)  #similar
plot(mara.all$an, mara.all$cn)
plot(mara.all$ca, mara.all$cn)  #cn is the outlier here

plot(mara.all$an, mara.all$cn,cex=0)
text(mara.all$an, mara.all$cn,mara.all$MOTIF, cex=0.8)
    #ZEB1
    #POU2F1 / POU3F1    are more important in atypical vs naive

plot(mara.all$ca, mara.all$cn,cex=0)
text(mara.all$ca, mara.all$cn,mara.all$MOTIF, cex=0.8)
#ZEB1
#POU2F1 / POU3F1    are more important in classic vs atypical
#posibly also: MYC  maybe maybe NFKB




##############################################
########## e-life ############################
##############################################

#samples in order for mara: a, c, naive

#in atypical, really different:
#ZEB1 
#POU2F1/POU3F1 
#EPAS/BCL3 maybe
#SP4/PML
#NFKB1
#LEF1
#MYBL1
#IKZF2
#STAT1_STAT3_BCL6
#RUNX3_BCL11A
#EPAS1_BCL3

#in classic
#SP3
#maybe: NRF1,TAF1, SIX5_SMARCC2_HCFC1, YBX1_FOS_NFYC_NFYA_NFYB_CEBPZ, 
#   ELK4_ETV5_ELK1_ELK3_ELF4, ZNF711_TFAP2A_TFAP2D, ZBTB33_CHD2 ?
#STAT5A ?

#naive:
#FLI1?










# ##############################################
# ########## Final plot ########################
# ##############################################
# 
# 
#want mattias data

mara.mattias <- read.table("mara/active_matrices.mattias")
colnames(mara.mattias) <- c("MOTIF","m")

mara.mattias <- read.table("mara/activity_table.mattias")
colnames(mara.mattias) <- c("MOTIF","m")

mara.mattias


motifs_interesting <- read.csv("interesting_motif.csv",sep="\t", stringsAsFactors = FALSE)

motifs_interesting

library(reshape2)

mara.all.normf <- rowMeans(mara.all[,-1])
mara.all.norm <- mara.all
for(i in 2:ncol(mara.all.norm)){
  mara.all.norm[,i] <- mara.all.norm[,i]/mara.all.normf
}
#melt(mara.all)

mara.all.m <- melt(mara.all.norm)
#nba.m <- ddply(nba.m, .(variable), transform, rescale = rescale(value))
#mara.all.m

p <- ggplot(mara.all.m[mara.all.m$MOTIF %in% motifs_interesting$tf,], aes(variable, MOTIF)) +
  geom_tile(aes(fill = value),colour = "white") +
  scale_fill_gradient(low = "white",high = "steelblue")
p
