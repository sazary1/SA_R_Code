##################################################################################################################################
# SA 
# May 15 2014
# LWRNASEQ analysis gene level
###################################################################################################################################

setwd ("/Users/saeedehazary/Documents/LWRNASEQ")

toc=read.delim("/Users/saeedehazary/Documents/LWRNASEQ/mouse_s1pr1_15-May-2014_TMM.txt", sep="\t", header=T)

###################################################################################################################################

### Creating the design matrix for row data

###################################################################################################################################

prefix="C_AB"

Group <- c('B','C','B','C','C','A','A','C','B','A','C','A','B','A','C','B')
Litter <- c('DO297','DO297','DO297','DO297','DO307','DO307','DO307','DO307','DO307','DO307','DO307','DO310','DO310','DO310','DO310','DO310')
design.matrix <- model.matrix(~ 0 + Group + Litter)
#design.matrix <- model.matrix(~ 0 + Group)

print(design.matrix)

toc.int <- apply(toc, 2, as.integer)
rownames(toc.int) <- rownames(toc)
toc <- toc.int

DGE <- DGEList(counts = toc)
DGE <- estimateGLMCommonDisp(DGE, design.matrix)
DGE <- estimateGLMTrendedDisp(DGE, design.matrix)
DGE <- estimateGLMTagwiseDisp(DGE, design.matrix)


fit <- glmFit(DGE, design.matrix)


################
### C vs A B

#################

lrt_C_AB <- glmLRT(fit, contrast = c(-0.5, -0.5, 1, 0, 0))

exp_C_AB=data.frame(topTags(lrt_C_AB, n=nrow(toc)))
topTags(lrt_C_AB, n=nrow(toc))

################
### A vs B 

#################

lrt_A_B <- glmLRT(fit, contrast = c(1, -1, 0, 0, 0))

exp_A_B=data.frame(topTags(lrt_A_B, n=nrow(toc))) 


##################

### B vas C 

##################

lrt_B_C <- glmLRT(fit, contrast = c(0, 1, -1, 0, 0))

exp_B_C=data.frame(topTags(lrt_B_C, n=nrow(toc))) 


##################

### A vs C clean

##################

lrt_A_C <- glmLRT(fit, contrast = c(1, 0, -1,0,0))

exp_A_C=data.frame(topTags(lrt_A_C, n=nrow(toc))) 


###################################################################################################################################

### Heatmap

##################################################################################################################################


########################################################

### ABC list

########################################################

exp_ABC=rbind(data.frame(subset(exp_A_B, exp_A_B$FDR<=0.05)),
              data.frame(subset(exp_A_C, exp_A_C$FDR<=0.05)) ,
              data.frame(subset(exp_B_C, exp_B_C$FDR<=0.05)))  



########################################################


heatmap.file=data.frame(subset(LW.gene, rownames(LW.gene) %in% rownames(exp_ABC)))


#rownames(heatmap.file)=heatmap.file[,1]
#heatmap.file=heatmap.file[, -1] 

install.packages("gplots")
library(gplots)
library(seriation)

f<-colorRampPalette(c("blue", "white", "firebrick"))
col.all <- f(256)

pdf("heatmap_LW_gene.pdf")

finalmat <- scale(t(cpm(heatmap.file, prior.count=2, log=TRUE)))
finalmat[finalmat > 1] = 1
finalmat[finalmat < -1] = -1
count_heatmap<- hmap(as.matrix(finalmat),   cexRow = 0.8, cexCol = 0.005, margins = c(10,5), options = list (prop=TRUE, col = col.all, main="gene LW", cex.main=0.5))

dev.off()


#######################

FC=2^exp_C_AB$logFC
exp_C_AB=data.frame(cbind(exp_C_AB, FC))

FC=2^exp_A_B$logFC
exp_A_B=data.frame(cbind(exp_A_B, FC))

FC=2^exp_A_C$logFC
exp_A_C=data.frame(cbind(exp_A_C, FC))

FC=2^exp_B_C$logFC
exp_B_C=data.frame(cbind(exp_B_C, FC))

for (i in 1:nrow(exp_C_AB)){
  
if (i==1){
Rank_A_B=which(rownames(exp_A_B)==rownames(exp_C_AB)[i])
logFC_A_B=exp_A_B[Rank_A_B, 1]
}
else {
r=which(rownames(exp_A_B)==rownames(exp_C_AB)[i])
l=exp_A_B[r, 1]
Rank_A_B=c(Rank_A_B,r)
logFC_A_B=c(logFC_A_B, l)
FC_A_B=2^ logFC_A_B
}
}



for (i in 1:nrow(exp_C_AB)){
  
if (i==1){
Rank_B_C=which(rownames(exp_B_C)==rownames(exp_C_AB)[i])
logFC_B_C=exp_B_C[Rank_B_C, 1]
}
else {
r=which(rownames(exp_B_C)==rownames(exp_C_AB)[i])
Rank_B_C=c(Rank_B_C,r)
l=exp_B_C[r, 1]
logFC_B_C=c(logFC_B_C,l)
FC_B_C=2^ logFC_B_C
}
}


for (i in 1:nrow(exp_C_AB)){
  
if (i==1){
Rank_A_C=which(rownames(exp_A_C)==rownames(exp_C_AB)[i])
logFC_A_C=exp_A_C[Rank_A_C, 1]
}
else {
r=which(rownames(exp_A_C)==rownames(exp_C_AB)[i])
Rank_A_C=c(Rank_A_C,r)
l=exp_A_C[r, 1]
logFC_A_C=c(logFC_A_C,l)
FC_A_C=2^ logFC_A_C
}
}

exp_C_AB_rank=data.frame(cbind(exp_C_AB, Rank_A_B, Rank_A_C,Rank_B_C, logFC_A_B, FC_A_B, logFC_A_C, FC_A_C, logFC_B_C, FC_B_C))

FC=2^ exp_C_AB$logFC