##################################################################################################################################
# SA 
# May 15 2014
# LWRNASEQ analysis gene level
###################################################################################################################################

setwd ("/Users/saeedehazary/Documents/LWRNASEQ")

toc=read.delim("/Users/saeedehazary/Documents/LWRNASEQ/mouse_s1pr1_15-May-2014_TMM.txt", sep="\t", header=T)

library(limma)
library(edgeR)
library(gplots)
library(seriation)


###################################################################################################################################

### Creating the design matrix

###################################################################################################################################

Group <- c('B','C','B','C','C','A','A','C','B','A','C','A','B','A','C','B')
Litter <- c('DO297','DO297','DO297','DO297','DO307','DO307','DO307','DO307','DO307','DO307','DO307','DO310','DO310','DO310','DO310','DO310')
design.matrix <- model.matrix(~ 0 + Group + Litter)
#design.matrix <- model.matrix(~ 0 + Group)

print(design.matrix)

######################
# edgeR dispersion

######################

toc.int <- apply(toc, 2, as.integer)
rownames(toc.int) <- rownames(toc)
toc <- toc.int

DGE <- DGEList(counts = toc)
DGE <- estimateGLMCommonDisp(DGE, design.matrix)
DGE <- estimateGLMTrendedDisp(DGE, design.matrix)
DGE <- estimateGLMTagwiseDisp(DGE, design.matrix)

fit <- glmFit(DGE, design.matrix)

###############
### C vs A B

###############

lrt_C_AB <- glmLRT(fit, contrast = c(-0.5, -0.5, 1, 0, 0))

exp_C_AB=data.frame(topTags(lrt_C_AB, n=nrow(toc)))

##############
### A vs B 

##############

lrt_A_B <- glmLRT(fit, contrast = c(1, -1, 0, 0, 0))

exp_A_B=data.frame(topTags(lrt_A_B, n=nrow(toc))) 

##############
### C vs B

##############

lrt_C_B <- glmLRT(fit, contrast = c(0, -1, 1, 0, 0))

exp_C_B=data.frame(topTags(lrt_C_B, n=nrow(toc))) 

###############
###  C vs A

###############

lrt_C_A <- glmLRT(fit, contrast = c(-1, 0, 1, 0, 0))

exp_C_A=data.frame(topTags(lrt_C_A, n=nrow(toc))) 

###################################################################################################################################

### Heatmap

##################################################################################################################################

########################################################

### ABC list 

########################################################

exp_ABC=rbind(data.frame(subset(exp_A_B, exp_A_B$FDR<=0.05)),
              data.frame(subset(exp_C_A, exp_C_A$FDR<=0.05)) ,
              data.frame(subset(exp_C_B, exp_C_B$FDR<=0.05)))  

exp_ABC=rbind(data.frame(subset(exp_A_B, exp_A_B$FDR<=0.5)),
              data.frame(subset(exp_C_A, exp_C_A$FDR<=0.5)) ,
              data.frame(subset(exp_C_B, exp_C_B$FDR<=0.5)))  

exp_ABC=rbind(data.frame(subset(exp_A_B, exp_A_B$FDR<=0.1)),
              data.frame(subset(exp_C_A, exp_C_A$FDR<=0.1)) ,
              data.frame(subset(exp_C_B, exp_C_B$FDR<=0.1)))  

exp_ABC=rbind(data.frame(subset(exp_A_B, exp_A_B$PValue<=0.05)),
              data.frame(subset(exp_C_A, exp_C_A$PValue<=0.05)) ,
              data.frame(subset(exp_C_B, exp_C_B$PValue<=0.05)))  

exp_ABC=rbind(data.frame(subset(exp_A_B, exp_A_B$FDR<=0.9)),
              data.frame(subset(exp_C_A, exp_C_A$FDR<=0.9)) ,
              data.frame(subset(exp_C_B, exp_C_B$FDR<=0.9)) )


#exp_ABC_all = rbind (exp_A_B,exp_A_C,exp_B_C) 
#exp_ABC_all=data.frame(cbind(rownames(exp_ABC_all), exp_ABC_all))
#exp_ABC_all_uni=data.frame(unique(exp_ABC_all$rownames.exp_ABC_all.))

########################################################


heatmap.file=data.frame(subset(toc, rownames(toc) %in% rownames(exp_ABC)))


f<-colorRampPalette(c("blue", "white", "firebrick"))
col.all <- f(256)

pdf("heatmap_LW_gene_0.05Pvalue.pdf")

finalmat <- scale(t(cpm(heatmap.file, prior.count=2, log=TRUE)))
finalmat[finalmat > 1] = 1
finalmat[finalmat < -1] = -1
hmap(as.matrix(finalmat),   cexRow = 0.8, cexCol = 0.2, margins = c(8,5), options = list (prop=TRUE, col = col.all, main="", cex.main=0.08))  # for FDR=0.5  cexCol = 0.2 

dev.off()


######################################################

#######################
### adding rank and FC

#######################

FC_C_AB=2^exp_C_AB$logFC
exp_C_AB=data.frame(cbind(exp_C_AB, FC_C_AB))

FC_A_B=2^exp_A_B$logFC
exp_A_B=data.frame(cbind(exp_A_B, FC_A_B)) 

FC_C_A=2^exp_C_A$logFC
exp_C_A=data.frame(cbind(exp_C_A, FC_C_A))

FC_C_B=2^exp_C_B$logFC
exp_C_B=data.frame(cbind(exp_C_B, FC_C_B))



#########
# Rank
#########

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
Rank_C_B=which(rownames(exp_C_B)==rownames(exp_C_AB)[i])
logFC_C_B=exp_C_B[Rank_C_B, 1]
}
else {
r=which(rownames(exp_C_B)==rownames(exp_C_AB)[i])
Rank_C_B=c(Rank_C_B,r)
l=exp_C_B[r, 1]
logFC_C_B=c(logFC_C_B,l)
FC_C_B=2^ logFC_C_B
}
}


for (i in 1:nrow(exp_C_AB)){ 
if (i==1){
Rank_C_A=which(rownames(exp_C_A)==rownames(exp_C_AB)[i])
logFC_C_A=exp_C_A[Rank_C_A, 1] 
}
else {
r=which(rownames(exp_C_A)==rownames(exp_C_AB)[i])
Rank_C_A=c(Rank_C_A,r)
l=exp_C_A[r, 1]
logFC_C_A=c(logFC_C_A,l)
FC_C_A=2^ logFC_C_A
}
}

exp_C_AB_rank=data.frame(cbind(exp_C_AB, Rank_A_B, Rank_C_A,Rank_C_B, FC_A_B,  FC_C_A,  FC_C_B))

######################################################

exp_C_AB=data.frame(cbind(rownames(exp_C_AB), exp_C_AB))
colnames(exp_C_AB)[1]="gene.ID"
rownames(exp_C_AB)=NULL

exp_C_A=data.frame(cbind(rownames(exp_C_A), exp_C_A))
colnames(exp_C_A)[1]="gene.ID"
rownames(exp_C_A)=NULL

exp_C_B=data.frame(cbind(rownames(exp_C_B), exp_C_B))
colnames(exp_C_B)[1]="gene.ID"
rownames(exp_C_B)=NULL

exp_A_B=data.frame(cbind(rownames(exp_A_B), exp_A_B))
colnames(exp_A_B)[1]="gene.ID"
rownames(exp_A_B)=NULL

whole.results=merge(exp_C_AB, exp_A_B ,by="gene.ID", all=T)
whole.results=merge(whole.results, exp_C_A,by="gene.ID", all=T)
whole.results=merge(whole.results, exp_C_B,by="gene.ID", all=T)

s_whole_results=whole.results[ ,c(1,7,5,6,13,11,12,19,17,18,25,23,24)]
colnames(s_whole_results)= c("gene..ID","FC..C_AB","PValue..C_AB","FDR..C_AB","FC..A_B", "PValue..A_B","FDR..A_B","FC..C_A","PValue..C_A","FDR..C_A","FC..C_B", "PValue..C_B","FDR..C_B")
s_whole_results=s_whole_results[order(s_whole_results$FDR..C_AB), ]
row.names(s_whole_results)=NULL