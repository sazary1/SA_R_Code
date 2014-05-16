##################################################################################################################################
# SA 
# May 8 2014
# LWRNASEQ analysis
###################################################################################################################################

setwd ("/Users/saeedehazary/Documents/LWRNASEQ")


###################################################################################################################################

### Creating list of IDs and combining the outputs of eXpress

###################################################################################################################################

LW_ID=read.table("/Users/saeedehazary/Documents/LWRNASEQ/LW_ID.txt", sep="\t", header=F)
LW_ID
colnames(LW_ID)[1]="ID"

for(j in LW_ID$ID ) { 
  if (j=="DO297-2"){
    a= sprintf("/Users/saeedehazary/Documents/LWRNASEQ/Sample_%s/results.xprs", j) 
    b=read.table(a, header=T, sep="\t" )
    colnames(b) [2]="Transcript_name"
    colnames(b) [8]=sprintf("%s_eff_counts",j)
    count_LW=b[ ,c(2,8)]
  } else{
    a= sprintf("/Users/saeedehazary/Documents/LWRNASEQ/Sample_%s/results.xprs", j) 
    b=read.table(a, header=T, sep="\t" )
    colnames(b) [2]="Transcript_name"
    colnames(b) [8]=sprintf("%s_eff_counts",j)
    b=b[ ,c(2,8)]
    
    count_LW<- merge(count_LW,b,by="Transcript_name", all=TRUE)
    print(j)
  }
  
}

dim(count_LW)
#82934    17


colnames(count_LW)=c("Transcript_name", "B_DO297.2" , "C_DO297.3" , "B_DO297.6", "C_DO297.7" , "C_DO307.1" , "A_DO307.2"
                     , "A_DO307.3" , "C_DO307.4" ,"B_DO307.5" , "A_DO307.6" , "C_DO307.7" , "A_DO310.1"
                     , "B_DO310.2" , "A_DO310.3" ,"C_DO310.4" , "B_DO310.5")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

### Write and READ

write.table(count_LW, file="/Users/saeedehazary/Documents/LWRNASEQ/count_LW.txt", sep="\t", row.names=F)
count_LW=read.table("/Users/saeedehazary/Documents/LWRNASEQ/count_LW.txt", sep="\t", header=T)
#82934 17

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

###################################################################################################################################

# Normalizing the data

###################################################################################################################################

#install.packages("edgeR")
#nstall.packages("limma")
library(limma)
library(edgeR)
#install.packages("cpm")
library(cpm)


count_LW_all <- count_LW
count_15_LW <- count_LW_all[,2:ncol(count_LW_all)]
rownames(count_15_LW) <- count_LW_all[,1]
keep <- rowSums(cpm(count_15_LW)>5)>=10
toc <- count_15_LW[keep,]
table(keep)
norm_factors <- calcNormFactors(as.matrix(toc), method = "TMM")
#1:16

colsum=colSums(toc, na.rm=TRUE)
scalefactor <- mean(colsum)/colsum
print(scalefactor)
#1:16

count_norm_LW <- data.frame(sweep(as.matrix(toc), MARGIN = 2, norm_factors*scalefactor, '*'))
#14360 16

count_norm_LW=cbind(rownames(count_norm_LW),count_norm_LW , row.names=NULL)
dim(count_norm_LW)
#14360 17

colnames(count_norm_LW)[1]="Transcript_ID" 

colnames(count_norm_LW)=c("Transcript_ID", "B_DO297.2" , "C_DO297.3" , "B_DO297.6", "C_DO297.7" , "C_DO307.1" , "A_DO307.2"
                          , "A_DO307.3" , "C_DO307.4" ,"B_DO307.5" , "A_DO307.6" , "C_DO307.7" , "A_DO310.1"
                          , "B_DO310.2" , "A_DO310.3" ,"C_DO310.4" , "B_DO310.5")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

### Write and READ

write.table(count_norm_LW, file="/Users/saeedehazary/Documents/LWRNASEQ/count_norm_LW.txt", sep="\t", row.names=F)
count_norm_LW=read.table( "/Users/saeedehazary/Documents/LWRNASEQ/count_norm_LW.txt", sep="\t",  header=T)
# 14360 17

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

###################################################################################################################################

### Creating the design matrix for row data

##################################################################################################################################

##############
#desigm matrix
#############
Group <- c('B','C','B','C','C','A','A','C','B','A','C','A','B','A','C','B')
Litter <- c('DO297','DO297','DO297','DO297','DO307','DO307','DO307','DO307','DO307','DO307','DO307','DO310','DO310','DO310','DO310','DO310')
design <- model.matrix(~ 0 + Group + Litter)
print(design)

###################################################################################################################################

### using edge R making contrast c-(A+B)  nonclean data

###################################################################################################################################
library(edgeR)
library(limma)


count_norm_LW_b=count_norm_LW 
row.names(count_norm_LW_b)=count_norm_LW_b[, 1]
count_norm_LW_b=count_norm_LW_b[, -1]

count_norm_LW_b.int <- apply(count_norm_LW_b, 2, as.integer)
rownames(count_norm_LW_b.int) <- rownames(count_norm_LW_b)
count_norm_LW_b <- count_norm_LW_b.int

D_raw <- DGEList(counts=count_norm_LW_b)
D_raw <- estimateGLMCommonDisp(D_raw , design)
D_raw <- estimateGLMTrendedDisp(D_raw, design)
D_raw <- estimateGLMTagwiseDisp(D_raw , design)


fit_raw_old<- glmFit(D_raw, design)


lrt_raw_old <- glmLRT(fit_raw_old, contrast = c(-0.5, -0.5, 1, 0, 0)) 

exp_T_A_BC=data.frame(topTags(lrt_raw_old, n=nrow(count_norm_LW)))
exp_T_A_BC_2=cbind(rownames(exp_T_A_BC),  exp_T_A_BC)
colnames(exp_T_A_BC_2)[1]="Ensembl.Transcript.ID" 
rownames(exp_T_A_BC_2)=NULL



##########################################################

### Adding gene symbols, we get the mart file from Ensembl

#########################################################

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_T_A_BC_an=merge(mart, exp_T_A_BC_2, by="Ensembl.Transcript.ID", all.y=T)

### sort

exp_T_A_BC_an_s<- exp_T_A_BC_an[order(exp_T_A_BC_an$FDR),]


################

### A vs B clean  with out cleaning 

#################

lrt_A_B_old <- glmLRT(fit_raw_old, contrast = c(1, -1, 0, 0, 0))

exp_T_A_B=data.frame(topTags(lrt_A_B_old, n=nrow(count_norm_LW))) #14382 5
exp_T_A_B_2=cbind(rownames(exp_T_A_B), exp_T_A_B)
colnames(exp_T_A_B_2)[1]="Ensembl.Transcript.ID"
row.names(exp_T_A_B_2) <- NULL


### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_T_A_B_an=merge(mart, exp_T_A_B_2, by="Ensembl.Transcript.ID", all.exp_T_A_B_2=all)

### sort

exp_T_A_B_an_s<- exp_T_A_B_an[order(exp_T_A_B_an$FDR),] #14384 8


##################

### B vas C clean 

##################

lrt_B_C_old <- glmLRT(fit_raw_old, contrast = c(0, 1, -1, 0, 0))

exp_T_B_C=data.frame(topTags(lrt_B_C_old, n=nrow(count_norm_LW)))
exp_T_B_C_2=cbind(rownames(exp_T_B_C), exp_T_B_C)
colnames(exp_T_B_C_2)[1]="Ensembl.Transcript.ID"
row.names(exp_T_B_C_2) <- NULL


### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_T_B_C_an=merge(mart,exp_T_B_C_2, by="Ensembl.Transcript.ID", all.exp_T_B_C_2=all)

### sort

exp_T_B_C_an_s <- exp_T_B_C_an[order(exp_T_B_C_an$FDR),]


##################

### A vs C clean

##################

lrt_A_C_old<- glmLRT(fit_raw_old, contrast = c(1, 0, -1, 0, 0))

exp_T_A_C=data.frame(topTags(lrt_A_C_old, n=nrow(count_norm_LW))) #13948 5
exp_T_A_C_2=cbind(rownames(exp_T_A_C), exp_T_A_C)
colnames(exp_T_A_C_2)[1]="Ensembl.Transcript.ID"
row.names(exp_T_A_C_2) <- NULL


### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_T_A_C_an=merge(mart, exp_T_A_C_2, by="Ensembl.Transcript.ID", all.exp_T_A_C_2=all)

### sort

exp_T_A_C_an_s <- exp_T_A_C_an[order(exp_T_A_C_an$FDR),]


###################################################################################################################################

### Heatmap

###################################################################################################################################


#########################################################

### ABC list

########################################################

exp_T_ABC=rbind(data.frame(subset(exp_T_A_B_an_s, exp_T_A_B_an_s$FDR<=0.05)),
                              data.frame(subset(exp_T_B_C_an_s, exp_T_B_C_an_s$FDR<=0.05)) ,
                              data.frame(subset(exp_T_A_C_an_s, exp_T_B_C_an_s$FDR<=0.05)))  #

goodlist_T <- data.frame(unique(exp_T_ABC$Ensembl.Transcript.ID)) #
#goodlist <- unique(exp_annot_clean_ABC$Associated.Gene.Name)

########################################################


colnames(goodlist_T)[1]="Ensembl.Transcript.ID"

heatmap.file.T=data.frame(subset(count_norm_LW, count_norm_LW$Transcript_ID %in% goodlist_T$Ensembl.Transcript.ID))


rownames(heatmap.file.T)=heatmap.file.T[,1]
heatmap.file.T=heatmap.file.T[, -1] 

install.packages("gplots")
library(gplots)
library(seriation)

f<-colorRampPalette(c("blue", "white", "firebrick"))
col.all <- f(256)

pdf("heatmap_litter_adj.pdf")

finalmat <- scale(t(cpm(heatmap.file.T, prior.count=2, log=TRUE)))
#finalmat <- scale(t(log(count_norm_LW_clean_4, base=2)))
finalmat[finalmat > 1] = 1
finalmat[finalmat < -1] = -1
count_heatmap<- hmap(as.matrix(finalmat),   cexRow = 0.8, cexCol = 0.005, margins = c(10,5), options = list (prop=TRUE, col = col.all, main="litter.adjustment", cex.main=0.5))

dev.off()
