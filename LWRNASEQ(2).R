######################################################################
# SA 
# May 8 2014
# LWRNASEQ analysis
#####################################################################

setwd ("/Users/saeedehazary/Documents/LWRNASEQ")


#### Creating list of IDs and combining the outputs of eXpress ##############################
# To avoid running the below for loop you can read the file 
count_LW=read.table("/Users/saeedehazary/Documents/LWRNASEQ/count_LW.txt", sep="\t", header=T)
###############################################################################################

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

#write.table(count_LW, file="/Users/saeedehazary/Documents/LWRNASEQ/count_LW.txt", sep="\t", row.names=F)

########################################################################################################

# Creating norm factor and scale factor using edge R and applying those on the count_LW data set

########################################################################################################
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
write.table(count_norm_LW, file="/Users/saeedehazary/Documents/LWRNASEQ/count_norm_LW.txt", sep="\t", row.names=F)

#########################################################################################################

### Creating the design matrix 

A=c(0,0,0,0,0,1,1,0,0,1,0,1,0,1,0,0)
B=c(1,0,1,0,0,0,0,0,1,0,0,0,1,0,0,1)
C=c(0,1,0,1,1,0,0,1,0,0,1,0,0,0,1,0)
litter=c("B","C","B","C","C","A","A", "C", "B", "A", "C", "A", "B", "A", "C", "B")
LW_ID=read.table("/Users/saeedehazary/Documents/LWRNASEQ/LW_ID.txt", sep="\t", header=F)
design=data.frame(cbind(LW_ID,A,B,C))
row.names(design)=design[,1]
design=design[,-1]
design_matrix=as.matrix(design)
design_matrix

########################################################################################################

### using edge R making contrast c-(A+B)

## instead of runnig below code read the file 
exp_list_b=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_b.txt", sep="\t", header=T)

count_norm_LW_b=count_norm_LW
row.names(count_norm_LW_b)=count_norm_LW_b[, 1]
count_norm_LW_b=count_norm_LW_b[, -1]

D <- DGEList(counts=count_norm_LW_b)
D <- estimateGLMCommonDisp(D, design)
D <- estimateGLMTagwiseDisp(D, design)
fit <- glmFit(D, design)
lrt <- glmLRT(fit, contrast = c(-0.5, -0.5, 1))
#lrt <- glmLRT(fit, contrast = c(-1, 1, 0))

exp_list=data.frame(topTags(lrt, n=nrow(count_norm_LW)))
exp_list_b=cbind(rownames(exp_list), exp_list)
rownames(exp_list_b)=NULL
colnames(exp_list_b)[1]="Ensembl.Transcript.ID"
#write.table(exp_list_b,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_b.txt", sep="\t", row.names=F)

#######################################################################################################

### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)
names(mart)
names(exp_list_b)

exp_list_annotate=merge(mart, exp_list_b, by="Ensembl.Transcript.ID", all.exp_list_b=all)

### sort

exp_list_annotate_sort <- exp_list_annotate[order(exp_list_annotate$FDR),]
write.table(exp_list_annotate_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_sort.txt", sep="\t", row.names=F)

##############################################################
### cluster the transcripts to seperate Hempglubin transcripts
##############################################################

# subset the count_norm_LW data set

sub_count_norm_LW=data.frame(subset(count_norm_LW, count_norm_LW$Transcript_ID %in% exp_list_annotate_sort$Ensembl.Transcript.ID[1:4000]))
#4000 17


#sub_count_norm_LW_all=data.frame(subset(count_norm_LW, count_norm_LW$Transcript_ID %in% exp_list_annotate_sort$Ensembl.Transcript.ID))

### as.matrix

sub_count_matrix <- as.matrix(sub_count_norm_LW[,2:ncol(sub_count_norm_LW)])

### Create row.names

rownames(sub_count_matrix) <- sub_count_norm_LW[,1]

### Run cluster command

library(mclust)

clust <- Mclust(sub_count_matrix, G=2)

### Make a data frame with two columns col1 has Ensembl.Transcript.ID and col2 has clusterID

clustsubset=cbind(rownames(data.frame(clust$classification)), data.frame(clust$classification))
colnames(clustsubset) <- c("Ensembl.Transcript.ID","clusterID")

### Adding cluster ID to exp_list_annotate_sort

exp_list_annotate_clust <- merge(exp_list_annotate_sort, clustsubset, by = "Ensembl.Transcript.ID", all.exp_list_annotate_sort=exp_list_annotate_sort)

exp_list_annotate_clust_sort <- exp_list_annotate_clust[order(exp_list_annotate_clust$FDR),]

#write.table(exp_list_annotate_clust_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clust_sort.txt", sep="\t", row.names=F)

exp_list_annotate_clust_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clust_sort.txt", sep="\t", header=T)

all_14360_mgi_symbol_list=data.frame(exp_list_annotate_clust_sort$Ensembl.Transcript.ID)
write.table(all_14360_mgi_symbol_list, file="/Users/saeedehazary/Documents/LWRNASEQ/all_14360_mgi_symbol_list.txt", sep="/t", row.names=F)

#############################################################
### exclude transcripts with cluster 2 from our main data set
#############################################################

exp_list_annotate_clust_sort_exclude=data.frame(subset(exp_list_annotate_clust_sort, exp_list_annotate_clust_sort$clusterID==2))
#412 9

Hgb_mgi_symbol_list=data.frame(exp_list_annotate_clust_sort_exclude$Ensembl.Transcript.ID)
write.table(Hgb_mgi_symbol_list, file="/Users/saeedehazary/Documents/LWRNASEQ/Hgb_mgi_symbol_list", sep="/t", row.names=F)

############################
### creating a clean data set
############################

count_LW_clean=data.frame(subset(count_LW, !(count_LW$Transcript_name  %in% exp_list_annotate_clust_sort_exclude$Ensembl.Transcript.ID)))
#82522 17

##########################################################

### repeatinfg the process for count_LW_clean data set

###########################################################

### Normalizing the clean data again

library(limma)
library(edgeR)
#install.packages("cpm")
library(cpm)


count_LW_clean_all <- count_LW_clean
count_16_clean_LW <- count_LW_clean_all[,2:ncol(count_LW_clean_all)]
rownames(count_16_clean_LW) <- count_LW_clean_all[,1]
keep_clean <- rowSums(cpm(count_16_clean_LW)>5)>=10
toc_clean <- count_16_clean_LW[keep_clean,]
table(keep_clean)
norm_factors_clean <- calcNormFactors(as.matrix(toc_clean), method = "TMM")
#1:16

colsum=colSums(toc_clean, na.rm=TRUE)
scalefactor_clean <- mean(colsum)/colsum
print(scalefactor_clean)
#1:16

count_norm_LW_clean <- data.frame(sweep(as.matrix(toc_clean), MARGIN = 2, norm_factors_clean*scalefactor_clean, '*'))
#14384 16

count_norm_LW_clean=cbind(rownames(count_norm_LW_clean),count_norm_LW_clean)
rownames(count_norm_LW_clean)=NULL
dim(count_norm_LW_clean)
#14384 17

colnames(count_norm_LW_clean)[1]="Transcript_ID" 

#########################################################################################################

### Creating the design matrix 

A=c(0,0,0,0,0,1,1,0,0,1,0,1,0,1,0,0)
B=c(1,0,1,0,0,0,0,0,1,0,0,0,1,0,0,1)
C=c(0,1,0,1,1,0,0,1,0,0,1,0,0,0,1,0)
LW_ID=read.table("/Users/saeedehazary/Documents/LWRNASEQ/LW_ID.txt", sep="\t", header=F)
design=data.frame(cbind(LW_ID,A,B,C))
row.names(design)=design[,1]
design=design[,-1]
design_matrix=as.matrix(design)
design_matrix


count_norm_LW_clean_b=count_norm_LW_clean 
row.names(count_norm_LW_clean_b)=count_norm_LW_clean_b[, 1]
count_norm_LW_clean_b=count_norm_LW_clean_b[, -1]

D_clean <- DGEList(counts=count_norm_LW_clean_b)
D_clean <- estimateGLMCommonDisp(D_clean, design)
D_clean <- estimateGLMTagwiseDisp(D_clean, design)
fit_clean <- glmFit(D_clean, design)

##############
### c - (A+B)
##############

lrt_C_AB <- glmLRT(fit_clean, contrast = c(-0.5, -0.5, 1))

exp_list_clean_C_AB=data.frame(topTags(lrt_C_AB, n=nrow(count_norm_LW_clean))) 
exp_list_clean_C_AB_2=cbind(rownames(exp_list_clean_C_AB), exp_list_clean_C_AB)
colnames(exp_list_clean_C_AB_2)[1]="Ensembl.Transcript.ID"
row.names(exp_list_clean_C_AB_2) <- NULL

### writting the file
#write.table(exp_list_clean_C_AB_2, "/Users/saeedehazary/Documents/LWRNASEQ/exp_list_clean_C_AB_2.txt", sep="\t", row.names=F)


### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_list_annotate_clean_C_AB=merge(mart, exp_list_clean_C_AB_2, by="Ensembl.Transcript.ID", all.exp_list_clean_C_AB_2=all)
names( count_norm_LW)
count_norm_LW_cc=count_norm_LW
colnames(count_norm_LW_cc)[1]="Ensembl.Transcript.ID"
count_norm_LW_annotated=merge(mart, count_norm_LW_cc, by="Ensembl.Transcript.ID", all.count_norm_LW=all)
write.table(count_norm_LW_annotated,file="/Users/saeedehazary/Documents/LWRNASEQ/count_norm_LW_annotated.txt", sep="\t", row.names=F)

### sort

exp_list_annotate_clean_C_AB_sort <- exp_list_annotate_clean_C_AB[order(exp_list_annotate_clean_C_AB$FDR),]


#write.table(exp_list_annotate_clean_C_AB_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_C_AB_sort.txt", sep="\t", row.names=F)
#exp_list_annotate_clean_C_AB_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_C_AB_sort.txt", sep="\t", header=T)


##############
### A vs B
##############

lrt_A_B <- glmLRT(fit_clean, contrast = c(1, -1, 0))

exp_list_clean_A_B=data.frame(topTags(lrt_A_B, n=nrow(count_norm_LW_clean))) #14382 5
exp_list_clean_A_B_2=cbind(rownames(exp_list_clean_A_B), exp_list_clean_A_B)
colnames(exp_list_clean_A_B_2)[1]="Ensembl.Transcript.ID"
row.names(exp_list_clean_A_B_2) <- NULL

### writting the file
# write.table(exp_list_clean_A_B_2, "/Users/saeedehazary/Documents/LWRNASEQ/exp_list_clean_A_B_2.txt", sep="\t", row.names=F)


### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_list_annotate_clean_A_B=merge(mart, exp_list_clean_A_B_2, by="Ensembl.Transcript.ID", all.exp_list_clean_A_B_2=all)

### sort

exp_list_annotate_clean_A_B_sort <- exp_list_annotate_clean_A_B[order(exp_list_annotate_clean_A_B$FDR),] #14384 8

#write.table(exp_list_annotate_clean_A_B_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_A_B_sort.txt", sep="\t", row.names=F)
#exp_list_annotate_clean_A_B_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_A_B_sort.txt", sep="\t", header=T)


##################
### B vas C
###################

lrt_B_C <- glmLRT(fit_clean, contrast = c(0, 1, -1))

exp_list_clean_B_C=data.frame(topTags(lrt_B_C, n=nrow(count_norm_LW_clean))) #13948 5
exp_list_clean_B_C_2=cbind(rownames(exp_list_clean_B_C), exp_list_clean_B_C)
colnames(exp_list_clean_B_C_2)[1]="Ensembl.Transcript.ID"
row.names(exp_list_clean_B_C_2) <- NULL

### writting the file
# write.table(exp_list_clean_B_C_2, "/Users/saeedehazary/Documents/LWRNASEQ/exp_list_clean_B_C_2.txt", sep="\t", row.names=F)


### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_list_annotate_clean_B_C=merge(mart, exp_list_clean_B_C_2, by="Ensembl.Transcript.ID", all.exp_list_clean_B_C_2=all)

### sort

exp_list_annotate_clean_B_C_sort <- exp_list_annotate_clean_B_C[order(exp_list_annotate_clean_B_C$FDR),]

#write.table(exp_list_annotate_clean_B_C_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_B_C_sort.txt", sep="\t", row.names=F)
#exp_list_annotate_clean_B_C_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_B_C_sort.txt", sep="\t", header=T)

##################
### A vs C
###################

lrt_A_C <- glmLRT(fit_clean, contrast = c(1, 0, -1))

exp_list_clean_A_C=data.frame(topTags(lrt_A_C, n=nrow(count_norm_LW_clean))) #13948 5
exp_list_clean_A_C_2=cbind(rownames(exp_list_clean_A_C), exp_list_clean_A_C)
colnames(exp_list_clean_A_C_2)[1]="Ensembl.Transcript.ID"
row.names(exp_list_clean_A_C_2) <- NULL

### writting the file
# write.table(exp_list_clean_A_C_2, "/Users/saeedehazary/Documents/LWRNASEQ/exp_list_clean_A_C_2.txt", sep="\t", row.names=F)

### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_list_annotate_clean_A_C=merge(mart, exp_list_clean_A_C_2, by="Ensembl.Transcript.ID", all.exp_list_clean_A_C_2=all)

### sort

exp_list_annotate_clean_A_C_sort <- exp_list_annotate_clean_A_C[order(exp_list_annotate_clean_A_C$FDR),]

#write.table(exp_list_annotate_clean_A_C_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_A_C_sort.txt", sep="\t", row.names=F)
#exp_list_annotate_clean_A_C_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_A_C_sort.txt", sep="\t", header=T)

########################################################

### cluster the samples

#######################################################


count_norm_LW_2=count_norm_LW


### kwwping the top 650 differentially expreesed genes

count_norm_LW_3=data.frame(subset(count_norm_LW_2, count_norm_LW_2$Transcript_ID %in% exp_list_annotate_clean_C_AB_sort$Ensembl.Transcript.ID[1:650]))
rownames(count_norm_LW_3)=NULL

rownames(count_norm_LW_3)=count_norm_LW_3[,1]
count_norm_LW_3=count_norm_LW_3[, -1]

### transpose the data set

t_count_norm_LW_3=t(count_norm_LW_3)

d <- dist(as.matrix(t_count_norm_LW_3))   # find distance matrix 
hc <- hclust(d)                # apply hirarchical clustering 
plot(hc) 


groups <- cutree(hc, k=3) # cut tree into 3 clusters

### draw dendogram with red borders around the 3 clusters 

rect.hclust(hc, k=3, border="red")


#########################################################

### Heatmap

########################################################

exp_list_annotate_clean_ABC=rbind(data.frame(subset(exp_list_annotate_clean_A_B, exp_list_annotate_clean_A_B$FDR<=0.05)),
                                  data.frame(subset(exp_list_annotate_clean_A_C, exp_list_annotate_clean_A_C$FDR<=0.05)) ,
                                  data.frame(subset(exp_list_annotate_clean_B_C, exp_list_annotate_clean_B_C$FDR<=0.05)))

goodlist <- unique(exp_list_annotate_clean_ABC$Ensembl.Transcript.ID)
goodlist <- unique(exp_list_annotate_clean_ABC$Associated.Gene.Name)
count_norm_LW_clean_2=count_norm_LW_clean

count_norm_LW_clean_4=data.frame(subset(count_norm_LW_clean_2, count_norm_LW_clean_2$Transcript_ID %in% exp_list_annotate_clean_A_B_sort$Ensembl.Transcript.ID[1:500]))

count_norm_LW_clean_4=data.frame(subset(count_norm_LW_clean_2, count_norm_LW_clean_2$Transcript_ID %in% goodlist))

rownames(count_norm_LW_clean_4)=count_norm_LW_clean_4[,1]
count_norm_LW_clean_4=count_norm_LW_clean_4[, -1]  # 2000 16
#count_norm_LW_clean_5=as.matrix(count_norm_LW_clean_4)
#names(count_norm_LW_clean_5)

myvars=c("DO297.3_eff_counts", "DO297.7_eff_counts", "DO307.1_eff_counts", "DO307.4_eff_counts", 
         "DO307.7_eff_counts", "DO310.4_eff_counts", "DO297.2_eff_counts", "DO297.6_eff_counts",
         "DO307.5_eff_counts", "DO310.2_eff_counts", "DO310.5_eff_counts", "DO307.2_eff_counts",
         "DO307.3_eff_counts", "DO307.6_eff_counts", "DO310.3_eff_counts", "DO310.1_eff_counts")

myvars_A_B=c( "DO297.2_eff_counts", "DO297.6_eff_counts",
         "DO307.5_eff_counts", "DO310.2_eff_counts", "DO310.5_eff_counts", "DO307.2_eff_counts",
         "DO307.3_eff_counts", "DO307.6_eff_counts", "DO310.3_eff_counts", "DO310.1_eff_counts")

myvars_B_C=c("DO297.3_eff_counts", "DO297.7_eff_counts", "DO307.1_eff_counts", "DO307.4_eff_counts", 
         "DO307.7_eff_counts", "DO310.4_eff_counts", "DO297.2_eff_counts", "DO297.6_eff_counts",
         "DO307.5_eff_counts", "DO310.2_eff_counts", "DO310.5_eff_counts")

myvars_C_A=c("DO297.3_eff_counts", "DO297.7_eff_counts", "DO307.1_eff_counts", "DO307.4_eff_counts", 
         "DO307.7_eff_counts", "DO310.4_eff_counts", "DO307.2_eff_counts",
         "DO307.3_eff_counts", "DO307.6_eff_counts", "DO310.3_eff_counts", "DO310.1_eff_counts")

count_norm_LW_clean_6=count_norm_LW_clean_4[, myvars_A_B]

count_norm_LW_clean_7=as.matrix(count_norm_LW_clean)

count_heatmap<- heatmap.2(as.matrix(cpm(count_norm_LW_clean_7, prior.count=2, log=TRUE)), col = heat.colors(256), scale="column", main="A vs B Transcripts by Sample ID")


t_count_norm_LW_clean_7=t(count_norm_LW_clean_7)
count_heatmap<- heatmap(t_count_norm_LW_clean_7, col = heat.colors(256), scale="row", margins=c(10,10),  main="A vs B Transcripts by Sample ID")


############################################################

install.packages("gplots")
library(gplots)

pdf("test_heatmap.pdf")
count_heatmap<- heatmap.2(as.matrix(cpm(count_norm_LW_clean_4, prior.count=2, log=TRUE)), margins=c(10,10), colsep=c(1:6),rowsep=(1:62),  sepwidth=c(0.05,0.05), sepcolor="white", trace="none",
                        Rowv=T,Colv=T, scale="row", dendrogram="both",key=F)

dev.off()
f<-colorRampPalette(c("blue", "white", "firebrick"))
col.all <- f(256)
pdf("test_heatmap_seriate.pdf")
finalmat <- scale(t(cpm(count_norm_LW_clean_4, prior.count=2, log=TRUE)))
finalmat[finalmat > 1] = 1
finalmat[finalmat < -1] = -1
count_heatmap<- hmap(as.matrix(finalmat),   cexRow = 0.8, cexCol = 0.00005, margins = c(7,7), options = list(prop=TRUE, col = col.all))
dev.off()
############################################################











