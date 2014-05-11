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
i#nstall.packages("limma")
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

exp_list=data.frame(topTags(lrt, n=nrow(count_norm_LW)))
exp_list_b=cbind(rownames(exp_list), exp_list, row.names=F)
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

##############################################################
###cluster the transcripts to seperate Hempglubin transcripts
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

#exp_list_annotate_clust_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clust_sort.txt", sep="\t", header=T)

### we want to exclude transcripts with cluster 2 from our main data set

exp_list_annotate_clust_sort_exclude=data.frame(subset(exp_list_annotate_clust_sort, exp_list_annotate_clust_sort$clusterID==2))
#412 9

count_norm_LW_clean=data.frame(subset(count_norm_LW, !(count_norm_LW$Transcript_ID  %in% exp_list_annotate_clust_sort_exclude$Ensembl.Transcript.ID)))
#13948 17

##########################################################

### repeatinfg the process for count_norm_LW_clean data set

###########################################################

### using edge R making contrast c-(A+B)

## read the below file instead of command
exp_list_clean_b=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_clean_b.txt", sep="\t", header=T)

count_norm_LW_clean_b=count_norm_LW_clean 
row.names(count_norm_LW_clean_b)=count_norm_LW_clean_b[, 1]
count_norm_LW_clean_b=count_norm_LW_clean_b[, -1]

D <- DGEList(counts=count_norm_LW_clean_b)
D <- estimateGLMCommonDisp(D, design)
D <- estimateGLMTagwiseDisp(D, design)
fit <- glmFit(D, design)
lrt <- glmLRT(fit, contrast = c(-0.5, -0.5, 1))

exp_list_clean=data.frame(topTags(lrt, n=nrow(count_norm_LW_clean))) #13948 5
exp_list_clean_b=cbind(rownames(exp_list_clean), exp_list_clean)
colnames(exp_list_clean_b)[1]="Ensembl.Transcript.ID"

### removing row anmes

row.names(exp_list_clean_b) <- NULL 

### writting the file

# write.table(exp_list_clean_b, "/Users/saeedehazary/Documents/LWRNASEQ/exp_list_clean_b.txt", sep="\t", row.names=F)


### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)
names(mart)
names(exp_list_clean_b)

exp_list_annotate_clean=merge(mart, exp_list_clean_b, by="Ensembl.Transcript.ID", all.exp_list_clean_b=all)

### sort

exp_list_annotate_clean_sort <- exp_list_annotate_clean[order(exp_list_annotate_clean$FDR),]

#write.table(exp_list_annotate_clean_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_sort.txt", sep="\t", row.names=F)


exp_list_annotate_clean_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_sort.txt", sep="\t", header=T)

########################################################

### cluster the samples

#######################################################


count_norm_LW_2=count_norm_LW


### kwwping the top 650 differentially expreesed genes

count_norm_LW_3=data.frame(subset(count_norm_LW_2, count_norm_LW_2$Transcript_ID %in% exp_list_annotate_clean_sort$Ensembl.Transcript.ID[1:650]))
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

count_norm_LW_2=count_norm_LW

###keep top 2000 transcripts 

count_norm_LW_4=data.frame(subset(count_norm_LW_2, count_norm_LW_2$Transcript_ID %in% exp_list_annotate_clean_sort$Ensembl.Transcript.ID[1:10]))
rownames(count_norm_LW_4)=count_norm_LW_4[,1]
count_norm_LW_4=count_norm_LW_4[, -1]  # 2000 16
count_norm_LW_5=as.matrix(count_norm_LW_4)
names(count_norm_LW_4)

myvars=c("DO297.3_eff_counts", "DO297.7_eff_counts", "DO307.1_eff_counts", "DO307.4_eff_counts", 
         "DO307.7_eff_counts", "DO310.4_eff_counts", "DO297.2_eff_counts", "DO297.6_eff_counts",
         "DO307.5_eff_counts", "DO310.2_eff_counts", "DO310.5_eff_counts", "DO307.2_eff_counts",
         "DO307.3_eff_counts", "DO307.6_eff_counts", "DO310.3_eff_counts", "DO310.1_eff_counts")
count_norm_LW_6=count_norm_LW_4[, myvars]

count_norm_LW_7=as.matrix(count_norm_LW_6)

count_heatmap<- heatmap(count_norm_LW_7, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))


t_count_norm_LW_7=t(count_norm_LW_7)
count_heatmap<- heatmap(t_count_norm_LW_7, Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(1,1))


############################################################

install.packages("gplots")
library(gplots)

count_heatmap<- heatmap.2(t_count_norm_LW_7, colsep=c(1:6),rowsep=(1:62),  sepwidth=c(0.05,0.05), sepcolor="white", trace="none",
                        Rowv=T,Colv=T, scale="none", dendrogram="both",key=F, 
                        lhei = c(0.05,5),margins=c(1,8))



############################################################











