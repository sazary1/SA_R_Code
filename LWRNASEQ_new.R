##################################################################################################################################
# SA 
# May 8 2014
# LWRNASEQ analysis
###################################################################################################################################

setwd ("/Users/saeedehazary/Documents/LWRNASEQ")


###################################################################################################################################

### Creating list of IDs and combining the outputs of eXpress

###################################################################################################################################

LW_ID=read.table("/Users/saeedehazary/Documents/LWRNASEQ/LW_ID.txt", sep="\t", header=T)
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

# Normalizing teh data

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

###################################################################################################################################


#############
# old version
#############

A=c(0,0,0,0,0,1,1,0,0,1,0,1,0,1,0,0)
B=c(1,0,1,0,0,0,0,0,1,0,0,0,1,0,0,1)
C=c(0,1,0,1,1,0,0,1,0,0,1,0,0,0,1,0)
LW_ID=read.table("/Users/saeedehazary/Documents/LWRNASEQ/LW_ID.txt", sep="\t", header=F)
design_old=data.frame(cbind(LW_ID,A,B,C))
row.names(design_old)=design_old[,1]
design_old=design_old[,-1]
design_old=as.matrix(design_old)
design_old

#############
# New version
#############

geno <- factor(c("B","C","B","C","C","A","A","C","B","A","C","A","B","A","C","B"), levels = c("A","B","C"))
litter=factor(c(rep("DO297",4),rep("DO307",7), rep("DO310", 5)), levels = c("DO297", "DO307", "DO310"))
ident <- data.frame(geno, litter)
design <- model.matrix(~ -1 + geno + litter, ident)
design


###################################################################################################################################

### using edge R making contrast c-(A+B)  nonclean data

###################################################################################################################################

library(limma)
library(edgeR)


count_norm_LW_b=count_norm_LW
row.names(count_norm_LW_b)=count_norm_LW_b[, 1]
count_norm_LW_b=count_norm_LW_b[, -1]
#14360 16

D_raw <- DGEList(counts=count_norm_LW_b)
D_raw <- estimateGLMCommonDisp(D_raw , design)
D_raw <- estimateGLMTagwiseDisp(D_raw , design)
#D_raw <- estimateGLMCommonDisp(D_raw , design_old)
#D_raw <- estimateGLMTagwiseDisp(D_raw , design_old)


fit_raw <- glmFit(D_raw, design)
#fit_raw <- glmFit(D_raw, design_old)
lrt_raw <- glmLRT(fit_raw, contrast = c(-0.5, -0.5, 1,0,0))  ## Error in glmfit$coefficients %*% contrast : non-conformable arguments 

exp_list=data.frame(topTags(lrt_raw, n=nrow(count_norm_LW)))
exp_list_b=cbind(rownames(exp_list), exp_list)
rownames(exp_list_b)=NULL
colnames(exp_list_b)[1]="Ensembl.Transcript.ID" 


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

### Write and READ

write.table(exp_list_b,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_b.txt", sep="\t", row.names=F)
exp_list_b=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_b.txt", sep="\t", header=T)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

##########################################################

### Adding gene symbols, we get the mart file from Ensembl

#########################################################

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_list_annotate=merge(mart, exp_list_b, by="Ensembl.Transcript.ID", all.exp_list_b=all)

### sort

exp_list_annotate_sort <- exp_list_annotate[order(exp_list_annotate$FDR),]

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

### Write and READ

write.table(exp_list_annotate_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_sort.txt", sep="\t", row.names=F)
exp_list_annotate_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_sort.txt", sep="\t", header=T)  

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




###################################################################################################################################


### RBC data 

###################################################################################################################################

SRR1106109.xprs=read.table("/Users/saeedehazary/Documents/LWRNASEQ/Erythro/SRR1106109/results.xprs", header=T, sep="\t" )
colnames(SRR1106109.xprs) [2]="Transcript_name"
colnames(SRR1106109.xprs) [8]="SRR1106109_eff_counts"
count_SRR1106109=SRR1106109.xprs[ ,c(2,8)]

SRR1106108.xprs=read.table("/Users/saeedehazary/Documents/LWRNASEQ/Erythro/SRR1106108/results.xprs", header=T, sep="\t" )
colnames(SRR1106108.xprs) [2]="Transcript_name"
colnames(SRR1106108.xprs) [8]="SRR1106108_eff_counts"
count_SRR1106108=SRR1106108.xprs[ ,c(2,8)]

SRR1106110.xprs=read.table("/Users/saeedehazary/Documents/LWRNASEQ/Erythro/SRR1106110/results.xprs", header=T, sep="\t" )
colnames(SRR1106110.xprs) [2]="Transcript_name"
colnames(SRR1106110.xprs) [8]="SRR1106110_eff_counts"
count_SRR1106110=SRR1106110.xprs[ ,c(2,8)]


count_erythro<- merge(count_SRR1106109, count_SRR1106108, by="Transcript_name", all=TRUE)
count_erythro<- merge(count_erythro, count_SRR1106110, by="Transcript_name", all=TRUE)

### add cases to this data set

myvars<- c("Transcript_name","C_DO297.3","C_DO297.7","C_DO307.1","C_DO307.4","C_DO307.7","C_DO310.4")
count_LW_cases=count_LW[,myvars] 

count_erythro_cases= merge(count_erythro, count_LW_cases, by="Transcript_name", all=TRUE)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#### Write and Read

write.table(count_erythro_cases, file="/Users/saeedehazary/Documents/LWRNASEQ/count_erythro_cases.txt", sep="\t", row.names=F)
count_erythro_cases=read.table("/Users/saeedehazary/Documents/LWRNASEQ/count_erythro_cases.txt", sep="\t", header=T)     # 82934  9

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

####################################

### Normalizing the data RBC

####################################

library(limma)
library(edgeR)
library(cpm)


count_erythro_cases_copy <- count_erythro_cases
count_erythro_cases_copy_noID <- count_erythro_cases_copy[,2:ncol(count_erythro_cases_copy)]
rownames(count_erythro_cases_copy_noID) <- count_erythro_cases_copy[,1]
keep_erythro <- rowSums(cpm(count_erythro_cases_copy_noID)>5)>=4
toc_erythro <- count_erythro_cases_copy_noID[keep_erythro,]
table(keep_erythro)
norm_factors_erythro <- calcNormFactors(as.matrix(toc_erythro), method = "TMM")

colsum_erythro=colSums(toc_erythro, na.rm=TRUE)
scalefactor_erythro <- mean(colsum_erythro)/colsum_erythro
print(scalefactor_erythro)

count_erythro_cases_norm <- data.frame(sweep(as.matrix(toc_erythro), MARGIN = 2, norm_factors_erythro*scalefactor_erythro, '*'))

count_erythro_cases_norm <-cbind(rownames(count_erythro_cases_norm),count_erythro_cases_norm , row.names=NULL)
dim(count_erythro_cases_norm)

colnames(count_erythro_cases_norm)[1]="Transcript_ID" 
#15437

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#### Write and Read

write.table(count_erythro_cases_norm, file="/Users/saeedehazary/Documents/LWRNASEQ/count_erythro_cases_norm.txt", sep="\t", row.names=F)
count_erythro_cases_norm=read.table("/Users/saeedehazary/Documents/LWRNASEQ/count_erythro_cases_norm.txt", sep="\t", header=T)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

####################################

### creating design_erythro matrix 

####################################

erythro <- factor(c(rep("erythro", 3), rep("LWcase", 6)), levels = c("erythro", "LWcase"))

design_erythro <- model.matrix(~ erythro -1)
colnames(design_erythro) <- c ("erythro", "LWcase")
design_erythro


####################################

### erythro vs LWcase contrast

####################################

count_erythro_cases_norm_copy=count_erythro_cases_norm

row.names(count_erythro_cases_norm_copy)=count_erythro_cases_norm_copy[, 1]
count_erythro_cases_norm_copy=count_erythro_cases_norm_copy[, -1]

D_erythro<- DGEList(counts=count_erythro_cases_norm_copy)
D_erythro<- estimateGLMCommonDisp (D_erythro, design_erythro)
D_erythro<- estimateGLMTagwiseDisp(D_erythro, design_erythro)
fit_erythro <- glmFit(D_erythro, design_erythro)

lrt_case_erythro <- glmLRT(fit_erythro, contrast = c(1,-1))

exp_list_erythro_cases=data.frame(topTags(lrt_case_erythro , n=nrow(count_erythro_cases_norm))) 

exp_list_erythro_cases_2=cbind(rownames(exp_list_erythro_cases), exp_list_erythro_cases)
colnames(exp_list_erythro_cases_2)[1]="Ensembl.Transcript.ID"
row.names(exp_list_erythro_cases_2) <- NULL


### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_list_erythro_cases_annotate=merge(mart, exp_list_erythro_cases_2, by="Ensembl.Transcript.ID", all.exp_list_erythro_cases_2=all)

### sort

exp_list_erythro_cases_annotate_sort <- exp_list_erythro_cases_annotate[order(exp_list_erythro_cases_annotate$FDR),]
#15911

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#### Write and Read

write.table(exp_list_erythro_cases_annotate_sort, file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_erythro_cases_annotate_sort.txt", sep="\t", row.names=F)
exp_list_erythro_cases_annotate_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_erythro_cases_annotate_sort.txt", sep="\t", header=T)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

####################################

### data preparation for enrichment

####################################

FDR_less5_rbc=data.frame(subset(exp_list_erythro_cases_annotate_sort,exp_list_erythro_cases_annotate_sort$FDR <=0.05 & exp_list_erythro_cases_annotate_sort$logFC >0)) #  are expressed sig more in RBC 3012
FDR_less5_heart=data.frame(subset(exp_list_erythro_cases_annotate_sort,exp_list_erythro_cases_annotate_sort$FDR <=0.05 & exp_list_erythro_cases_annotate_sort$logFC<1)) #  are expressed sig more in heart
colnames(Hgb_mgi_symbol_list)[1]="Ensembl.Transcript.ID"


write.table(exp_list_erythro_cases_annotate_sort[,1], "/Users/saeedehazary/Documents/LWRNASEQ/all.for.118.txt", quote=F, col=F, row=F)

exp_list_annotate_sort_855= data.frame(subset(exp_list_annotate_sort, exp_list_annotate_sort$FDR<=0.05)) # new design used


rbc.in855 <- data.frame(intersect(FDR_less5_rbc$Ensembl.Transcript.ID, exp_list_annotate_sort_855$Ensembl.Transcript.ID)) #118




write.table(rbc.in855, "/Users/saeedehazary/Documents/LWRNASEQ/rbc.in855.txt", row.name=F, sep="\t")



colnames(rbc.in855)[1]="Ensembl.Transcript.ID"
rbc.in855.annot=merge(mart, rbc.in855, by="Ensembl.Transcript.ID", all.y=T)


write.table(rbc.in855.annot, "/Users/saeedehazary/Documents/LWRNASEQ/rbc.in855.annot.txt", row.name=F, sep="\t")


###################################################################################################################################

### creating a clean data set   

###################################################################################################################################


count_LW_clean=data.frame(subset(count_LW, !(count_LW$Transcript_name %in% rbc.in855$Ensembl.Transcript.ID)))

#82816

rownames(count_LW_clean) <-count_LW_clean[,1]

rownames(count_LW_clean) <- NULL

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

### Write and READ

write.table(count_LW_clean, file="/Users/saeedehazary/Documents/LWRNASEQ/count_LW_clean.txt", sep="\t", row.names=F)
count_LW_clean  =read.table(      "/Users/saeedehazary/Documents/LWRNASEQ/count_LW_clean.txt", header=T, sep="\t")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

###################################################################################################################################

### Normalizing count_LW_clean data set

###################################################################################################################################

library(limma)
library(edgeR)
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


count_norm_LW_clean=cbind(rownames(count_norm_LW_clean),count_norm_LW_clean)
rownames(count_norm_LW_clean)=NULL
dim(count_norm_LW_clean)
 #14251

colnames(count_norm_LW_clean)[1]="Transcript_ID" 

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

### Write and READ

write.table(count_norm_LW_clean, file="/Users/saeedehazary/Documents/LWRNASEQ/count_norm_LW_clean.txt", sep="\t", row.names=F)
count_norm_LW_clean=read.table("/Users/saeedehazary/Documents/LWRNASEQ/count_norm_LW_clean.txt", header=T, sep="\t")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

###################################################################################################################################

###  Making contrast in clean data set

###################################################################################################################################

count_norm_LW_clean_b=count_norm_LW_clean 
row.names(count_norm_LW_clean_b)=count_norm_LW_clean_b[, 1]
count_norm_LW_clean_b=count_norm_LW_clean_b[, -1]

D_clean <- DGEList(counts=count_norm_LW_clean_b)
D_clean <- estimateGLMCommonDisp(D_clean, design)  # design is the new version 
D_clean <- estimateGLMTagwiseDisp(D_clean, design)
fit_clean <- glmFit(D_clean, design)

###################

### c - (A+B) clean

###################

lrt_C_AB <- glmLRT(fit_clean, contrast = c(-0.5, -0.5, 1, 0, 0)) # Error in glmfit$coefficients %*% contrast : non-conformable arguments

exp_list_clean_C_AB=data.frame(topTags(lrt_C_AB, n=nrow(count_norm_LW_clean))) 
exp_list_clean_C_AB_2=cbind(rownames(exp_list_clean_C_AB), exp_list_clean_C_AB)
colnames(exp_list_clean_C_AB_2)[1]="Ensembl.Transcript.ID"
row.names(exp_list_clean_C_AB_2) <- NULL

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

### Write and READ

write.table(exp_list_clean_C_AB_2, "/Users/saeedehazary/Documents/LWRNASEQ/exp_list_clean_C_AB_new.txt", sep="\t", row.names=F)
read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_clean_C_AB_new.txt", sep="\t", header=T)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_list_annotate_clean_C_AB=merge(mart, exp_list_clean_C_AB_2, by="Ensembl.Transcript.ID", all.exp_list_clean_C_AB_2=all)
names( count_norm_LW)
count_norm_LW_cc=count_norm_LW
colnames(count_norm_LW_cc)[1]="Ensembl.Transcript.ID"
count_norm_LW_annotated=merge(mart, count_norm_LW_cc, by="Ensembl.Transcript.ID", all.count_norm_LW=all)

write.table(count_norm_LW_annotated,file="/Users/saeedehazary/Documents/LWRNASEQ/count_norm_LW_annotated.txt", sep="\t", row.names=F)

### sort

exp_list_annotate_clean_C_AB_sort <- exp_list_annotate_clean_C_AB[order(exp_list_annotate_clean_C_AB$FDR),] #15911

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

### Write and READ



write.table(exp_list_annotate_clean_C_AB_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_C_AB_sort.txt", sep="\t", row.names=F)

exp_list_annotate_clean_C_AB_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_C_AB_sort.txt", sep="\t", header=T)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

################################

### preparing file for enrichment

#################################

names(exp_list_annotate_clean_C_AB_sort)

s_exp_list_annotate_clean_C_AB_sort=data.frame(subset(exp_list_annotate_clean_C_AB_sort, exp_list_annotate_clean_C_AB_sort$FDR<= 0.05)) #first step design matrix on raw data was new one and created 818 sig and then we excluded 118 rbc  from main data frame and run again using new matrix

write.table(s_exp_list_annotate_clean_C_AB_sort[,1], "exp_annotate_clean_C_AB_sort_enrich.txt", quote=F, row=F, col=F)



####################
### pos vs neg

##################

s_exp_annot_clean_C_AB_sort_p=data.frame(subset(s_exp_list_annotate_clean_C_AB_sort, s_exp_list_annotate_clean_C_AB_sort$logFC>0))

s_exp_annot_clean_C_AB_sort_n=data.frame(subset(s_exp_list_annotate_clean_C_AB_sort, s_exp_list_annotate_clean_C_AB_sort$logFC<0))



write.table(s_exp_annot_clean_C_AB_sort_p, "s_exp_annot_clean_C_AB_sort_p.txt", row.names=F, sep="\t") 
write.table(s_exp_annot_clean_C_AB_sort_n, "s_exp_annot_clean_C_AB_sort_n.txt", row.names=F, sep="\t")




s_exp_annot_clean_C_AB_sort_p=read.table("/Users/saeedehazary/Documents/LWRNASEQ/s_exp_annot_clean_C_AB_sort_p.txt", sep="\t", header=T)
s_exp_annot_clean_C_AB_sort_n=read.table("/Users/saeedehazary/Documents/LWRNASEQ/s_exp_annot_clean_C_AB_sort_n.txt", sep="\t", header=T)

write.table(s_exp_annot_clean_C_AB_sort_p[,1], "/Users/saeedehazary/Documents/LWRNASEQ/enrich2.p.txt", quote=F, row=F, col=F)
write.table(s_exp_annot_clean_C_AB_sort_n[,1], "/Users/saeedehazary/Documents/LWRNASEQ/enrich2.n.txt", quote=F, row=F, col=F)
write.table(exp_list_annotate_clean_C_AB_sort[,1], "/Users/saeedehazary/Documents/LWRNASEQ/enrich.all.txt", quote=F, row=F, col=F)

exp_list_annotate_clean_C_AB_sort[,1]

write.table(exp_list_annotate_sort_855, "/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_sort_855.txt", row.names=F, sep="\t")





################
### A vs B clean

#################

lrt_A_B <- glmLRT(fit_clean, contrast = c(1, -1, 0, 0, 0))

exp_list_clean_A_B=data.frame(topTags(lrt_A_B, n=nrow(count_norm_LW_clean))) #14251 5
exp_list_clean_A_B_2=cbind(rownames(exp_list_clean_A_B), exp_list_clean_A_B)
colnames(exp_list_clean_A_B_2)[1]="Ensembl.Transcript.ID"
row.names(exp_list_clean_A_B_2) <- NULL



### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_list_annotate_clean_A_B=merge(mart, exp_list_clean_A_B_2, by="Ensembl.Transcript.ID", all.exp_list_clean_A_B_2=all)

### sort

exp_list_annotate_clean_A_B_sort <- exp_list_annotate_clean_A_B[order(exp_list_annotate_clean_A_B$FDR),] #14251


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
### Write and Read

write.table(exp_list_annotate_clean_A_B_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_A_B_sort.txt", sep="\t", row.names=F)
exp_list_annotate_clean_A_B_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_A_B_sort.txt", sep="\t", header=T)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

##################

### B vas C clean 

##################

lrt_B_C <- glmLRT(fit_clean, contrast = c(0, 1, -1, 0, 0))

exp_list_clean_B_C=data.frame(topTags(lrt_B_C, n=nrow(count_norm_LW_clean))) #13948 5
exp_list_clean_B_C_2=cbind(rownames(exp_list_clean_B_C), exp_list_clean_B_C)
colnames(exp_list_clean_B_C_2)[1]="Ensembl.Transcript.ID"
row.names(exp_list_clean_B_C_2) <- NULL



### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_list_annotate_clean_B_C=merge(mart, exp_list_clean_B_C_2, by="Ensembl.Transcript.ID", all.exp_list_clean_B_C_2=all)

### sort

exp_list_annotate_clean_B_C_sort <- exp_list_annotate_clean_B_C[order(exp_list_annotate_clean_B_C$FDR),]

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
### Write and Read

write.table(exp_list_annotate_clean_B_C_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_B_C_sort.txt", sep="\t", row.names=F)
exp_list_annotate_clean_B_C_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_B_C_sort.txt", sep="\t", header=T)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

##################

### A vs C clean

##################

lrt_A_C <- glmLRT(fit_clean, contrast = c(1, 0, -1,0,0))

exp_list_clean_A_C=data.frame(topTags(lrt_A_C, n=nrow(count_norm_LW_clean))) #13948 5
exp_list_clean_A_C_2=cbind(rownames(exp_list_clean_A_C), exp_list_clean_A_C)
colnames(exp_list_clean_A_C_2)[1]="Ensembl.Transcript.ID"
row.names(exp_list_clean_A_C_2) <- NULL



### Adding gene symbols, we get the mart file from Ensembl

mart=read.table("/Users/saeedehazary/Documents/LWRNASEQ/mart_export-2.txt", sep="\t", header=T)

exp_list_annotate_clean_A_C=merge(mart, exp_list_clean_A_C_2, by="Ensembl.Transcript.ID", all.exp_list_clean_A_C_2=all)

### sort

exp_list_annotate_clean_A_C_sort <- exp_list_annotate_clean_A_C[order(exp_list_annotate_clean_A_C$FDR),]

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
### Write and Read

write.table(exp_list_annotate_clean_A_C_sort,file="/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_A_C_sort.txt", sep="\t", row.names=F)
exp_list_annotate_clean_A_C_sort=read.table("/Users/saeedehazary/Documents/LWRNASEQ/exp_list_annotate_clean_A_C_sort.txt", sep="\t", header=T)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




###################################################################################################################################

### Heatmap

##################################################################################################################################


#########################################################

### ABC list

########################################################

exp_annot_clean_ABC=rbind(data.frame(subset(exp_list_annotate_clean_A_B, exp_list_annotate_clean_A_B$FDR<=0.05)),
                          data.frame(subset(exp_list_annotate_clean_A_C, exp_list_annotate_clean_A_C$FDR<=0.05)) ,
                          data.frame(subset(exp_list_annotate_clean_B_C, exp_list_annotate_clean_B_C$FDR<=0.05)))  #2350

goodlist <- data.frame(unique(exp_annot_clean_ABC$Ensembl.Transcript.ID)) #1642
#goodlist <- unique(exp_annot_clean_ABC$Associated.Gene.Name)

########################################################

colnames(goodlist)[1]="Ensembl.Transcript.ID"

heatmap.file=data.frame(subset(count_norm_LW_clean, count_norm_LW_clean$Transcript_ID %in% goodlist$Ensembl.Transcript.ID))


rownames(heatmap.file)=heatmap.file[,1]
heatmap.file=heatmap.file[, -1] 

install.packages("gplots")
library(gplots)
library(seriation)

f<-colorRampPalette(c("blue", "white", "firebrick"))
col.all <- f(256)

pdf("heatmap_postclean_postlitter.pdf")

finalmat <- scale(t(cpm(heatmap.file, prior.count=2, log=TRUE)))
#finalmat <- scale(t(log(count_norm_LW_clean_4, base=2)))
finalmat[finalmat > 1] = 1
finalmat[finalmat < -1] = -1
count_heatmap<- hmap(as.matrix(finalmat),   cexRow = 0.8, cexCol = 0.005, margins = c(10,5), options = list (prop=TRUE, col = col.all, main="Post.cleaning_Post.litter.adjustment", cex.main=0.5))

dev.off()







