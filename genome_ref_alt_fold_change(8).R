#############################################################################################################
# SA
### Regularized linear model 
### predictors : snps in  100 kb
### dependent variable:  ASE fold changefrom DegSeq
### Population: 66 CEU individuals

##############################################################################################################

####################################################################

# Files you need to read in, in this code

# Derhap_CEU_ID.txt
# DEGSeq outputs for all the individuals (we get this from eXpress)
# All the hapmap phase files from chr1 to chr22
# refFlat_hg18.txt
# ID_66.txt

####################################################################

setwd ("/Users/saeedehazary/Documents/RNAseq/src")

##########################Reading allele counts for all transcripts in 66 CEU population#########################

# To avoid running this for loop you can read the file 
count=read.table("/Users/saeedehazary/Documents/R_files/count.txt", sep="\t", header=T)

#################################################################################################################

#### Creating list of CEU IDs

list2=read.table("/Users/saeedehazary/Documents/trio/Derhap_CEU_ID.txt", sep="\t", header=F)
list2
colnames(list2)[1]="ID" 

#### Reading the Degseq out puts for each individuals in CEU population

for(j in list2$ID ) { 
  if (j=="NA06985"){
    a= sprintf("Der_%s_CEU/Degseq_from_bowtie/output_score.txt", j) 
    b=read.table(a, header=T, sep="\t" )
    colnames(b) [1]="Transcript_name"
    colnames(b) [2]=sprintf("%s_allele1_counts",j)
    colnames(b) [3]=sprintf("%s_allele0_counts",j)
    #colnames(b) [5]=sprintf("%s_log2.Fold_change_normalized",j)
    count=b[ ,c(1:3)]
    
  } else{
    
    a= sprintf("Der_%s_CEU/Degseq_from_bowtie/output_score.txt", j) 
    b=read.table(a, header=T, sep="\t" )
    colnames(b) [1]="Transcript_name"
    colnames(b) [2]=sprintf("%s_allele1_counts",j)
    colnames(b) [3]=sprintf("%s_allele0_counts",j)
    #colnames(b) [5]=sprintf("%s_log2.Fold_change_normalized",j)
    b=b[ ,c(1:3)]
    
    count<- merge(count,b,  by="Transcript_name", all=TRUE)
  }
}

dim(count)
#15935 133

colnames(count)[1]="Transcript_ID"

#write.table(count, file="/Users/saeedehazary/Documents/R_files/count.txt", sep="\t", row.names=F)

########################## Reading fold change for all transcripts in 66 CEU population ############################

# To avoid running this for loop you can read the file 
count_fold_change=read.table("/Users/saeedehazary/Documents/R_files/count_fold_change.txt", sep="\t", header=T)

#####################################################################################################################

#### Creating list of CEU IDs

list2=read.table("/Users/saeedehazary/Documents/trio/Derhap_CEU_ID.txt", sep="\t", header=F)
list2
colnames(list2)[1]="ID"

#### Reading the Degseq out puts for each individuals in CEU population

for(j in list2$ID ) { 
  if (j=="NA06985"){
    a= sprintf("Der_%s_CEU/Degseq_from_bowtie/output_score.txt", j) 
    b=read.table(a, header=T, sep="\t" )
    colnames(b) [1]="Transcript_name"
    colnames(b) [4]=sprintf("%s_log2.Fold_change",j)
    count_fold_change=b[ ,c(1,4)]
    
  } else{
    
    a= sprintf("Der_%s_CEU/Degseq_from_bowtie/output_score.txt", j) 
    b=read.table(a, header=T, sep="\t" )
    colnames(b) [1]="Transcript_name"
    colnames(b) [4]=sprintf("%s_log2.Fold_change",j)
    b=b[ ,c(1,4)]
    
    count_fold_change<- merge(count_fold_change,b,  by="Transcript_name", all=TRUE)
  }
}

colnames(count_fold_change)[1]="Transcript_ID"

#### Writting the data

#write.table(count_fold_change, file="/Users/saeedehazary/Documents/R_files/count_fold_change.txt", sep="\t", row.names=F)


###########################################################################################################

# Normalizing the count data

###########################################################################################################


#################### Reading counts for all transcripts (non ASE) in 66 CEU population####################

#To avoid running this for loop you can read the file 
count_standardize=read.table("/Users/saeedehazary/Documents/R_files/count_standardize.txt", sep="\t", header=T)

##########################################################################################################


list2=read.table("/Users/saeedehazary/Documents/trio/Derhap_CEU_ID.txt", sep="\t", header=F)
list2
colnames(list2)[1]="ID"

#### Reading the eXpress out puts for each individuals in CEU population

for(j in list2$ID ) { 
  if (j=="NA06985"){
    a= sprintf("/Users/saeedehazary/Documents/RNASeq/src/Der_all_hap/Der_%s_CEU/results.xprs", j) 
    b=read.table(a, header=T, sep="\t" )
    colnames(b) [2]="Transcript_name"
    colnames(b) [8]=sprintf("%s_eff_counts",j)
    count_standardize=b[ ,c(2,8)]
    
  } else{
    
    a= sprintf("/Users/saeedehazary/Documents/RNASeq/src/Der_all_hap/Der_%s_CEU/results.xprs", j) 
    b=read.table(a, header=T, sep="\t" )
    colnames(b) [2]="Transcript_name"
    colnames(b) [8]=sprintf("%s_eff_counts",j)
    b=b[ ,c(2,8)]
    
    count_standardize<- merge(count_standardize,b,by="Transcript_name", all=TRUE)
    print(j)
    
  }
}

dim(count_standardize)
#[1] 46568    67

#write.table(count_standardize, file="/Users/saeedehazary/Documents/R_files/count_standardize.txt", sep="\t", row.names=F)

########################################################################################################

# Creating norm factor and scale factor using edge R and applying those on the count data set

########################################################################################################

library(edgeR)

counts_all <- count_standardize
counts <- counts_all[,2:ncol(counts_all)]
rownames(counts) <- counts.all[,1]
keep <- rowSums(cpm(counts)>2)>=60
toc <- counts[keep,]
table(keep)
norm_factors <- calcNormFactors(as.matrix(toc), method = "TMM")

norm_factors=data.frame(norm_factors)
norm_factors_dup=norm_factors[rep(1:nrow(norm_factors), each=2), ]
print(norm_factors_dup)
#1:132


count_132=count[,2:ncol(count)]
rownames(count_132)=count[,1]
colsum=colSums(count_132, na.rm=TRUE)
scalefactor <- mean(colsum)/colsum
print(scalefactor)
#1:132

count_norm <- data.frame(sweep(as.matrix(count_132), MARGIN = 2,norm_factors_dup*scalefactor, '*'))
#15935 132

count_norm=cbind(rownames(count_norm),count_norm , row.names=NULL)
dim(count_norm)
#15935 133

colnames(count_norm)[1]="Transcript_ID" 

############################################################################################################################

# Genome file preparation

############################################################################################################################

#### readinfg the refFlat file(pulling out TSS)

refFlat=read.table("refFlat_hg18.txt",header=F, sep="\t")
colnames(refFlat) [2]="Transcript_ID"
colnames(refFlat) [3]="chr"

#### creating new variable using TSS+_100000

nmax = length(refFlat$V4)
TSSmax<-rep(NA,nmax)
TSSmin<-rep(NA,nmax)

for (i in 1:nmax) {
  
  if (refFlat$V4[i] =="+")  {
    TSSmax[i]=refFlat$V5[i] + 100000
    TSSmin[i]=refFlat$V5[i] - 100000
  } else{
    TSSmax[i]=refFlat$V6[i] + 100000
    TSSmin[i]=refFlat$V6[i] - 100000
  }
}

#### Run one of the below , not all#################################################################################

refFlat.new_100 <- data.frame(cbind(refFlat[,1:6],TSSmax,TSSmin))
refFlat.new_250 <- data.frame(cbind(refFlat[,1:6],TSSmax,TSSmin))
refFlat.new_50 <- data.frame(cbind(refFlat[,1:6],TSSmax,TSSmin))

########## Run one of the below, not all,  Excluding some transcripts from refFlat file##############################

#excluding chrX transcripts and duplicte transcripts

refFlat.new=subset(refFlat.new_250, refFlat.new_250$chr!="chrX")
#46181

refFlat.new=subset(refFlat.new_100, refFlat.new_100$chr!="chrX")
#46181

refFlat.new=subset(refFlat.new_50, refFlat.new_50$chr!="chrX")

################################################# more exclusion on refFlat file ####################################

refFlat.new=subset(refFlat.new, refFlat.new$chr!="chrY")
#45823

refFlat.new_dup=subset(refFlat.new, duplicated(refFlat.new$Transcript_ID))

refFlat.new_nodup=subset(refFlat.new, !(refFlat.new$Transcript_ID %in% refFlat.new_dup$Transcript_ID))
#43115 8

#### checking the data

c=subset(refFlat.new_nodup, duplicated(refFlat.new_nodup$Transcript_ID))
grep("NM_001127389", refFlat.new_nodup$Transcript_ID)
grep("NM_173844", refFlat.new_nodup$Transcript_ID)
grep("NM_173565", refFlat.new_nodup$Transcript_ID)
#some chtr column have random and it creates an error in the later for loop and we need to exclude those
random=grep("random", refFlat.new_nodup$chr, value=TRUE)
random

#### we need to change the chr variable from integer to character

typeof(refFlat.new_nodup$chr)
#integer 
refFlat.new_nodup<- data.frame(lapply(refFlat.new_nodup, as.character), stringsAsFactors=FALSE)
typeof(refFlat.new_nodup$chr)
#character

#### excluding the transcript with string "random" in the chr column

random_file=refFlat.new_nodup[grepl("random",refFlat.new_nodup$chr), ]
refFlat.new_nodup_b=subset(refFlat.new_nodup, !(refFlat.new_nodup$Transcript_ID %in% random_file$Transcript_ID))
dim(refFlat.new_nodup_b)
#43070 8

################################ subseting the refFlat ###################################################################

#### we need to include only 15700 transcripts since these 15700 are in both genome and count file

#### before merge we need to remove "_" from TRanscript_Id in refFlat file

k=data.frame(sub ('_', '', refFlat.new_nodup_b$Transcript_ID))
refFlat.new_nodup_b=cbind(k, refFlat.new_nodup_b)
colnames(refFlat.new_nodup_b)[3]="T_ID"
colnames(refFlat.new_nodup_b)[1]="Transcript_ID"
#43070 9

#### only those TRanscript_IDs in refFlat file that we have data for that in count

#refFlat.new_nodup_c=merge(refFlat.new_nodup_b, count_fold_change, by="Transcript_ID", all.count_fold_change=all)
refFlat.new_nodup_c=merge(refFlat.new_nodup_b, count, by="Transcript_ID", all.count=all)
dim(refFlat.new_nodup_c)
#15700 141

refFlat.new_nodup_c=refFlat.new_nodup_c[order(refFlat.new_nodup_c$Transcript_ID),]

#### checking the data 

y=count_norm$Transcript_ID[(count_norm$Transcript_ID %in% refFlat.new_nodup_b$Transcript_ID)]
#15700 (These are are in both file)
t=count_norm$Transcript_ID[!(count_norm$Transcript_ID %in% refFlat.new_nodup_b$Transcript_ID)]
#235 (these are existed in st_count_d but not in refFlatfile)
t

#### Checked___I checked these 235 in the refflat file they are the transcripts with non clear entry and duplicated

#### sort the refFlat.new_nodup_c file

refFlat.new_nodup_c=refFlat.new_nodup_c[order(refFlat.new_nodup_c$Transcript_ID),]
names(refFlat.new_nodup_c)#Transcript_ID is the first column
#15700 141 

##### change Transcrippt_ID to character

p=as.character(refFlat.new_nodup_c$Transcript_ID)
typeof(p)
refFlat.new_nodup_c["Transcript_ID_char"]=p

dim(refFlat.new_nodup_c)
refFlat.new_nodup_c=refFlat.new_nodup_c[, c(142,2:141)]

colnames(refFlat.new_nodup_c) [1]="Transcript_ID"
typeof(refFlat.new_nodup_c$Transcript_ID)

#################################### Preparing coun_norm_data #########################################################

#### we need as.character to be able to sort the data 

typeof(count_norm$Transcript_ID)
l=as.character(count_norm$Transcript_ID)
count_norm['Transcript_ID']=l

#### sorting the data by ID

count_norm_ordered=count[order(count_norm$Transcript_ID),]


#### We want to exclude those 235 transcripts that are excluded from refflat file because of error creation and we donot have genome file for them 
#### and make count file to have 15700 Transcripts 

count_norm_b=count_norm_ordered[(count_norm_ordered$Transcript_ID %in% refFlat.new_nodup_c$Transcript_ID), ]
dim(count_norm_b)
#15700 133 (These are are in both file)

#### removing row.names

count_norm_b=data.frame(count_norm_b, row.names=NULL)
#15700 133

########################################## Preparing count_fc data ##################################################

#### we need to change to as.character to be able to sort the data 

typeof(count_fold_change$Transcript_ID)
l=as.character(count_fold_change$Transcript_ID)
count_fold_change['Transcript_ID']=l

#### sorting the data by ID

count_fold_change_ordered=count_fold_change[order(count_fold_change$Transcript_ID),]

#### We want to exclude those 235 transcripts that are excluded from refflat file because of error creation and we donot have genome similarity matrix for them 
#### and make count file to have  15700 Transcripts 

count_fc=count_fold_change_ordered[(count_fold_change_ordered$Transcript_ID %in% refFlat.new_nodup_c$Transcript_ID), ]
###count_fold_change=data.frame(count_fold_change, row.names=NULL)

dim (count_fc)
#15700 67  (These are are in both file, genome and count_fc)

count_fc=data.frame(count_fc, row.names=NULL)

#######################################################################################################################################

#### Running forloop to create a genome data frame in which we compare alleles intrapersonal and create codes as

# 0 : both alleles have refrence snp
# 1 : allele hap0 has refrence snp and allele hap1 has alternative snp
# 2 : allele hap1 has refrence snp and allele hap0 has refrence snp
# 3 : both alleles have alternative snps

################################################# forloop to crate a genome file per transcript ######################################################

#### packages neede for forloop

#install.packages("matrixStats")
library(matrixStats)
library(gdata)
#install.packages("dummies")
library(dummies)

#### counting number of enries in count_fc data set


#### reading the refrence , alternative alleles file

ref_alt_allchr_consensus_CEU=read.table("/Users/saeedehazary/Documents/RNASeq/src/ref_alt/ref_alt_allchr_consensus_CEU.txt", sep="\t", header=T)

#### forloop to crate a genpme file per transcript

list_386_678=list_678[386:678]

for (t in list_386_678 ) {
  
print(t)

refFlat_a= refFlat.new_nodup_c[t, ]
x=refFlat_a$chr
u=read.table(sprintf("~/Dropbox/DeoLabShared/RNAseq/src/python_scripts_hap_files/hapmap3_r2_b36_fwd.consensus.qc.poly.%s_ceu_clean.phased", x), sep="\t", header=TRUE)
k= subset(u, refFlat_a$TSSmin<=u$phys_position & u$phys_position<=refFlat_a$TSSmax )
  
  if (nrow(k)==0) {
    next
  }
  
  else {
    
    k_char= data.frame(lapply(k, as.character), stringsAsFactors=FALSE)
    
    for (g in 3:ncol(k_char)){
      
      if (g==3)
        {
        k_char_split= data.frame(do.call(rbind, strsplit(as.vector(k_char[,g]), split = " ")))
        colnames(k_char_split) <- c(colnames(k_char)[g],colnames(k_char)[g])
      }
      
      else
      {
        b= data.frame(do.call(rbind, strsplit(as.vector(k_char[,g]), split = " ")))
        colnames(b) <- c(colnames(k_char)[g],colnames(k_char)[g])
        k_char_split= data.frame(k_char_split,b)  
      }
      
    }
      
 ID_66=read.table("ID_66.txt" , header=F, sep="\t")
 myvars <- names(k_char_split) %in% ID_66$V1 
 k_char_split=k_char_split[myvars]
#132 columns
      
#### adding the name of snps to first column
    
genome_split=cbind(k$rsID, k_char_split)
colnames(genome_split)[1]="rsID"
    
#### adding the ref and alt column
    
ref_alt_genome=merge (ref_alt_allchr_consensus_CEU, genome_split, by="rsID", all.genome_split=all)
      
     
#### We need names (which includes 66 names of CEU ID individuals) to create the column name for our new columns later in line 462
     
names= names(ref_alt_genome)[seq(4, 134, 2)] 

ref_alt_genome_b=ref_alt_genome
ref_alt_genome_b<- data.frame(lapply(ref_alt_genome_b, as.character), stringsAsFactors=FALSE)
   
 
#### for loop for creating 0,1,2,3 codes for intrapersonal alleles comparision

ptm <- proc.time()
    for (j in seq(4, 134, 2)) {
      
      for (i in 1: nrow(ref_alt_genome_b)) {
        
        if      (ref_alt_genome_b[i,j]==ref_alt_genome_b[i,2] && ref_alt_genome_b[i,j+1]==ref_alt_genome_b[i,2]) {ref_alt_genome_b$g_b[i]=0}  #ref is in column 2
        else if (ref_alt_genome_b[i,j]==ref_alt_genome_b[i,2] && ref_alt_genome_b[i,j+1]==ref_alt_genome_b[i,3]) {ref_alt_genome_b$g_b[i]=1}
        else if (ref_alt_genome_b[i,j]==ref_alt_genome_b[i,3] && ref_alt_genome_b[i,j+1]==ref_alt_genome_b[i,2]) {ref_alt_genome_b$g_b[i]=2}
        else if (ref_alt_genome_b[i,j]==ref_alt_genome_b[i,3] && ref_alt_genome_b[i,j+1]==ref_alt_genome_b[i,3]) {ref_alt_genome_b$g_b[i]=3}
        else {ref_alt_genome_b$g_b[i]="NA"}
        
      }
      names(ref_alt_genome_b)[names(ref_alt_genome_b)=='g_b']<-sprintf("compare%d", j)
      print(j)
    }
print(proc.time() - ptm)
    
colnames(ref_alt_genome_b)[136:201]= paste("code", names, sep="_")
   
#### Substing data set to keep snp ID and new columns 

ref_alt_genome_c=ref_alt_genome_b[,c(1, 136:201)]
a3a=refFlat_a$Transcript_ID

#### writting the final file

write.table(ref_alt_genome_c, file=sprintf("/Users/saeedehazary/Documents/R_files/%s_ref_alt_genome", a3a), sep="\t", row.names=F)

  }    
}

###################################################forloop to crate glmnet output using genome and count_fc data #################

# these are the data that we need to use
### sub_count_fc
### ref_alt_genome_c

###################################################################################################################################

list_107_382=list_1_382[107:382]
list_107_382[32]
list_33=list_107_382[33:length(list_107_382)]
list_93=list_33[93:length(list_33)]
list_1_494=list_678[1:494]

library(glmnet)

for (yy in list_1_494){
  
  print (yy)
  
  if (yy==8320) { 
  next
  }
  else {
  
  #if (yy==156|yy==258| yy==441 | yy==442 | yy==443 | yy==934|yy==1309 |yy==1318|yy==1365 |yy==2388 |yy==2403|yy==2825|
       # yy==2900 |yy==3268 |yy==3380 |yy==3528| yy==5269 |yy==5377| yy==6161| yy==6162 | yy==6189){
  #next
#}

#else {
    
sub_count_fc=data.frame(subset(count_fc, count_fc$Transcript_ID==count_fc[yy,1]), row.names=NULL)
# 1 67

z3z=sub_count_fc$Transcript_ID

if (!file.exists(sprintf("/Users/saeedehazary/Documents/R_files/%s_ref_alt_genome", z3z ))) {
  next
}

else{
  
ref_alt_genome_c=read.table(sprintf("/Users/saeedehazary/Documents/R_files/%s_ref_alt_genome", z3z ), sep="\t", header=T)

sub_count_fc=sub_count_fc[ ,c(2:67)]
# 1 66
mean=rowMeans (sub_count_fc , na.rm=TRUE)
sd=rowSds(sub_count_fc, na.rm=TRUE)
x=mean-(2*sd)
y=mean+(2*sd)
#excluding outliers
#values below mean_2sd changes to NA
sub_count_fc[sub_count_fc<x]<-NA
#dim(sub_count_fc)
#1 66
#ci=(ni-mean)/sd
#sub_count_fold_change<- (sub_count_fold_change-mean)/sd

mean2=rowMeans(abs(sub_count_fc),na.rm=TRUE )

#### Transposing genome data set to use snps as predictor in linear model

t_ref_alt_genome_c=as.data.frame(t(ref_alt_genome_c))

#### Keeping the rownames as new column

t_ref_alt_genome_c$rsID=factor(row.names(t_ref_alt_genome_c))

#### as.character

t_ref_alt_genome_c= data.frame(lapply(t_ref_alt_genome_c, as.character), stringsAsFactors=FALSE)

#### duplicating first row to colnames

colnames(t_ref_alt_genome_c)=t_ref_alt_genome_c[1,]

#### deleting first row
t_ref_alt_genome_c=t_ref_alt_genome_c[-1,]

#### removing individual IDs to first column
col_1=ncol(t_ref_alt_genome_c)-1
t_ref_alt_genome_c=t_ref_alt_genome_c[, c(ncol(t_ref_alt_genome_c),1:col_1)]

#### changing rsID to ID
colnames(t_ref_alt_genome_c)[1]="ID"

#### as.character
cc=as.character(list2$ID)
list2["ID"]=cc
t_ref_alt_genome_c$fc="NA"

#### adding the fc counts to transposed genome ,  predictors (snps) and dependent variable (fc) in one data

for (p in list2$ID) {
  
  a=which(t_ref_alt_genome_c$ID==sprintf("code_%s_A.%s_B", p, p) )
  #print(a)
  b=which(colnames(sub_count_fc)==sprintf("%s_log2.Fold_change", p) )
  #print(b)
  t_ref_alt_genome_c$fc[a]=sub_count_fc[,b]
  
}

head(t_ref_alt_genome_c) #fc added to the end of data frame

#### saving the number of columns in data set before creating dummy variables

s=ncol(t_ref_alt_genome_c)

#### changing the order of column : col1 is ID , col2 is fold change and col3 to end is snps 

t_ref_alt_genome_d=t_ref_alt_genome_c[, c(1, ncol(t_ref_alt_genome_c), 2:(ncol(t_ref_alt_genome_c)-1))]

#######creating dummies using dummies package

t_ref_alt_genome_e=t_ref_alt_genome_d
t_ref_alt_genome_f=t_ref_alt_genome_e[, 3:ncol(t_ref_alt_genome_e)]

t_ref_alt_genome_f=dummy.data.frame(t_ref_alt_genome_f,  sep="_" , dummy.class="ALL", drop=FALSE)

t_ref_alt_genome_g=data.frame(cbind(as.numeric(t_ref_alt_genome_d$fc), t_ref_alt_genome_f))
colnames(t_ref_alt_genome_g) [1]= "fc"
t_ref_alt_genome_h<- t_ref_alt_genome_g[!is.na(t_ref_alt_genome_g$fc), ]
#as.integer


#### glmnet

if (yy==list_1_494[1]) {

thresh = ncol(t_ref_alt_genome_h)

fit <- glmnet(as.matrix(t_ref_alt_genome_h[,2:thresh]), as.numeric(t_ref_alt_genome_h[,1]), alpha = 1)
cvfit=cv.glmnet(as.matrix(t_ref_alt_genome_h[,2:thresh]), as.numeric(t_ref_alt_genome_h[,1]), alpha = 1)
#cvfit
mlambda = cvfit$lambda.min
Coefficients <- coef(fit, s = cvfit$lambda.min)
#Coefficients
Active.Index.min <- which(Coefficients != 0)
Active.Coefficients.min <- Coefficients[Active.Index.min]
#Active.Index.min # identifies the covariates that are active in the model and
n_coef_min=length(Active.Coefficients.min) # shows the coefficients of those covariates
cvm_mlambda=cvfit$cvm[which(cvfit$lambda==cvfit$lambda.min)]
cvsd_mlambda=cvfit$cvsd[which(cvfit$lambda==cvfit$lambda.min)]
cvm_mlambda_std=cvm_mlambda/mean2
T_ID=sprintf("%s", z3z )
n_snps=nrow(ref_alt_genome_c)

mylist<-list(T_ID,n_snps, mlambda, n_coef_min, cvm_mlambda, cvsd_mlambda, cvm_mlambda_std)
names(mylist)<- c("Transcript_ID", "number_of_snps", "min.lambda", "number_coef_min", "cvm_mlambda", "cvsd_mlambda", "cvm_mlambda_std" )

#### z3z is the name of transcript 

#assign(sprintf("%s_glmnet_results", z3z ), mylist)
#write.list(mylist, filename=sprintf("/Users/saeedehazary/Documents/R_files/%s_glmnet_results", z3z), append=FALSE, closefile = TRUE, outfile)

LS.df_f = as.data.frame(do.call(cbind, mylist))

}

else{
  
  thresh = ncol(t_ref_alt_genome_h)
  fit <- glmnet(as.matrix(t_ref_alt_genome_h[,2:thresh]), as.numeric(t_ref_alt_genome_h[,1]), alpha = 1)
  cvfit=cv.glmnet(as.matrix(t_ref_alt_genome_h[,2:thresh]), as.numeric(t_ref_alt_genome_h[,1]), alpha = 1)
  #cvfit
  mlambda = cvfit$lambda.min
  Coefficients <- coef(fit, s = cvfit$lambda.min)
  #Coefficients
  Active.Index.min <- which(Coefficients != 0)
  Active.Coefficients.min <- Coefficients[Active.Index.min]
  #Active.Index.min # identifies the covariates that are active in the model and
  n_coef_min=length(Active.Coefficients.min) # shows the coefficients of those covariates
  cvm_mlambda=cvfit$cvm[which(cvfit$lambda==cvfit$lambda.min)]
  cvsd_mlambda=cvfit$cvsd[which(cvfit$lambda==cvfit$lambda.min)]
  cvm_mlambda_std=cvm_mlambda/mean2
  T_ID=sprintf("%s", z3z )
  n_snps=nrow(ref_alt_genome_c)
  
  mylist<-list(T_ID, n_snps, mlambda, n_coef_min, cvm_mlambda, cvsd_mlambda, cvm_mlambda_std)
  names(mylist)<- c("Transcript_ID", "number_of_snps", "min.lambda", "number_coef_min", "cvm_mlambda", "cvsd_mlambda", "cvm_mlambda_std" )
  
  #assign(sprintf("%s_glmnet_results", z3z ), mylist)
  #write.list(mylist, filename=sprintf("/Users/saeedehazary/Documents/R_files/%s_glmnet_results", z3z), append=FALSE, closefile = TRUE, outfile)
  
  LS.df2 = as.data.frame(do.call(cbind, mylist))
  LS.df_f= rbind(LS.df_f, LS.df2)
}

}
}
}

write.table(LS.df_f, file="/Users/saeedehazary/Documents/R_files/LS.df_f.txt", sep="\t", row.names=F)




install.packages("covTest")
library(covTest)

qq=lars.en(as.matrix(t_ref_alt_genome_h[,2:thresh]), as.numeric(t_ref_alt_genome_h[,1]),lambda2=1)
qq=lars.glm(as.matrix(t_ref_alt_genome_h[,2:thresh]), as.numeric(t_ref_alt_genome_h[,1]),family="binomial")
covTest(qq, as.matrix(t_ref_alt_genome_h[,2:thresh]), as.numeric(t_ref_alt_genome_h[,1]))


fit#### z3z is the name of transcript 
















