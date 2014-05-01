setwd("/Users/saeedehazary/Documents/RNAseq/src")

####################################################################
#Files you need to read in to run this script
####################################################################
#Derhap_CEU_ID.txt
#DEGSeq outputs for all the individuals (we get this from eXpress)
#All the hapmap phase files from chr1 to chr22
#refFlat_hg18.txt
#ID_66.txt



########################################################################################################################################
#Preparing input file to create count diff matrix
########################################################################################################################################

##########################Reading allele counts for all transcripts in 66 CEU population############################

#Creating list of CEU IDs

list2=read.table("/Users/saeedehazary/Documents/trio/Derhap_CEU_ID.txt", sep="\t", header=F)
list2
colnames(list2)[1]="ID"


#Reading the Degseq out puts for each individuals in CEU population

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

#count_all<- merge(count,list, all.count=TRUE, by="Transcript_name")
#colnames(data_frame_all)[200:289]="snp"
#names(data_frame_all)

###########################################################################################################
                                       #Normalizing the count data
###########################################################################################################
#Rahul_SA

#Normalizing using edge R

#Get the final file from "Standardize.R" script  

# The final file is "count_norm" and use this in the fit forloop


#########################################
#Creating the dataframe from all CEU individuals and their transcripts with eff_counts in one file
##################################################################################################


list2=read.table("/Users/saeedehazary/Documents/trio/Derhap_CEU_ID.txt", sep="\t", header=F)
list2
colnames(list2)[1]="ID"

#Reading the eXpress out puts for each individuals in CEU population

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


####################
#Rahul code
####################
library(edgeR)
counts_all <- count_standardize
counts <- counts_all[,2:ncol(counts_all)]
rownames(counts) <- counts.all[,1]
keep <- rowSums(cpm(counts)>2)>=60
toc <- counts[keep,]
table(keep)
norm_factors <- calcNormFactors(as.matrix(toc), method = "TMM")


###Me
norm_factors=data.frame(norm_factors)
norm_factors_dup=norm_factors[rep(1:nrow(norm_factors), each=2), ]
print(norm_factors_dup)
#1:132

########count data is from genome_count_fit.R code which is allele specific counts matrix with 133 columns , first column is transcript ID and 132 columns are 66 CEU pop
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
################################################################################################################################
################################################################################################################################


########################################################################################################################################
#preparing input files to create genome similarity matrix
########################################################################################################################################

#Reading phased files, Using assign to create a dataset 

for(i in 1:22 ) { 
  a= sprintf("~/Dropbox/DeoLabShared/RNAseq/src/python_scripts_hap_files/hapmap3_r2_b36_fwd.consensus.qc.poly.chr%d_ceu_clean.phased", i)
  assign(sprintf("chr%d_CEU", i), read.table(a, sep="\t", header=TRUE))
}
######################################################################################################################

#readinfg the refFlat file(pulling out TSS)

refFlat=read.table("refFlat_hg18.txt",header=F, sep="\t")
colnames(refFlat) [2]="Transcript_ID"
colnames(refFlat) [3]="chr"


######################################################################################################################

#creating new variable using TSS+_250000

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
refFlat.new_100 <- data.frame(cbind(refFlat[,1:6],TSSmax,TSSmin))
refFlat.new_250 <- data.frame(cbind(refFlat[,1:6],TSSmax,TSSmin))
refFlat.new_50 <- data.frame(cbind(refFlat[,1:6],TSSmax,TSSmin))

########################################## Excluding some transcripts from refFlat file###############################

#excluding chrX transcripts and duplicte transcripts

refFlat.new=subset(refFlat.new_250, refFlat.new_250$chr!="chrX")
#46181

refFlat.new=subset(refFlat.new_100, refFlat.new_100$chr!="chrX")
#46181

refFlat.new=subset(refFlat.new_50, refFlat.new_100$chr!="chrX")

######################################################################################################################


refFlat.new=subset(refFlat.new, refFlat.new$chr!="chrY")
#45823
refFlat.new_dup=subset(refFlat.new, duplicated(refFlat.new$Transcript_ID))
refFlat.new_dup

refFlat.new_nodup=subset(refFlat.new, !(refFlat.new$Transcript_ID %in% refFlat.new_dup$Transcript_ID))
#43115 8

#checking the data
c=subset(refFlat.new_nodup, duplicated(refFlat.new_nodup$Transcript_ID))
grep("NM_001127389", refFlat.new_nodup$Transcript_ID)
grep("NM_173844", refFlat.new_nodup$Transcript_ID)
grep("NM_173565", refFlat.new_nodup$Transcript_ID)
#some chtr column have random and it creates an error in the later for loop and we need to exclude those
random=grep("random", refFlat.new_nodup$chr, value=TRUE)
random

#we need to change the chr variable from integer to character
typeof(refFlat.new_nodup$chr)
#integer 
refFlat.new_nodup<- data.frame(lapply(refFlat.new_nodup, as.character), stringsAsFactors=FALSE)
typeof(refFlat.new_nodup$chr)
#character

#excluding the transcript with string "random" in the chr column
random_file=refFlat.new_nodup[grepl("random",refFlat.new_nodup$chr), ]
refFlat.new_nodup_b=subset(refFlat.new_nodup, !(refFlat.new_nodup$Transcript_ID %in% random_file$Transcript_ID))
dim(refFlat.new_nodup_b)
#43070 8

################################subseting the refFlat###################################################################

#we need to just run for 15700 transcripts ince these 15700 are in both genome similarity and count file

#before merge we need to remove "_" from TRanscript_Id in refFlat file
k=data.frame(sub ('_', '', refFlat.new_nodup_b$Transcript_ID))
refFlat.new_nodup_b=cbind(k, refFlat.new_nodup_b)
colnames(refFlat.new_nodup_b)[3]="T_ID"
colnames(refFlat.new_nodup_b)[1]="Transcript_ID"

#only those TRanscript_IDs in refFlat file that we have data for that in count

refFlat.new_nodup_c=merge(refFlat.new_nodup_b, count, by="Transcript_ID", all.count=all)
#15700 141


#checking
y=count_norm$Transcript_ID[(count_norm$Transcript_ID %in% refFlat.new_nodup_b$Transcript_ID)]
#15700 (These are are in both file)
t=count_norm$Transcript_ID[!(count_norm$Transcript_ID %in% refFlat.new_nodup_b$Transcript_ID)]
#235 (these are existed in st_count_d but not in refFlatfile)
t
#Checked___I checked these 235 in the refflat file they are the transcripts with non clear entry and duplicated

#sort the refFlat.new_nodup_c file
refFlat.new_nodup_c=refFlat.new_nodup_c[order(refFlat.new_nodup_c$Transcript_ID),]
names(refFlat.new_nodup_c)#Transcript_ID is the first column
#15700 141 

#change Transcrippt_ID to character
p=as.character(refFlat.new_nodup_c$Transcript_ID)
typeof(p)
refFlat.new_nodup_c["Transcript_ID_char"]=p

dim(refFlat.new_nodup_c)
refFlat.new_nodup_c=refFlat.new_nodup_c[, c(142,2:141)]


colnames(refFlat.new_nodup_c) [1]="Transcript_ID"
typeof(refFlat.new_nodup_c$Transcript_ID)

#############################################################################################################


#sub_count_norm=merge(refFlat.new_nodup_b,count_norm, by="Transcript_ID", all.count_norm=all)

#we do not need column2:9
#sub_count_norm=sub_count_norm[, c(1,10:141)]
#15700 133



#checking to make sure sure merging went well
#grep("NM000017", st_count_d$Transcript_ID)
#st_count_d[13367, ]

#grep("NM000017", count$Transcript_name)
#count[13367, ]
#st_count[13367,]

#grep("NM000017", st_count_e$Transcript_ID)
#st_count_e[1,]


#checking
###grep('NM000017', count_norm$Transcript_ID)
###grep('NM000022', count_norm$Transcript_ID)

# we need as.character to be able to sort the data 
typeof(count_norm$Transcript_ID)
l=as.character(count_norm$Transcript_ID)
count_norm['Transcript_ID']=l

#sorting the data by ID
count_norm_ordered=count[order(count_norm$Transcript_ID),]


#We want to exclude those 235 transcripts that are excluded from refflat file because of error creation and we donot have genome similarity matrix for them 
#and make count file to have  15700 Transcripts 

count_norm_b=count_norm_ordered[(count_norm_ordered$Transcript_ID %in% refFlat.new_nodup_c$Transcript_ID), ]
dim(count_norm_b)
#15700 133 (These are are in both file)

# removing row.names
count_norm_b=data.frame(count_norm_b, row.names=NULL)
dim(count_norm_b)
#15700 133


########################################################################################################################################
#Running forloop to create genome similarity matrix and count diff matrix and fit linear model between two triangle matrix
########################################################################################################################################
#for (t in 1:nrow(refFlat.new_nodup_c))


install.packages("matrixStats")
library(matrixStats)
library(gdata)

#x1=runif (10, 1, 15700)
#x2=runif (10, 1, 15700)
#set.seed(10)
for (t in 1:1000) {
  print(t)
  #refFlat_a= subset(refFlat.new_nodup_b, refFlat.new_nodup_b$Transcript_ID==(x))
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
      
      
      ID_66=read.table("ID_66.txt" , header=F, sep="\t")
      myvars <- names(k_char_split) %in% ID_66$V1 
      k_char_split=k_char_split[myvars]
      #132 columns
    }
    
    # if (t==394) {
    #   next
    #}
    # else{
    
    genome=matrix(data=NA, nrow=132, ncol=132, byrow=TRUE)
    for (i in 1:132){
      
      for( j in 1:132) {
        
        
        a=sum(k_char_split[,i]==k_char_split[,j])/nrow(k_char_split)
        genome [i, j]=a
        
      }
    }
    
    genome_triangle=upperTriangle(genome, diag=FALSE)
    #s=assign(sprintf("genome_%s_250", refFlat.new_nodup_c[t,1]), upperTriangle(genome, diag=FALSE))
    #save(s, file=sprintf("/Users/saeedehazary/Documents/R_files/permanent_files/genome_%s_250", refFlat.new_nodup_c[t,1]))
    
    
    
    #subseting a transcript
    sub_count_norm_b=data.frame(subset(count_norm_b, count_norm_b$Transcript_ID==count_norm_b[t,1]), row.names=NULL)
    #keeping the count columns to further calculation 132 columns
    sub_count_norm_b=sub_count_norm_b[ ,c(2:133)]
    #calculate mean ,sd, mean+2sd and mean-2sd per row for each transcript
    mean=rowMeans ( sub_count_norm_b, na.rm=TRUE)
    sd=rowSds(sub_count_norm_b, na.rm=TRUE)
    x=mean-(2*sd)
    y=mean+(2*sd)
    #excluding outliers
    #values above man+2sd and below mean_2sd changes to NA
    sub_count_norm_b[sub_count_norm_b<x]<-NA
    sub_count_norm_b[sub_count_norm_b>y]<-NA
    #ci=(ni-mean)/sd
    sub_count_norm_b<- (sub_count_norm_b-mean)/sd
    
    
    
    #Creating count matrix
    count_matrix=matrix(data=NA, nrow=132, ncol=132, byrow=TRUE)
    for (i in 1:132)
    {
      
      for( j in 1:132)  
      {
        #|ci-c2| #pairwise comparison of ci
        a=abs(sub_count_norm_b[i]-sub_count_norm_b[j])
        count_matrix [i, j]=as.numeric(a)
        #132 132
      }
    }
    
    
    count_triangle=upperTriangle(count_matrix, diag=FALSE)
    
    
    #Some of the count matrics have no data and all are NA and they create an error
    if (sum(is.na(count_triangle))=="8646") { 
      
      next}
    
    else{
      
      
      if (t==1) {
        fit= lm (genome_triangle ~ count_triangle)
        rsquared <- summary(fit)$r.squared
        P_value<- summary(fit)$coefficients[2,4]
        t_value<-summary(fit)$coefficients[2,3]
        Beta<- summary(fit)$coefficients[2,1]
        df<-fit$df.residual
        Transcript_ID=refFlat.new_nodup_c[t,1]
      } 
      
      else {
        fit= lm (genome_triangle ~ count_triangle)
        rsquared_b <- summary(fit)$r.squared
        
        
        #skipping non entry outputs , since they create an error
        if ( nrow(summary(fit)$coefficients)==1 ) 
        {next
        }
        else {
          
          P_value_b<- summary(fit)$coefficients[2,4]
          
          Beta_b<- summary(fit)$coefficients[2,1]
          df_b<-fit$df.residual  
          t_value_b<-summary(fit)$coefficients[2,3]
          t_value=c(t_value_b, t_value)
          rsquared=c (rsquared, rsquared_b)
          P_value= c (P_value_b, P_value)
          Beta= c(Beta_b, Beta)
          df= c(df_b, df)
          Transcript_ID_b=refFlat.new_nodup_c[t,1]
          Transcript_ID=c(Transcript_ID_b, Transcript_ID)
          fit3_100bp_1_1000=data.frame(Transcript_ID, Beta,t_value,  P_value, rsquared, df )
          
          
          ###############################################
          #Checking the linear model assumptions
          ###############################################
          #install.packages("car")
          #library(car)
          #We are checking if the residuals are normally distributed
          #qqPlot(fit)
          #we rae checking if constant error variance holds
          
          #library(car) 
          #ncvTest(fit)
          #spreadLevelPlot(fit)
          
          #############################################
          
          
          #capture.output(u, file=sprintf("/Users/saeedehazary/Documents/R_files/fit_250/%s_fit", refFlat.new_nodup_c[t,1]))
          
          
          
          print(t)
          
          #Data Control
          #if (refFlat.new_nodup_c[t,1]!=st_count_e[t,1]) { 
          #print("There is an Error in data")
          #}
          #else {
          
          # }
          
        }
      }
    }
  }
}









write.table(fit_new_edgeR_250bp_1_1000, file="/Users/saeedehazary/Documents/R_files/fit_250/fit_new_edgeR_250bp_1_1000.txt", sep="\t", row.names=F)






