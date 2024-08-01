###################################################################################
####### this script applies GWAS on traits followed by meta-analysis###############
###################################################################################

#function for scaling matrices
KinshipTransform <- function(Matrix) {
  nn <- ncol(Matrix)
  return((sum(diag(Matrix))-as.numeric(matrix(1,1,nn) %*% Matrix %*% matrix(1,nn,1))/nn)/(nn-1))
}

##############################
#load packages#
library(tidyverse);library(MM4LMM)
library(metaGE);library(corrplot)
##################################################
# loading phenotypic data and its transformation##
keep <- read.csv("gfblup_lines_ids.csv") # these are the names of the hybrids in three genetic groups
keep <- keep$GENOTYPE
phenogy <-  read.csv("grain.csv") #hybrids in rows and traits in columns
rownames(phenogy) <- phenogy$GENOTYPE
phenogy <- phenogy[,-1] #remove the column containing hybrid IDs
phenogy <- phenogy[keep,] #keep the hybrids from three genetic groups
good_lines <- setdiff(rownames(phenogy),rownames(phenogy)[rowSums(is.na(phenogy))>20]) #we removed lines with too many missing values as it affects power
good_envs <- setdiff(colnames(phenogy),colnames(phenogy)[colSums(is.na(phenogy))>20]) #we removed experiments with too many missing values
phenogy <- phenogy[good_lines,good_envs]
keep <- rownames(phenogy) # hybrids with too many missing values are now removed.
Ngeno <- length(keep)

#traits/experiments are renamed as they had too long names
phenoenvs <- NULL
phenoenvs$Experiment <- colnames(phenogy)
phenoenvs <- data.frame(phenoenvs)
phenoenvs$Experiment <- gsub(".grain_yield_15","",phenoenvs$Experiment)
phenoenvs <- phenoenvs %>% mutate(env=paste(substr(Experiment,1,3), 
                                            substr(Experiment,str_locate(Experiment,"20")[,1]+2,str_locate(Experiment,"20")[,1]+3),
                                            substr(Experiment,str_locate(Experiment,"20")[,1]+5,str_locate(Experiment,"20")[,1]+10),sep=""))

phenoenvs$env[which(substr(phenoenvs$env,6,10)=="OPT")] = paste(substr(phenoenvs$env[which(substr(phenoenvs$env,6,10)=="OPT")],1,5),"W",sep="")
phenoenvs$env[which(substr(phenoenvs$env,6,10)=="WD")] = paste(substr(phenoenvs$env[which(substr(phenoenvs$env,6,10)=="WD")],1,5),"R",sep="")

colnames(phenogy) <- phenoenvs$env #shortened names are now added to the phenogy 

#load genomic kinship matrix
kinship <- readRDS("K.rds")
kinship <- kinship[keep,keep]
kinship <- kinship/KinshipTransform(kinship) #scaled kinship matrix

# Genotyping data
# The genotyping data are necessary to estimate the genetic correlation btw envts for metaGE

geno <- readRDS("genomat.rds")
#geno=geno[seq(1,nrow(geno),by=100),] # Keep only one out of 100 SNP
dim(geno)
geno = geno[,match(keep,colnames(geno))]
snp_name =rownames(geno)

freq=apply(geno,1,mean)/2
summary(freq)

remove=which(freq<0.05|freq>0.95) #removing alleles that are more or less fixed
length(remove)

geno=geno[-remove,]


het=apply(geno,1,function(a) sum(a==1))/Ngeno
summary(het)

remove=which(het>0.15) 
length(remove)

if(length(remove)>0){geno=geno[-remove,]}
  


##################################################################################
# Association mapping (necessary to identify the SNP not associated, to estimate the genetic correlation btw envts, see next step)
##################################################################################


pval_save = matrix(NA,nrow(geno),ncol(phenogy))
rownames(pval_save)=rownames(geno)
colnames(pval_save)=colnames(phenogy)
beta_save = pval_save

for (i in 1:ncol(phenogy))  {
  print(i)  
  
  test= !is.na(phenogy[,i])
  Y <- phenogy[test,i]
  X <- t(geno[,test])
  VL <- list(G=kinship[test,test],Err=diag(sum(test)))
  Res <- MMEst(Y=Y, X=X, VarList = VL, Cofactor = NULL)
  
  Gwas <- map(AnovaTest(Res),~ .x[2,]) %>%
    bind_rows() %>%
    mutate(Marker = names(Res)) %>%
    mutate(Beta = map_dbl(Res, ~ .x$Beta[length(.x$Beta)])) %>%
    mutate(Beta_SE = map_dbl(Res, ~ .x$VarBeta[length(.x$Beta),length(.x$Beta)] %>% sqrt)) %>%
    rename(WaldStatistic = `Wald (Type III)`) %>%
    select(Marker,Beta,Beta_SE,WaldStatistic,df,pval)
  
  pval_save[,i] = Gwas$pval[match(rownames(pval_save),Gwas$Marker)]
  beta_save[,i] = Gwas$Beta[match(rownames(beta_save),Gwas$Marker)]
  
}

##################################################################################
# loading SNP map file
##################################################################################

snp_map <- readRDS("map_v4.rds")
head(snp_map)
dim(snp_map)
snp_map <- snp_map[match(rownames(geno),snp_map$snp),]



################
## experiment wise GWAS results are saved in the directory
freq=apply(geno,1,mean)/2

File.list=NULL

for (i in 1:ncol(phenogy))  {
  print(i)  
  
  data =tibble(Marker_Name=rownames(geno),Chromosome=snp_map$chr,Marker_Position=snp_map$pos,Maf = sapply(freq,function(a) min(a,1-a)),SNP_Weight = beta_save[,i], Pvalue = pval_save[,i] , Allele1 = 1, Allele2 = 2)
  
  
  write_delim(data, file = paste(colnames(phenogy)[i],".txt",sep=""))
  
}


## load GWAS results for meta-analysis

File.list = paste0(colnames(phenogy),".txt")

names(File.list) = colnames(phenogy)

## First provide the list of variable names. these variables will be fetched from each file
Names.list <- list(MARKER="Marker_Name",
                   CHR="Chromosome",
                   POS="Marker_Position",
                   FREQ="Maf",
                   EFFECT="SNP_Weight",
                   PVAL="Pvalue",
                   ALLELE0="Allele1",
                   ALLELE1="Allele2")
MinFreq <- 0.05
## Now fetch the variables specified
MetaData <- metaGE.collect(File.list, Names.list,MinFreq = MinFreq,DropDuplicates = F)



##################################################################################
# Run metaGE on the gwas results to get an estimate of the btw envt genetic correlation matrix
##################################################################################


Threshold <- 0.6 #threshold on posteriors to be H1

MatCorr <- metaGE.cor(MetaData$Data,Threshold = Threshold)

corrplot(MatCorr,order = 'hclust',title = "GWAS")

##############################################
##fixed effect procedure of meta-analysis.

FeDF <- metaGE.fit(MetaData$Data, MatCorr, Method = "Fe",NA.omit = T)

metaGE.pvalplot(FeDF$PVALUE, Main='GWAS: Fixed Effect Procedure')

saveRDS(FeDF,file = "GY_FeDF.rds")

##################################################################################
##random effect procedure of meta-analysis.

ReDF <- metaGE.fit(MetaData$Data, MatCorr, Method = "Re",NA.omit = T)

metaGE.pvalplot(ReDF$PVALUE, Main='GWAS: Random Effect Procedure')


saveRDS(ReDF,file = "GY_ReDF.rds")

####end####

