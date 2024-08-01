##############################################################################################
####################### this script will fit GBLUP and GFBLUP models##########################
##############################################################################################

# kinship transform function
KinshipTransform <- function(Matrix) {
  nn <- ncol(Matrix)
  return((sum(diag(Matrix))-as.numeric(matrix(1,1,nn) %*% Matrix %*% matrix(1,nn,1))/nn)/(nn-1))
}

###############################
####load necessary packages####
library(sommer)
library(caret)
library(tidyverse)

###############################

####load data#### 
go.start <- 1 
go.end <- 2 #go start and go end will help grab the correct kinship data if multiple scripts were used to calculate kinships.
load(paste0(go.start,"_",go.end,".Rdata")) #kinships
cyverse.go.ids <- names(selected.kinships) #to store the names of go terms
grain_yield <-  read.csv("grain.csv") #hybrid IDs in rows and phenotypes in columns.
platform_traits <- read.csv("platform_traits.csv")  #hybrid IDs in rows and phenotypes in columns.
gen.groups <- read.csv("gengroups.csv") #the file contains information about hybrids and their probable genetic groups.
K <- readRDS("K.rds") #genomic kinship matrix with the whole SNP panel.
cyverse.go <- readRDS("cyversego.rds") #cyverse.go contains coordinates of genes and their corresponding go terms.
###############################
####assign phenotypic data####

pheno <- grain_yield #here we will chose the traits we are interested in. platform or field grain yield.
keep=intersect(pheno$GENOTYPE,gen.groups$GENOTYPE) # keep has the common hybrids between pheno and gen.groups
gen.groups_new <- gen.groups[match(keep,gen.groups$GENOTYPE),] #only hybrids in keep are kept.
rownames(gen.groups_new)<- 1:nrow(gen.groups_new)
selected_pheno <- pheno[match(keep,pheno$GENOTYPE),] #only hybrids in keep are kept.
rownames(selected_pheno) <- 1:nrow(selected_pheno)
K <- K[match(keep,rownames(K)),match(keep,rownames(K))]
K <- K/KinshipTransform(K)   # this is your scaled matrix
####################################
####creating files for processing and storing results####

final_results <- data.frame(matrix(data = NA, nrow = 1,ncol=6)) # to store the final results with all go terms, trait columns, and three genetic groups.
colnames(final_results) <- c("GBLUP","GF_BLUP","Experiment","Gain","GO_Term","pred_group")


###############################################
####start of different loops####
loop.start <-  1
loop.end <-   length(names(selected.kinships))
## i will guide GO terms selection
for (i in loop.start:loop.end){
  selected.go.id <- cyverse.go.ids[i]
  #selected.go.annot <- cyverse.go[cyverse.go$goterm==selected.go.id,]
      ###storing results
    
    KS <- selected.kinships[[selected.go.id]]
    KS <- KS[match(keep,rownames(KS)),match(keep,rownames(KS))]
    KS <- KS/KinshipTransform(KS)   # this is your scaled matrix
    
    KR <- remain.kinships[[selected.go.id]]
    KR <- KR[match(keep,rownames(KR)),match(keep,rownames(KR))]
    KR <- KR/KinshipTransform(KR)   # this is your scaled matrix
    
    predictive_abilities_all_pheno_col <- data.frame(matrix(data = NA, nrow = 1,ncol=6)) #to store results for all columns of phenotypes
    colnames(predictive_abilities_all_pheno_col) <- c("GBLUP","GF_BLUP","Experiment","Gain","GO_Term","pred_group")
    
    ## a will help choosing different columns of phenotypes one by one
  for (a in 2:ncol(selected_pheno)) {
      phen <- data.frame(selected_pheno[,1],selected_pheno[,1],selected_pheno[,a]) #sommer requires two ID columns to fit two random effects (f and r)
      colnames(phen) <- c("id","id2","y")
      
      
      predictive_abilities <- data.frame(matrix(data = NA, nrow = 3,ncol=6))
      colnames(predictive_abilities) <- c("GBLUP","GF_BLUP","Experiment","Gain","GO_Term","pred_group")
      ## b will help in navigating through different population groups as well as storing predictive abilities
      for (b in 1:3) {
        
        phen_na <- phen
        index <- gen.groups_new[gen.groups_new$GROUP_NUM==b,] #the index of hybrids that are to be masked.
        index <- as.numeric(rownames(index))
        phen_na$y[index] <- NA
        
        
        mod_gblup <- mmer(y~1,random = ~vsr(id,Gu=K),rcov = ~units,data = phen_na, verbose = F,tolParInv = 1)
        if(is.vector(mod_gblup)){predictive_abilities$GF_BLUP[b] <- NA} else {
          pred_gblup <- mod_gblup$U$`u:id`$y
          pred_gblup <- pred_gblup[match(as.character(phen$id),as.character(names(pred_gblup)))]#sommer tends to reorder characters. we do to this make sure the order is good.
          predictive_abilities$GBLUP[b] <- cor(pred_gblup[index],phen[index,3],use = "na.or.complete") #pearson cor=predictive ability
          }
        
        mod_gfblup <- mmer(y~1,random = ~vsr(id,Gu=KS)+vsr(id2,Gu=KR),rcov = ~units,data = phen_na, verbose = F,tolParInv = 1)
        if(is.vector(mod_gfblup)){predictive_abilities$GF_BLUP[b] <- NA} else {
          pred_gfblup <- mod_gfblup$U$`u:id`$y + mod_gfblup$U$`u:id2`$y
          pred_gfblup <- pred_gfblup[match(as.character(phen$id),as.character(names(pred_gfblup)))]
          predictive_abilities$GF_BLUP[b] <- cor(pred_gfblup[index],phen[index,3],use = "na.or.complete")
          }
        rm(mod_gblup,pred_gblup,mod_gfblup,pred_gfblup)
        
        predictive_abilities$Experiment[b] <- colnames(selected_pheno)[a]  #to store the environment name/trait name in the output
        predictive_abilities$Gain[b] <- predictive_abilities$GF_BLUP[b]-predictive_abilities$GBLUP[b] #to store the gain
        predictive_abilities$GO_Term[b] <- selected.go.id # to store GO id
        predictive_abilities$pred_group[b] <- paste0("Predicted-Group:",b) # to store the id of the genetic group predicted.
        }
      
      predictive_abilities_all_pheno_col <- rbind(predictive_abilities_all_pheno_col,predictive_abilities) #to store results of all pheno columns for a given go term
      
    }
    predictive_abilities_all_pheno_col <- predictive_abilities_all_pheno_col[-1,] #the first row was created with NA to make this variable. now it is removed.
    
    ####running models with cross validation####
    # k will help select columns in phenotypic data
    
    final_results <- rbind(final_results,predictive_abilities_all_pheno_col)
 
      print(paste("Status:",i))
  
}


final_results <- final_results[-1,] #the first row was created with NA to make this variable. now it is removed.

saveRDS(final_results,file = paste0("res_gfblup_",go.start,"_",go.end,".rds"))

####end####















