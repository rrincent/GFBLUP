#################################################################################
######### this script will be used to apply Wilcoxon Mann-Whitney Test###########
#################################################################################
#load libraries
library(tidyverse)
library(dplyr)
###############################################################################
#load Data
godes <- data.table::fread("GO Terms and Descriptions from Ensembl Plants.txt")
all_pop_results <- readRDS("all_pop_results_supp.rds")
all_pop_results <- all_pop_results[all_pop_results$GO_Term %in% godes$`GO term accession`,]

######################################################
#create variables for processing and storing results
terms <- unique(godes$`GO term accession`)
terms <- intersect(terms,all_pop_results$GO_Term)
p.results <- data.frame(matrix(data = NA,nrow = 1,ncol = 2))
colnames(p.results) <- c("GO_term","p.value")

##############################
#running loops to apply MWT

for (i in 1:length(terms)) { #i will guide through all go terms one by one
  sel.term <- terms[i] 
  sel.data <- all_pop_results[all_pop_results$GO_Term %in% sel.term,] # get predictive abilities of a given term
  
  results <- data.frame(matrix(data = NA,nrow = 1,ncol = 2)) #an empty data frame to store results
  colnames(results) <- c("GO_term","p.value")
  
  results$GO_term <- sel.term #store the name of the go term
  
  #formatting data for the test
  test.data <- data.frame(sel.data$GBLUP,sel.data$GF_BLUP) 
  test.data <- na.omit(test.data)
  colnames(test.data) <- c("GBLUP","GF_BLUP")
  
  #Wilcox mann whitney test
  test.result <- wilcox.test(x = test.data$GBLUP,y=test.data$GF_BLUP, paired=T,alternative=c("two.sided"))
  results$p.value <- test.result$p.value
  p.results <- rbind(p.results,results)
  
  print(i)  
}

p.results <- p.results[-1,] #removing the empty line created in the start

term_gains <- aggregate(Gain~GO_Term,all_pop_results,mean) # to get mean gain per term
term_gains <- term_gains[match(p.results$GO_term,term_gains$GO_Term),] #ensures the same order of p.results and terms_gains variables
p.results$Gain <- term_gains$Gain
sel_godes <- godes[match(p.results$GO_term,godes$`GO term accession`),]
p.results$GO_Description <- sel_godes$`GO term name`

#plot pvalues
ggplot(p.results,aes(x = p.value))+geom_histogram(color="black",fill="skyblue3",bins = 50)+
  ggtitle("Grain Yield")+ylab("Frequency")+xlab("Raw P-values")+
  theme(text = element_text(face = "bold"))



#function for multiple testing correction of pvalues
CorrectingPValues <- function(p){
  
  ## Get a p0 estimate
  p0 = min(2*sum(p > 0.5,na.rm = T)/length(p),1-1/length(p)) 
  
  ## Get the corrected p-values
  CorrectedPval <-  p.adjust(p, method = 'BH')*p0
  
  return(CorrectedPval)
}

Alpha <- 0.05
Signif <- CorrectingPValues(p.results$p.value) %>% 
  `<`(Alpha) %>% 
  which

Threshold <- p.results[Signif,]$p.value %>% max() #5% fdr threshold



sig_terms <- p.results[Signif,] #significant terms
saveRDS(sig_terms,file = "sig_terms.rds")

####end####