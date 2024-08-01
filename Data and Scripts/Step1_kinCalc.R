#########################################################################################
##### this script calculates genomic kinship matrices for features and remaining SNPs####
#########################################################################################
#loading libraties
library(statgenGWAS)

#####################
#loading input data
cyverse.go <- readRDS("cyversego.rds") #file containing gene coordinates and the go terms they belong to.
map_v4 <- readRDS("map_v4.rds") #coordinates of SNPs on AGPv4 assembly
geno <- readRDS("genomat.rds") #matrix containing hybrid IDs in rows and SNPs in columns

#####################
#creating variables for running and storing kinship matrices
cyverse.go.ids <- unique(cyverse.go$goterm)

go.start <- 1
go.end <- 2#length(cyverse.go.ids) #go.start and go.end will guide the start and end of the loop. The indexes can specified to calculate kinship matrices for a go term of interest.
#note: It is better to run multiple scripts, because if all the go terms are used in a single script, it will take too long and too much space probably breaking R.

selected.kinships <- list() #list to store kinships from feature SNPs
remain.kinships <- list() #list to store kinships from remaining SNPs
no.snp.found <- vector() #this vector will store IDs of go terms for which no feature SNP is found.
terms.less.than.10.snps <- vector()
snp.count <- data.frame(matrix(data=NA, nrow=go.end,ncol=2)) #snp.count will store the number of feature SNPs for each go term.
colnames(snp.count) <- c("GO","Count")
window_size <- 5000 #the upstream and downstream window of gene

#####################
####start of different loops####
## i will guide GO terms selection and j will guide through each gene linked to a go term
for (i in go.start:go.end){ 
  
  selected.go.id <- cyverse.go.ids[i]
  selected.go.annot <- cyverse.go[cyverse.go$goterm==selected.go.id,]
  snps_in_data <- map_v4[1,]
  snps_in_data[1,] <- NA
  
  #finding SNPs for each gene
  for(j in 1:nrow(selected.go.annot)){
    step1 <- map_v4[map_v4$chr==selected.go.annot$chr[j],]
    snps_present <- step1[(step1$pos>=selected.go.annot$start[j]-window_size)& (step1$pos<=selected.go.annot$end[j]+window_size),]
    snps_in_data <- rbind(snps_in_data,snps_present)
  }
  snps_in_data <- na.omit(snps_in_data)
  
  ##estimating kinship matrices
  
  selected_snps <- unique(snps_in_data$snp)
  snp.count$GO[i] <- selected.go.id
  snp.count$Count[i] <- length(selected_snps)
  
  if (length(selected_snps)==0) {no.snp.found[i]=selected.go.id} else{
    if(length(selected_snps)<10) {terms.less.than.10.snps[i]=selected.go.id} else{
    
      remain_snps <- setdiff(as.character(rownames(geno)),selected_snps)
      selected.kinships[[selected.go.id]] <- statgenGWAS::kinship(t(geno[match(selected_snps,rownames(geno)),]),method = "vanRaden")
      remain.kinships[[selected.go.id]] <- statgenGWAS::kinship(t(geno[match(remain_snps,rownames(geno)),]),method = "vanRaden")
    
    print(i)}
  }
}

no.snp.found <- na.omit(no.snp.found)
terms.less.than.10.snps <- na.omit(terms.less.than.10.snps)
save(no.snp.found,selected.kinships,remain.kinships,snp.count,file = paste0(go.start,"_",go.end,".Rdata"))
####end####