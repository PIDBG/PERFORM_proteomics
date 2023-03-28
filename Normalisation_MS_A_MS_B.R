############################################################################################################################################################
# ------------- This code normalises and filters the untargeted MS-based proteomic datasets (MS-A and MS-B)
############################################################################################################################################################

library(ggplot2)
library(grid)
library(gridExtra)
library(gplots)
library(stringr)
library(imputeLCMD)
library(limma)
library(stringr)
library(pcaMethods)

# read MS-A cohort 
protein_A <- read.delim(file = 'proteinGroups-A.txt')
rownames(protein_A) <- protein_A$Protein.IDs
A.idx <- read.delim(file = 'protein_info_A.csv')

# remove potential contaminants
protein_A <- protein_A[!(protein_A$Potential.contaminant== "+"),]
A.idx$LFQ <- paste("LFQ.intensity.", A.idx$proteomic.id, sep = "")

# log MS-A data
lfq.d <- protein_A[,1190:1356]	
lfq.d[lfq.d == 0] <- NA

# read MS-B cohort 
protein_B <- read.delim(file = 'proteinGroups.txt')
rownames(protein_B) <- protein_B$Protein.IDs

# remove contaminants 
protein_B <- protein_B[!(protein_B$Potential.contaminant=="+"),]

# log MS-B data
lfq.v <- protein.val[,1554:1772]	
lfq.v[lfq.v == 0] <- NA

# how many proteins are found in both? 
A.pros <- rownames(protein_A)
B.pros <- rownames(protein_B)

length(A.pros[A.pros %in% B.pros])

# remove features used in the 12 depletion column 
# top 12: Albumin, IgG, alpha1-acid glycoprotein, alpha1-antitrypsin, alpha2-macroglobulin, apoliproteins A1 and A2, 
# fibrinogen, haptoglobin, IgA, IgM, Transferrin 
uni.12 <- c('P02768',
            'P01859', 'P01857', 
            'P02763', 
            'P01009',
            'P01023', 
            'P02647', 
            'P02652', 
            'P02679', 'P02671', 'P02675', 
            'P00738', 
            'P01876', 'P01877', 
            'P01871', 
            'P02787')

uni.A <- lapply(uni.12, function(x){
  pro.name <- rownames(lfq.d)[grep(pattern = x, x = rownames(lfq.d))]
  return(pro.name)
})
uni.A <- unlist(uni.A)
lfq.d <- lfq.d[!(rownames(lfq.d) %in% uni.A),]

uni.B <- lapply(uni.12, function(x){
  pro.name <- rownames(lfq.v)[grep(pattern = x, x = rownames(lfq.v))]
  return(pro.name)
})
uni.B <- unlist(uni.B)
lfq.v <- lfq.v[!(rownames(lfq.v) %in% uni.B),]


# frequency of missingness per feature - MS-A 
A.na <- lapply(rownames(lfq.d), function(x){
  exp <- lfq.d[x,1:150]
  na.len <- length(exp[is.na(exp)])
  prop <- 0.90*ncol(lfq.d[,1:150])
  if(na.len>prop){return(x)}
})
A.na <- unlist(A.na)

# remove features with over 90% missingness
# also remove pool samples 
lfq.d.filter <- lfq.d[!(rownames(lfq.d) %in% A.na),1:150]

# frequency of missingness per feature - MS-B 
lfq.v <- lfq.v[,c(1:135,137:202)]
B.na <- lapply(rownames(lfq.v), function(x){
  exp <- lfq.v[x,]
  na.len <- length(exp[is.na(exp)])
  prop <- 0.90*ncol(lfq.v)
  if(na.len>prop){return(x)}
})
B.na <- unlist(B.na)

# remove features with over 90% missingness
lfq.v.filter <- lfq.v[!(rownames(lfq.v) %in% B.na),]

length(rownames(lfq.d.filter)[rownames(lfq.d.filter) %in% rownames(lfq.v.filter)])
dim(lfq.d.filter)
dim(lfq.v.filter)

# pca to identify sample outliers 
# log data first 
# MS-A
log.d <- log2(lfq.d.filter)
# MS-B
log.v <- log2(lfq.v.filter)

# calculate PCs 
# MS-A
A.idx <- A.idx[match(colnames(log.d), A.idx$LFQ),]
pca.d <- pca(t(log.d), nPcs = 10, center = TRUE, scale = "uv")
pca.d <- data.frame(pca.d@scores[,1:2], A.idx$broadDiagnosis)
pca.plot.d <- ggplot(pca.d, aes(x = PC1, y = PC2, color = A.idx.broadDiagnosis))+geom_point()+
  geom_text(label=rownames(pca.d))+stat_ellipse()
pca.plot.d

# remove 2 outliers 
out.d <- rownames(pca.d)[pca.d$PC1 < -30]

# calculate PCs 
# MS-B
colnames(log.v)[13] <- "LFQ.intensity.11"
validation.idx <- validation.idx[match(colnames(log.v), validation.idx$LFQ),]
pca.v <- pca(t(log.v), nPcs = 10, center = TRUE, scale = "uv")
pca.v <- data.frame(pca.v@scores[,1:2],validation.idx$gp)
pca.plot.v <- ggplot(pca.v, aes(x = PC1, y = PC2, color = validation.idx.gp))+geom_point()+
  geom_text(label=rownames(pca.v))+stat_ellipse()
pca.plot.v

