############################################################################################################################################################
# ------------- Code for the discovery of protein biomarkers for validation in PERFORM study 
############################################################################################################################################################

library(ggplot2)
library(grid)
library(stringr)
library(limma)
library(clipr)
library(nnet)
library(caTools)

top_signature_iterate <- function(sig_matrix){
  # take off AUCs
  col.1 <- ncol(sig_matrix)-1
  # collapse rows into one string 
  concat <- apply(sig_matrix[,1:col.1], 1, function(x){
    paste(x, collapse = "_")
  })
  # calculate how many times each string appears
  sig.max <- str_split(names(which(table(concat) == max(table(concat)))), "\\_")
  # return most frequent signature
  sig.max <- data.frame(ID = colnames(sig_matrix)[1:col.1], Present = data.frame(sig.max)[[1]])
  sig.max <- sig.max[sig.max$Present == 1,]
  return(sig.max)
}
# load MS-A and MS-B
load(file = "MSA_MSB.RData", verbose = T)
# load SomaScan
load(file = 'qcExp.RData', verbose = T)
# load SomaScan protein information
load(file = 'proteinData.RData', verbose = T)

# path to fspls 
path.fspls <- paste("fspls_lars_multi_meta.R")

############################################################################################################################################################
# ------------- Limma for differential abundance analysis
############################################################################################################################################################

# MS-A
# design matrix 
mat.A <- model.matrix(~0 + factor(broadDiagnosis) + factor(sex) +  as.numeric(ageMonths), data = meta.a)
colnames(mat.A) <- c("DB", "DV", "M", "age")

# remove CRP 
log.d <- log.d[!rownames(log.d) == 'sp|P02741|CRP_HUMAN',]

limma.a.all <- limma.fun(data = t(log.d), p.thresh = 1, c("DB", "DV"), start.col = 1,  
                                 end.col = nrow(log.d), n.features = nrow(log.d), model.mat = mat.dis)

limma.a.all$Color <- 'black'
limma.a.all$Color[abs(limma.a.all$logFC) > 1] <- 'limegreen'
limma.a.all$Color[limma.a.all$adj.P.Val < 0.05] <- 'gold'
limma.a.all$Color[limma.a.all$adj.P.Val < 0.05 & 
                            abs(limma.a.all$logFC) > 1] <- 'red'

limma.a.all$Color <- factor(limma.a.all$Color, levels = c("black", 'limegreen', "gold", 'red'))

# volcano plot for MS-A
volcano_MSA <- ggplot(limma.a.all, aes(x = logFC, y = -log10(adj.P.Val), color = Color))+
  geom_point(size = 2)+
  theme_bw()+
  scale_x_continuous(limits = c(-4, 4))+
  scale_color_manual(values = c("black", 'limegreen', 'gold', 'red'), 
                     labels = c('P-value > 0.05',
                                expression('absolute log'[2]*'FC > 1'),
                                "P-value < 0.05", 
                                expression('P-value < 0.05 and absolute log'[2]*'FC > 1')))+
  labs(y= "-log10 adjusted p-value", 
       color = "", 
       title = "MS-A cohort")+
  geom_label_repel(
    data = subset(limma.a.all, Color == 'red' | 
                    Color == 'gold' & -log10(adj.P.Val) > 2.5),
    aes(label = gene),
    fill = 'white',
    size = 2.5,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.25, "lines"), 
    color = 'black',
    max.overlaps = 15)+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = -1, linetype = 'dashed')+
  theme(legend.text.align = 0, 
        legend.position = 'right', 
        plot.title = element_text(hjust = 0.5))

# MS-B
# design matrix 
mat.b <- model.matrix(~0 + factor(gp) + factor(sex) + as.numeric(age), data = b.idx)

colnames(mat.b) <- c("DB", "DV", "M", "age")

log.v <- log.v[!rownames(log.v) == 'sp|P02741|CRP_HUMAN',]

limma.b.all <- limma.fun(data = t(log.v), p.thresh = 1, c("DB", "DV"), start.col = 1,  
                           end.col = nrow(log.v), n.features = nrow(log.v), model.mat = mat.val)

# volcano plot for MS-B
limma.b.all$Color <- 'black'
limma.b.all$Color[abs(limma.b.all$logFC) > 1] <- 'limegreen'
limma.b.all$Color[limma.b.all$adj.P.Val < 0.05] <- 'gold'
limma.b.all$Color[limma.b.all$adj.P.Val < 0.05 & 
                      abs(limma.b.all$logFC) > 1] <- 'red'

limma.b.all$Color <- factor(limma.b.all$Color, levels = c("black", 'limegreen', "gold", 'red'))

volcano_MSB <- ggplot(limma.b.all, aes(x = logFC, y = -log10(adj.P.Val), color = Color))+
  geom_point(size = 2)+
  theme_bw()+
  scale_x_continuous(limits = c(-3, 3))+
  scale_color_manual(values = c("black", 'limegreen', 'gold', 'red'), 
                     labels = c('NS',
                                expression('log'[2]*'FC'),
                                "P-value", 
                                expression('P-value and log'[2]*'FC')))+
  labs(y= "-log10 adjusted p-value", 
       color = "", 
       title = "MS-B cohort")+
  geom_label_repel(
    data = subset(limma.b.all, Color == 'red' | 
                    Color == 'gold' & -log10(adj.P.Val) > 3),
    aes(label = gene),
    fill = 'white',
    size = 2.5,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.25, "lines"), 
    color = 'black',
    max.overlaps = 15)+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = -1, linetype = 'dashed')+
  theme(legend.text.align = 0, 
        legend.position = 'none', 
        plot.title = element_text(hjust = 0.5))


# DE proteins in Somascan 
model.soma <- model.matrix(~0 + as.factor(qcExp$SampleGroup) + as.numeric(as.character(qcExp$ageMonths)) + as.factor(qcExp$Gender) + as.factor(qcExp$PlateId))
colnames(model.soma) <- c("DB", 'DV', 'age', 'M', 'plate2', 'plate3')

limma_all <- limma.fun(data = qcExp, p.thresh = 0.1, c("DB", "DV"), start.col = 33,  
                       end.col = 1332, n.features = ncol(qcExp), model.mat = model.soma)

# add uniprot info 
limma$uniprot <- uniprot.df(limma)
limma$name <- names.df(limma)

# volcano plot for somascan
limma_all$gene <- proteinData$EntrezGeneSymbol[match(rownames(limma_all), proteinData$SeqId)]

limma_all$Color <- 'black'
limma_all$Color[abs(limma_all$logFC) > 1] <- 'limegreen'
limma_all$Color[limma_all$adj.P.Val < 0.05] <- 'gold'
limma_all$Color[limma_all$adj.P.Val < 0.05 & 
                  abs(limma_all$logFC) > 1] <- 'red'

limma_all$Color <- factor(limma_all$Color, levels = c("black", "gold", 'red'))

# volcano plot 
volcano_somascan <- ggplot(limma_all, aes(x = logFC, y = -log10(adj.P.Val), color = Color))+
  geom_point(size = 2)+
  theme_bw()+
  scale_x_continuous(limits = c(-3, 3))+
  scale_color_manual(values = c("black", 'gold', 'red'), 
                     labels = c('NS',
                                expression('log'[2]*'FC'),
                                "P-value", 
                                expression('P-value and log'[2]*'FC')))+
  labs(y= "-log10 adjusted p-value", 
       color = "", 
       title = "SomaScan cohort")+
  geom_label_repel(
    data = subset(limma_all, Color == 'red' | 
                    Color == 'gold' & -log10(adj.P.Val) > 3 ),
    aes(label = gene),
    fill = 'white',
    size = 2.5,
    box.padding = unit(0.25, "lines"),
    point.padding = unit(0.25, "lines"), 
    color = 'black',
    max.overlaps = 15)+
  geom_hline(yintercept = -log10(0.05), linetype = 'dashed')+
  geom_vline(xintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = -1, linetype = 'dashed')+
  theme(legend.text.align = 0, 
        legend.position = 'none', 
        plot.title = element_text(hjust = 0.5))

###################### what is the intersection between datasets?

# summarise the overlaps between different datasets 
# overlap between somascan and untargeted cohorts
soma.a <- lapply(limma$uniprot, function(x){
  a.match <- rownames(limma.a.sig)[grep(pattern = x, x = rownames(limma.a.sig))]
  print(a.match)
})
soma.a <- unlist(soma.a)
length(soma.a)

soma.b <- lapply(limma$uniprot, function(x){
  b.match <- rownames(limma.b.sig)[grep(pattern = x, x = rownames(limma.b.sig))]
  print(b.match)
})
soma.b <- unlist(soma.b)
length(soma.b)

soma.a.b <- intersect(soma.a, soma.b)
length(soma.a.b)


# for the overlapping proteins, make sure the LFC goes in the same direction 

# take out stats for proteins found in all 3 datasets
soma.a.b.df <- data.frame(limma.a.sig[match(soma.a.b, rownames(limma.a.sig)),][,c(1,5)],
                              limma.b.sig[match(soma.a.b, rownames(limma.b.sig)),][,c(1,5)])
colnames(soma.a.b.df) <- c('A.LFC', 'A.Adj.P', 'B.LFC', 'B.Adj.P')

soma.a.b.df$Soma.LFC <- limma_all$LFC[match(soma.a.b.df$Uniprot, limma_all$Uniprot)]
soma.a.b.df$Soma.P <- limma_all$padj[match(soma.a.b.df$Uniprot, limma_all$Uniprot)]

dim((soma.a.b.df$Soma.LFC > 0 & soma.a.b.df$A.LFC > 0 & soma.a.b.df$B.LFC > 0) |
      (soma.a.b.df$Soma.LFC < 0 & soma.a.b.df$A.LFC < 0 & soma.a.b.df$B.LFC < 0))


############################################################################################################################################################
# ------------- Run FSPLS to identify signatures
############################################################################################################################################################

# SomaScan first
soma.exp <- qcExp[qcExp$SampleGroup=="DB" | qcExp$SampleGroup == "DV",]
soma.response <- class.ind(soma.exp$SampleGroup)
soma.exp <- soma.exp[,33:ncol(soma.exp)]
load.packages()
fspls.somascan <- fspls.iterate(n.iterations = 100, seed = 1, path.fspls = path.fspls, 
                                expression = soma.exp, response = soma.response[,1], 
                                p.thresh = 0.001, max = 5, beam = 10, split = T, split.ratio = 0.7)
# upset plot 
signatures.soma <- fspls.somascan$table
signatures.soma <- upset.soma[,c(1:3,5:ncol(upset.soma), 4)]
signatures.soma[is.na(signatures.soma)] <- 0
olnames(signatures.soma) <- proteinData$EntrezGeneSymbol[match(colnames(signatures.soma), proteinData$SeqId)]

top_soma <- top_signature_iterate(signatures.soma)
colSums(signatures.soma)[top_soma$ID]

################################################################################################
# MS-A Discovery 
meta.a.bv <- meta.a[meta.a$broadDiagnosis == "DB" | meta.a$broadDiagnosis == "DV",]
a.exp <- data.frame(t(log.d[,colnames(log.d) %in% meta.a.bv$LFQ]))
a.response <- class.ind(meta.a.bv$broadDiagnosis)

# FSPLS with iterations for MS-A cohort 
fspls.a <- fspls.iterate(n.iterations = 100, seed = 1, path.fspls = path.fspls, 
                                 expression = a.exp, response = a.response[,1], 
                                 p.thresh = 0.001, max = 5, beam = 10, split = T, split.ratio = 0.7)

upset.a <- fspls.discovery$table
upset.a <- upset.a[,c(1:3, 5:ncol(upset.a), 4)]
upset.a[is.na(upset.a)] <- 0

top_a <- top_signature_iterate(upset.a)
colSums(upset.a)[top_a$ID]

################################################################################################
# MS-B
b.idx.bv <- b.idx[b.idx$gp == "DB" | b.idx$gp == "DV",]
b.exp <- data.frame(t(log.v[,colnames(log.v) %in% b.idx$LFQ]))
b.response <- class.ind(b.idx.bv$gp)

# FSPLS with iterations for MS-B cohort 
fspls.b <- fspls.iterate(n.iterations = 100, seed = 123, path.fspls = path.fspls, 
                                  expression = b.exp, response = b.response[,2], 
                                  p.thresh = 0.001, max = 5, beam = 10, split = T, split.ratio = 0.7)

upset.b <- fspls.validation$table
upset.b <- upset.b[,c(1:4, 6:ncol(upset.b), 5)]
upset.b[is.na(upset.b)] <- 0

top_b <- top_signature_iterate(upset.b)
colSums(upset.b)[top_b$ID]


################################################################
# robustness of the proteins 
somascan.proteins <- unlist(str_split(fspls.somascan$sig, "__"))
somascan.proteins <- str_remove_all(somascan.proteins, "\\.test")
frequency.somascan <-data.frame(table(somascan.proteins))
frequency.somascan <- frequency.somascan[order(frequency.somascan$Freq),]

colnames(upset.soma) <- proteinData$SeqId[match(colnames(upset.soma), proteinData$EntrezGeneSymbol)]
pasted <- apply(upset.soma[,1:ncol(upset.soma)-1], 1, function(x){
  pasted <- paste(as.character(x), collapse = '-')
  print(pasted)
})

dups <- data.frame(pasted[duplicated(pasted)])
tab <- table(pasted[pasted %in% dups$pasted.duplicated.pasted..])

top.soma <- names(tab)[which(tab == max(tab))]
soma.auc <- data.frame(auc = upset.soma[,ncol(upset.soma)][pasted %in% top.soma], paste = pasted[pasted %in% top.soma])

if(length(top.soma) > 1){
  means.soma <- lapply(top.soma, function(x){
    return(mean(soma.auc$auc[soma.auc$paste == x]))
  })
  top.soma.sig <- top.soma[which(unlist(means.soma) == max(unlist(means.soma)))]}else(top.dis.sig <- top.soma)

upset.top.hit.soma <- upset.soma[which(pasted == top.soma.sig),]
features.soma <- names(upset.top.hit.soma[1,upset.top.hit.soma[1,] == 1])
soma.features.freq <- frequency.somascan[frequency.somascan$somascan.proteins %in% features.soma,]
soma.features.freq$name <- proteinData$TargetFullName[match(soma.features.freq$somascan.proteins, proteinData$SeqId)]

# MS-A frequency 
a.proteins.freq <- unlist(str_split(fspls.a$sig, "__"))
a.proteins.freq <- str_remove_all(a.proteins.freq, "\\.test")
frequency.a <-data.frame(table(a.proteins.freq))
frequency.a <- frequency.dis[order(frequency.a$Freq),]

colnames(upset.a) <- names.a.proteins
pasted.a <- apply(upset.a[,1:ncol(upset.a)-1], 1, function(x){
  pasted <- paste(as.character(x), collapse = '-')
  print(pasted)
})

dups.a <- data.frame(pasted.a[duplicated(pasted.a)])
tab.a <- table(pasted.a[pasted.a %in% dups.a$pasted.a.duplicated.pasted.a..])
top.a <- names(tab.a)[which(tab.a == max(tab.a))]
a.auc <- data.frame(auc = upset.a[,ncol(upset.a)][pasted.a %in% top.a], paste = pasted.a[pasted.a %in% top.a])

if(length(a.auc$auc) > 1){
  means.a <- lapply(top.a, function(x){
    return(mean(a.auc$auc[a.auc$paste == x]))
  })
  top.a.sig <- top.a[which(unlist(means.a) == max(unlist(means.a)))]}else(top.a.sig <- top.a)

upset.top.hit.a <- upset.a[which(pasted.a == top.a.sig),]
features.a <- names(upset.top.hit.a[1,upset.top.hit.a[1,] == 1])
a.features.freq <- frequency.a[frequency.a$a.proteins.freq %in% features.a,]
a.features.freq$name <- a.proteins$ID[match(a.features.freq$a.proteins.freq, rownames(a.proteins))]

# MS-B frequency 
b.proteins.freq <- unlist(str_split(fspls.b$sig, "__"))
b.proteins.freq <- str_remove_all(b.proteins.freq, "\\.test")
frequency.b <-data.frame(table(b.proteins.freq))
frequency.b <- frequency.vb[order(frequency.b$Freq),]

colnames(upset.b) <- names.b.proteins
pasted.b <- apply(upset.b[,1:ncol(upset.b)-1], 1, function(x){
  pasted <- paste(as.character(x), collapse = '-')
  print(pasted)
})

dups.b <- data.frame(pasted.b[duplicated(pasted.b)])
tab.b <- table(pasted.b[pasted.b %in% dups.b$pasted.b.duplicated.pasted.b..])
top.b <- names(tab.b)[which(tab.b == max(tab.b))]
b.auc <- data.frame(auc = upset.b[,ncol(upset.b)][pasted.b %in% top.b], paste = pasted.b[pasted.b %in% top.b])

if(length(b.auc$auc) > 1){
  means.b <- lapply(top.b, function(x){
    return(mean(b.auc$auc[b.auc$paste == x]))
  })
  top.b.sig <- top.b[which(unlist(means.b) == max(unlist(means.b)))]}else(top.b.sig <- top.b)

upset.top.hit.b <- upset.b[which(pasted.b == top.b.sig),]
features.b <- names(upset.top.hit.val[1,upset.top.hit.val[1,] == 1])
b.features.freq <- frequency.b[frequency.b$b.proteins.freq %in% features.b,]
b.features.freq$name <- b.proteins$ID[match(b.features.freq$b.proteins.freq, rownames(b.proteins))]

# export frequencies 
write.csv(file = "somascan_frequency.csv", soma.features.freq)
write.csv(file = "MS_A_frequency.csv", a.features.freq)
write.csv(file = "MS_B_frequency.csv", b.features.freq)