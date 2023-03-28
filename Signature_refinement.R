############################################################################################################################################################
# ------------- Signature refinement stage
# ------------- Analyse levels of proteins quantified using Luminex or ELISA (proteins from high-throughput datasets + literature)
############################################################################################################################################################

load.packages()
library(missMDA)
library(sva)
library(tidyverse)
library(lubridate)
library(sp)
library(caret)
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

roc_plot <- function(roc, title, colour){
  ci <- ci.se(roc, specificities=seq(0, 1, l=25))
  ci <- data.frame(x = as.numeric(rownames(ci)),
                   lower = ci[, 1],
                   upper = ci[, 3])
  p <- ggroc(roc) + theme_bw() + geom_abline(slope=1, intercept = 1, linetype = "dashed", alpha=0.7, color = "grey")+
    labs(x = '1-Specificity', y = "Sensitivity")+
    geom_ribbon(
      data = ci,
      aes(x = x, ymin = lower, ymax = upper),
      alpha = 0.3,
      fill = colour,
      inherit.aes = F) +
    annotate("text", x=0.4, y=0.1, label= paste("AUC: ", 
                                                round(auc(roc), 3)*100, 
                                                "% (95% CI: ",
                                                round(ci(roc)[1], 3)*100,
                                                '%-', 
                                                round(ci(roc)[3], 3)*100, "%)",
                                                sep = ''))+
    theme(plot.title = element_text(size=12))+
    ggtitle(title)
}

# read in metadata
meta_orig <- read.csv(file = 'Meta.csv', header= T, stringsAsFactors = F)

# read in Luminex protein abundances
lumi <- read.csv(file = 'All_Original_Proteins_analysis_ready.csv', header = T)
meta_orig <- meta_orig[match(lumi$Sample.ID, meta_orig$ID),]

# read in ELISA protein abundances
elisas <- read.csv(file = 'ELISA_proteins_analysis_ready.csv', header = T)
elisas <- elisas[elisas$X %in% meta_orig$ID,]

# convert all to numeric
lumi[,3:ncol(lumi)] <- data.frame(apply(lumi[,3:ncol(lumi)], 2, function(x) as.numeric(as.character(x))))

# merge the proteins 
all_proteins <- cbind(lumi,
                      elisas[match(lumi$Sample.ID, 
                                   elisas$X),4:ncol(elisas)])

all_proteins[,3:ncol(all_proteins)] <- data.frame(apply(all_proteins[,3:ncol(all_proteins)], 2, function(x) as.numeric(as.character(x))))


############################################################################################################
# PCA of the data 

# Luminex proteins 
lumi_proteins_pc <- estim_ncpPCA(log2(all_proteins[,c(5:25)]+1), scale = F, method = "Regularized")
lumi_proteins_pc$ncp
res_lumi_proteins <- imputePCA(log2(all_proteins[,c(5:25)]+1), ncp = lumi_proteins_pc$ncp, scale = F)

res_lumi_proteins <- prcomp(res_lumi_proteins$completeObs, scale. = F, rank = 6)
loadings_lumi <- res_lumi_proteins$rotation
loadings_lumi <- loadings_lumi[order(abs(loadings_lumi[,1]), decreasing = T),]

s_lumi <- summary(res_lumi_proteins)

res_lumi_proteins <- data.frame(pc1 = res_lumi_proteins$x[,1], pc2 = res_lumi_proteins$x[,2],
                                pc3 = res_lumi_proteins$x[,3], pc4 = res_lumi_proteins$x[,4],
                                plate = meta_orig$Plate[match(all_proteins$Sample.ID, meta_orig$ID)], 
                                group  = meta_orig$Group[match(all_proteins$Sample.ID, meta_orig$ID)],
                                sex = meta_orig$sex[match(all_proteins$Sample.ID, meta_orig$ID)],
                                age = meta_orig$age[match(all_proteins$Sample.ID, meta_orig$ID)])

# scale the proteins 
norm_cols <- cbind(all_proteins[,1:2], data.frame(log2(all_proteins[,3:ncol(all_proteins)])+1))
norm_cols <- norm_cols[norm_cols$Phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL"),]

#recalculate PCA
lumi_proteins_pc <- estim_ncpPCA(norm_cols[,c(5:25)], scale = F, method = "Regularized")
lumi_proteins_pc$ncp
res_lumi_proteins <- imputePCA(norm_cols[,c(5:25)], ncp = lumi_proteins_pc$ncp, scale = F)

res_lumi_proteins <- prcomp(res_lumi_proteins$completeObs, scale. = F, rank = 6)
loadings_lumi <- res_lumi_proteins$rotation
loadings_lumi <- loadings_lumi[order(abs(loadings_lumi[,1]), decreasing = T),]

s_lumi <- summary(res_lumi_proteins)

res_lumi_proteins <- data.frame(pc1 = res_lumi_proteins$x[,1], pc2 = res_lumi_proteins$x[,2],
                                pc3 = res_lumi_proteins$x[,3], pc4 = res_lumi_proteins$x[,4],
                                plate = meta_orig$Plate[match(norm_cols$Sample.ID, meta_orig$ID)],
                                group  = meta_orig$Group[match(norm_cols$Sample.ID, meta_orig$ID)],
                                sex = meta_orig$sex[match(norm_cols$Sample.ID, meta_orig$ID)],
                                age = meta_orig$age[match(norm_cols$Sample.ID, meta_orig$ID)])

ggplot(res_lumi_proteins, aes(x = pc1, y = pc2, color = as.factor(plate)))+geom_point()+theme_bw()+stat_ellipse()
ggplot(res_lumi_proteins, aes(x = pc1, y = pc2, color = as.factor(group)))+geom_point()+theme_bw()+stat_ellipse()

# remove plate effects
mod <- model.matrix(~0 + Phenotype, data = norm_cols)[,1]

norm_cols[,c(5:25)] <- t(ComBat(t(norm_cols[,c(5:25)]), 
                                batch = meta_orig$Plate[match(norm_cols$Sample.ID, meta_orig$ID)],
                                mod = mod))
# looks better, now look at PCA for all proteins

# Luminex + ELISA
all_proteins_pc <- estim_ncpPCA(norm_cols[,3:ncol(norm_cols)], scale = F, method = "Regularized")
all_proteins_pc$ncp
res_all_proteins <- imputePCA(norm_cols[,3:ncol(norm_cols)], ncp = all_proteins_pc$ncp, scale = F)

res_all_proteins <- prcomp(res_all_proteins$completeObs, scale. = F, rank = 6)
loadings_all <- res_all_proteins$rotation

s_all <- summary(res_all_proteins)

res_all_proteins <- data.frame(pc1 = res_all_proteins$x[,1], pc2 = res_all_proteins$x[,2],
                               pc3 = res_all_proteins$x[,3], pc4 = res_all_proteins$x[,4],
                               plate = meta_orig$Plate[match(norm_cols$Sample.ID, meta_orig$ID)], 
                               group  = meta_orig$Group[match(norm_cols$Sample.ID, meta_orig$ID)],
                               sex = meta_orig$sex[match(norm_cols$Sample.ID, meta_orig$ID)],
                               age = meta_orig$age[match(norm_cols$Sample.ID, meta_orig$ID)], 
                               patho = meta_orig$Causative.Bacteria[match(norm_cols$Sample.ID, meta_orig$ID)])

ggplot(res_all_proteins[res_all_proteins$group %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL"),], aes(x = pc1, y = pc2, color = group))+
  geom_point()+
  theme_bw()+
  stat_ellipse()
# good separation except for outliers 
# calculate 95% ellipse to remove outliers 
pc_a <- pcaMethods::pca(norm_cols[,3:ncol(norm_cols)], nPcs = 10, scale = 'none')
el <- pcaMethods:::simpleEllipse(pc_a@scores[,1], pc_a@scores[,2], alfa = 0.99, nrow(pc_a@scores))

ggplot(data.frame(pc_a@scores), aes(x = PC1, y = PC2))+geom_point()+
  annotate("path",
           x=el[,1],
           y=el[,2])

table(point.in.polygon(pc_a@scores[,1], pc_a@scores[,2], el[,1], el[,2]))

out_a <- meta_orig[meta_orig$ID %in% rownames(pc_a@scores)[point.in.polygon(pc_a@scores[,1], pc_a@scores[,2], el[,1], el[,2]) == 0],]
norm_cols_out <- norm_cols[!rownames(norm_cols) %in% out_a$ID,]

# Calculate p-values for DB vs DV:
norm_cols_out <- norm_cols_out[norm_cols_out$Phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL"),]
names <- colnames(norm_cols_out)[3:ncol(norm_cols_out)]

wilcox_res <- data.frame(variables = colnames(norm_cols_out)[3:ncol(norm_cols_out)])
wilcox_res$p_value <- NA
tests <- list()

for (i in 3:ncol(norm_cols_out)) {
  df <- data.frame(Group=norm_cols_out$Phenotype, 
                   Protein=norm_cols_out[,i])
  wc <- wilcox.test(df[,2]~df[,1])
  tests[[i-2]] <- wc
  wilcox_res$p_value[i-2] <- wc$p.value
}

wilcox_res$adjusted <- p.adjust(wilcox_res$p_value, method= 'BH', n = nrow(wilcox_res))
wilcox_res$print <- NA
wilcox_res$print[wilcox_res$adjusted > 0.05] <- 'ns'
wilcox_res$print[wilcox_res$adjusted < 0.05] <- '*'
wilcox_res$print[wilcox_res$adjusted <= 0.01] <- '**'
wilcox_res$print[wilcox_res$adjusted <= 0.001] <- '***'
wilcox_res$print[wilcox_res$adjusted <= 0.0001] <- '****'

############################################################################################################
# run FS-PLS 
############################################################################################################
path_fspls <- paste("fspls_lars_multi_meta.R")

response <- class.ind(norm_cols_out$Phenotype)[,1]

colnames(norm_cols_out) <- str_remove_all(colnames(norm_cols_out), "\\-")
luminex_exp <- norm_cols_out[,5:25]
elisa_exp <- norm_cols_out[,4:ncol(norm_cols_out)]

##### run FS-PLS on the Luminex proteins 
fspls_res_luminex <- fspls.iterate(n.iterations = 25, 
                                   seed = 105,
                                   path.fspls = path_fspls, 
                                   expression = luminex_exp, 
                                   response = response, 
                                   p.thresh = 0.001, 
                                   max = 5, 
                                   beam = 10, 
                                   split = T, 
                                   split.ratio = 0.7)

upset_luminex <- fspls_res_luminex$table
upset_luminex <- upset_luminex[,c(1:5, 7:ncol(upset_luminex), 6)]
upset_luminex[is.na(upset_luminex)] <- 0

colSums(upset_luminex)[top_signature_iterate(upset_luminex)$ID]
signature_luminex <- luminex_exp[,top_signature_iterate(upset_luminex)$ID]

# estimate model weights using glm 
signature_luminex$group <- norm_cols_out$Phenotype[match(rownames(signature_luminex), norm_cols_out$Sample.ID)]
signature_luminex$group[signature_luminex$group=="DEFINITE BACTERIAL"] <- 1
signature_luminex$group[signature_luminex$group=="DEFINITE VIRAL"] <- 0
signature_luminex$group <- as.numeric(signature_luminex$group)

luminex_weights <- summary(glm(group ~., data = signature_luminex[,1:6]))$coefficients[2:nrow(summary(glm(group ~., data = signature_luminex))$coefficients),'Estimate']

# weights saved so they can be used in validation 
save(file = 'Luminex_Signature_original_weights.RData', luminex_weights)

signature_luminex$DRS <- as.vector(as.matrix(signature_luminex[,1:length(top_signature_iterate(upset_luminex)$ID)]) %*% 
                                     as.matrix(luminex_weights[match(colnames(signature_luminex)[1:length(top_signature_iterate(upset_luminex)$ID)], names(luminex_weights))]))

luminex_signature_performance <- roc(as.numeric(signature_luminex$group), as.numeric(signature_luminex$DRS), plot = T, print.auc = T, ci = T)

auc(luminex_signature_performance, partial.auc = c(0.9,1.0), partial.auc.focus = 'sensitivity')*100
auc(luminex_signature_performance, partial.auc = c(0.9,1.0), partial.auc.focus = 'specificity')*100

auc(luminex_signature_performance, partial.auc = c(0.95,1.0), partial.auc.focus = 'sensitivity')*100
auc(luminex_signature_performance, partial.auc = c(0.95,1.0), partial.auc.focus = 'specificity')*100

max(luminex_signature_performance$sensitivities[luminex_signature_performance$specificities > 0.95])*100
max(luminex_signature_performance$specificities[luminex_signature_performance$sensitivities > 0.95])*100

coords_lumi <- coords(roc=luminex_signature_performance, x = "all", transpose = FALSE)
coords_lumi <- coords_lumi[coords_lumi$sensitivity >= .90, ]
coords_lumi <- coords_lumi[which(coords_lumi$sensitivity==min(coords_lumi$sensitivity)),]
coords_lumi <- coords_lumi$threshold[coords_lumi$specificity == max(coords_lumi$specificity)]

signature_luminex$class <- 0
signature_luminex$class[signature_luminex$DRS > coords_lumi] <- 1

confmat_lumi <- confusionMatrix(as.factor(signature_luminex$group), 
                                as.factor(signature_luminex$class))

confmat_lumi <- confmat_lumi$byClass

################################################################################################################
#  run FS-PLS including the ELISA proteins 
fspls_res_elisa <- fspls.iterate(n.iterations = 25, 
                                 seed = 105,
                                 path.fspls = path_fspls, 
                                 expression = elisa_exp, 
                                 response = response, 
                                 p.thresh = 0.001, 
                                 max = 5, 
                                 beam = 10, 
                                 split = T, 
                                 split.ratio = 0.7)

upset_elisa <- fspls_res_elisa$table
upset_elisa <- upset_elisa[,c(1:5, 7:ncol(upset_elisa), 6)]
upset_elisa[is.na(upset_elisa)] <- 0

colSums(upset_elisa)[top_signature_iterate(upset_elisa)$ID]

# calculate AUC 
weights_elisa <- summary(glm(group ~.,  data = data.frame(group = response, 
                                                          elisa_exp[,top_signature_iterate(upset_elisa)$ID])))$coefficients[,'Estimate']
weights_elisa <- weights_elisa[2:length(weights_elisa)]

elisa_signature <- data.frame(elisa_exp[,top_signature_iterate(upset_elisa)$ID], 
                              group = response)

elisa_signature[is.na(elisa_signature)] <- 0 

elisa_signature$DRS <- as.vector(as.matrix(elisa_signature[,match(names(weights_elisa), colnames(elisa_signature))]) %*% as.matrix(weights_elisa))

elisa_signature_per <- roc(elisa_signature$group, elisa_signature$DRS, ci = T, auc = T, plot = T, print.auc = T)

auc(elisa_signature_per, partial.auc = c(0.9,1.0), partial.auc.focus = 'sensitivity')*100
auc(elisa_signature_per, partial.auc = c(0.9,1.0), partial.auc.focus = 'specificity')*100

auc(elisa_signature_per, partial.auc = c(0.95,1.0), partial.auc.focus = 'sensitivity')*100
auc(elisa_signature_per, partial.auc = c(0.95,1.0), partial.auc.focus = 'specificity')*100

max(elisa_signature_per$sensitivities[elisa_signature_per$specificities > 0.9])*100
max(elisa_signature_per$specificities[elisa_signature_per$sensitivities > 0.9])*100

coords_elisa <- coords(roc=elisa_signature_per, x = "all", transpose = FALSE)
coords_elisa <- coords_elisa[coords_elisa$sensitivity >= .90, ]
coords_elisa <- coords_elisa[which(coords_elisa$sensitivity==min(coords_elisa$sensitivity)),]
coords_elisa <- coords_elisa$threshold[coords_elisa$specificity == max(coords_elisa$specificity)]

elisa_signature$class <- 0
elisa_signature$class[elisa_signature$DRS > coords_elisa] <- 1

confmat_elisa <- confusionMatrix(as.factor(elisa_signature$group), 
                                 as.factor(elisa_signature$class))

confmat_elisa <- confmat_elisa$byClass

# compare the ROC models 
roc.test(luminex_signature_performance, elisa_signature_per)

# plot the ROCs
roc_plot_lumi <- roc_plot(luminex_signature_performance, title = "A: Luminex signature", 'darkgrey')
roc_plot_elisa <- roc_plot(elisa_signature_per, title = "B: Luminex + ELISA signature", "limegreen")

roc_plots_all <- ggarrange(roc_plot_lumi, roc_plot_elisa, 
                           ncol=2, nrow=1, common.legend = TRUE,legend="bottom")

# plot both on one plot 
ci_lumi <- ci.se(luminex_signature_performance, specificities=seq(0, 1, l=25))
ci_lumi <- data.frame(x = as.numeric(rownames(ci_lumi)),
                      lower = ci_lumi[, 1],
                      upper = ci_lumi[, 3])

ci_elisa <- ci.se(elisa_signature_per, specificities=seq(0, 1, l=25))
ci_elisa <- data.frame(x = as.numeric(rownames(ci_elisa)),
                       lower = ci_elisa[, 1],
                       upper = ci_elisa[, 3])

overlaid_roc <- ggroc(list("Luminex" = luminex_signature_performance, "Luminex + ELISA" = elisa_signature_per))+
  scale_color_manual(values = c("black", "limegreen"))+
  theme_bw()+
  geom_ribbon(
    data = ci_lumi,
    aes(x = x, ymin = lower, ymax = upper),
    alpha = 0.3,
    fill = 'darkgrey',
    inherit.aes = F) +
  geom_ribbon(
    data = ci_elisa,
    aes(x = x, ymin = lower, ymax = upper),
    alpha = 0.2,
    fill = 'limegreen',
    inherit.aes = F)+
  labs(x = '1-Specificity', y = "Sensitivity", color = 'Signature')


####### generate confusion matrices 
signature_luminex$group
confmat_lumi <- confusionMatrix(reference = as.factor(signature_luminex$group), 
                                data = as.factor(signature_luminex$class))

confmat_lumi$table

confmat_elisa <- confusionMatrix(reference = as.factor(elisa_signature$group),
                                 data = as.factor(elisa_signature$class))

confmat_elisa$table 

###################### differential expression analysis to determine whether there are proteins that work well with CRP is lower 
meta_orig$CRP <- final_phenos$CRP_FIRST_RESEARCH_BLOOD[match(meta_orig$ID, final_phenos$UNIQUE_PATIENT_ID)]
meta_orig$CRP <- str_remove_all(meta_orig$CRP, '<')
meta_orig$CRP <- str_remove_all(meta_orig$CRP, '>')
meta_orig$CRP <- str_replace_all(meta_orig$CRP, '\\,', '\\.')
meta_orig$CRP <- as.numeric(meta_orig$CRP)

table(meta_orig$Group[is.na(meta_orig$CRP)])

meta_orig$CRP[is.na(meta_orig$CRP)] <- final_phenos$CRP_BASELINE[match(meta_orig$ID[is.na(meta_orig$CRP)], final_phenos$UNIQUE_PATIENT_ID)]
meta_orig$CRP[is.na(meta_orig$CRP)] <- final_phenos$CRP_SECOND_RESEARCH_BLOOD[match(meta_orig$ID[is.na(meta_orig$CRP)], final_phenos$UNIQUE_PATIENT_ID)]
meta_orig$CRP[is.na(meta_orig$CRP)] <- final_phenos$CRP_THIRD_RESEARCH_BLOOD[match(meta_orig$ID[is.na(meta_orig$CRP)], final_phenos$UNIQUE_PATIENT_ID)]
meta_orig$CRP <- str_remove_all(meta_orig$CRP, '<')
meta_orig$CRP <- str_remove_all(meta_orig$CRP, '>')
meta_orig$CRP <- str_replace_all(meta_orig$CRP, '\\,', '\\.')

meta_orig$CRP <- as.numeric(meta_orig$CRP)

meta_low_crp <- subset(meta_orig, CRP <= 60 & Group == "DEFINITE BACTERIAL" |
                         Group == 'DEFINITE VIRAL')
meta_low_crp <- meta_low_crp[meta_low_crp$ID %in% rownames(luminex_exp),]

# 31 x DB, 113 x DV
abundances_low_crp <- luminex_exp[rownames(luminex_exp) %in% meta_low_crp$ID,]
abundances_low_crp <- abundances_low_crp[match(meta_low_crp$ID, rownames(abundances_low_crp)),]

mm <- model.matrix(~0 + as.factor(Group) + as.factor(sex) + age, data = meta_low_crp)
colnames(mm) <- c("DB", 'DV', 'M', 'age')

# remove proteins in signature already 
# run Limma
limma_low_crp <- limma.fun(data = abundances_low_crp, 
                           p.thresh = 0.05, 
                           comparisons = c("DB", "DV"), 
                           start.col = 1,
                           end.col = ncol(abundances_low_crp),
                           n.features = ncol(abundances_low_crp),
                           model.mat = mm)