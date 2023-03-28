###################### 28th March 2023
### Validation Phase of Protein Signature Study

library(lubridate)
library(ggpubr)
library(viridis)
library(tidyverse)
library(clipr)
load.packages()
library(caret)

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

recalculate_auc <- function(expression, directions, n_features, group){
  rocs <- list()
  for(i in 1:n_features){
    if(i == 1){
      drs <- expression[,i]
      roc = roc(response = group, predictor = drs)
      rocs[[i]] <- roc
    }
    else{
      drs <- calculate_drs(as.data.frame(expression[,1:i]), directions[1:i])
      roc <- roc(response = group, predictor = drs$DRS, ci = T)
      rocs[[i]] <- roc
    }
  }
  return(rocs)
}

calculate_drs <- function(expression, weights){
  drs <- data.frame(DRS = as.matrix(expression) %*% as.matrix(weights))
  return(drs)
}


#############################################################################
# read in protein measurements 
protein_res <- read.csv(file = 'data/protein_measurements.csv', header = T)

# viusalise 
protein_res_b <- pivot_longer(protein_res, cols = 4:ncol(protein_res))

ggplot(protein_res_b,aes(x = name, y = value))+
  theme_bw()+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# log2 transform data
protein_res_out <- data.frame(log2(protein_res_out+1))

########################################################################
# read in model weights
load(file = 'Luminex_Signature_original_weights.RData', verbose = T)

signature_5 <- protein_res_out[,names(luminex_weights)]

# calculate disease risk scores 
signature_5[is.na(signature_5)] <- 0
# weighted DRS using original weights
signature_5$drs <- as.vector(as.matrix(signature_5) %*% luminex_weights)

signature_5$phenotype <- patients_match$FINAL_PHENOTYPE_DIAGNOSIS[match(rownames(signature_5), patients_match$UNIQUE_PATIENT_ID)]

signature_5$db_dv <- 0
signature_5$db_dv[signature_5$phenotype == 'DEFINITE BACTERIAL'] <- 1

# estimate performance using original model weights 
original_weights_per <- roc(signature_5$db_dv[signature_5$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], 
                            signature_5$drs[signature_5$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], plot = T, print.auc=T, ci = T)

auc(original_weights_per, partial.auc = c(0.9,1.0), partial.auc.focus = 'sensitivity')*100
auc(original_weights_per, partial.auc = c(0.9,1.0), partial.auc.focus = 'specificity')*100

auc(original_weights_per, partial.auc = c(0.95,1.0), partial.auc.focus = 'sensitivity')*100
auc(original_weights_per, partial.auc = c(0.95,1.0), partial.auc.focus = 'specificity')*100

## calculate sensitivity and specificity using Youden's 
coords_original <- coords(roc=original_weights_per, input = 'threshold', transpose = FALSE, x = 'best')

coords_original <- coords_original$threshold[coords_original$specificity == max(coords_original$specificity)]

signature_5$original_class <- 0
signature_5$original_class[signature_5$drs > coords_original] <- 1

confmat_original <- confusionMatrix(reference = as.factor(signature_5$db_dv[signature_5$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")]), 
                                    data = as.factor(signature_5$original_class[signature_5$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")]))

confmat_original$byClass

# in order to make predictions, set a limit on sensitivity 
coords_original_predictions <- coords(roc=original_weights_per, x = "all", transpose = FALSE)
coords_original_predictions <- coords_original_predictions[coords_original_predictions$sensitivity >= .90, ]
thresh_original_sens <- coords_original_predictions[which(coords_original_predictions$sensitivity==min(coords_original_predictions$sensitivity)),]
thresh_original_sens <- thresh_original_sens$threshold[thresh_original_sens$specificity == max(thresh_original_sens$specificity)]

signature_5$original_class_sens <- 0
signature_5$original_class_sens[signature_5$drs > thresh_original_sens] <- 1

confmat_original_sens <- confusionMatrix(reference = as.factor(signature_5$db_dv[signature_5$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                                              'DEFINITE VIRAL')]), 
                                         data = as.factor(signature_5$original_class_sens[signature_5$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                                                       'DEFINITE VIRAL')]))

confmat_original_sens$table

#### use simple DRS
signature_5$simple_drs <- (rowSums(signature_5[,names(luminex_weights)[luminex_weights > 0]])-rowSums(signature_5[,names(luminex_weights)[luminex_weights < 0]]))

simple_per <- roc(signature_5$db_dv[signature_5$phenotype %in% c("DEFINITE BACTERIAL", 'DEFINITE VIRAL')], 
                  signature_5$simple_drs[signature_5$phenotype %in% c("DEFINITE BACTERIAL", 'DEFINITE VIRAL')], plot= T, print.auc = T, ci = T)

auc(simple_per, partial.auc = c(0.9,1.0), partial.auc.focus = 'sensitivity')*100
auc(simple_per, partial.auc = c(0.9,1.0), partial.auc.focus = 'specificity')*100

auc(simple_per, partial.auc = c(0.95,1.0), partial.auc.focus = 'sensitivity')*100
auc(simple_per, partial.auc = c(0.95,1.0), partial.auc.focus = 'specificity')*100

## calculate sensitivity and specificity using Youden's 
coords_simple <- coords(roc=simple_per, input = 'threshold', transpose = FALSE, x = 'best')

coords_simple <- coords_simple$threshold[coords_simple$specificity == max(coords_simple$specificity)]

signature_5$simple_class <- 0
signature_5$simple_class[signature_5$simple_drs > coords_simple] <- 1

confmat_simple <- confusionMatrix(reference = as.factor(signature_5$db_dv[signature_5$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")]), 
                                  data = as.factor(signature_5$simple_class[signature_5$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")]))

confmat_simple$byClass

# reclssify samples using sensitivity threshold
coords_simple_pred <- coords(roc=simple_per, x = "all", transpose = FALSE)
coords_simple_pred <- coords_simple_pred[coords_simple_pred$sensitivity >= .90, ]
thresh_simple_sens <- coords_simple_pred[which(coords_simple_pred$sensitivity==min(coords_simple_pred$sensitivity)),]
thresh_simple_sens <- thresh_simple_sens$threshold[thresh_simple_sens$specificity == max(thresh_simple_sens$specificity)]

signature_5$simple_class_sens <- 0
signature_5$simple_class_sens[signature_5$simple_drs > thresh_simple_sens] <- 1

confmat_simple_sens <- confusionMatrix(reference = as.factor(signature_5$db_dv[signature_5$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                                            'DEFINITE VIRAL')]), 
                                       data = as.factor(signature_5$simple_class_sens[signature_5$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                                                   'DEFINITE VIRAL')]))

confmat_simple_sens$table

# estimate DRS with new weights 
# re-estimate weights 
new_weights <- summary(glm(db_dv ~., data = signature_5[signature_5$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL"),
                                                        c(1:5, 8)]))$coefficients[,'Estimate']
new_weights <- new_weights[2:length(new_weights)]

signature_5$retrained_drs <- as.vector(as.matrix(signature_5[,match(names(new_weights), colnames(signature_5))]) %*% as.matrix(new_weights))

retrained_per <- roc(signature_5$db_dv[signature_5$phenotype %in% c("DEFINITE BACTERIAL", 'DEFINITE VIRAL')], 
                     signature_5$retrained_drs[signature_5$phenotype %in% c("DEFINITE BACTERIAL", 'DEFINITE VIRAL')], plot= T, print.auc = T, ci = T)

auc(retrained_per, partial.auc = c(0.9,1.0), partial.auc.focus = 'sensitivity')*100
auc(retrained_per, partial.auc = c(0.9,1.0), partial.auc.focus = 'specificity')*100

auc(retrained_per, partial.auc = c(0.95,1.0), partial.auc.focus = 'sensitivity')*100
auc(retrained_per, partial.auc = c(0.95,1.0), partial.auc.focus = 'specificity')*100


## calculate sensitivity and specificity using Youden's 
coords_retrained <- coords(roc=retrained_per, input = 'threshold', transpose = FALSE, x = 'best')

coords_retrained <- coords_retrained$threshold[coords_retrained$specificity == max(coords_retrained$specificity)]

signature_5$retrained_class <- 0
signature_5$retrained_class[signature_5$retrained_drs > coords_retrained] <- 1

confmat_retrained <- confusionMatrix(reference = as.factor(signature_5$db_dv[signature_5$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")]), 
                                     data = as.factor(signature_5$retrained_class[signature_5$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")]))

confmat_retrained$byClass

# reclassify samples using sensitivity threshold
coords_retrained_pred <- coords(roc=retrained_per, x = "all", transpose = FALSE)
coords_retrained_pred <- coords_retrained_pred[coords_retrained_pred$sensitivity >= .90, ]
thresh_retrained_sens <- coords_retrained_pred[which(coords_retrained_pred$sensitivity==min(coords_retrained_pred$sensitivity)),]
thresh_retrained_sens <- thresh_retrained_sens$threshold[thresh_retrained_sens$specificity == max(thresh_retrained_sens$specificity)]

signature_5$ret_class_sens <- 0
signature_5$ret_class_sens[signature_5$retrained_drs > thresh_retrained_sens] <- 1

confmat_retrained_sens <- confusionMatrix(reference = as.factor(signature_5$db_dv[signature_5$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                                               'DEFINITE VIRAL')]), 
                                          data = as.factor(signature_5$ret_class_sens[signature_5$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                                                   'DEFINITE VIRAL')]))

confmat_retrained_sens$table

############ add on final protein 
additional_proteins <- protein_res_out[,!colnames(protein_res_out) %in% colnames(signature_5)]
additional_proteins[is.na(additional_proteins)] <- 0
additional_proteins[additional_proteins==-Inf] <- 0

signature_6 <- cbind(signature_5[,c(1:5,7,8)], additional_proteins[,c('final_protein')])

# direction of final protein is negative
weights_6 <- c(luminex_weights, -0.706)

signature_6$simple_drs <- (rowSums(signature_6[,names(weights_6)[weights_6 > 0]])-rowSums(signature_6[,names(weights_6)[weights_6 < 0]]))

six_simple_drs <- roc(signature_6$db_dv[signature_6$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], 
                      signature_6$simple_drs[signature_6$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], ci = T)

roc.test(six_simple_drs, simple_per)

# retrain weights 
weights_6_retrained <- summary(glm(db_dv ~., data = signature_6[signature_6$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL"), c(1:5, 7:8)]))$coefficients
weights_6_retrained <- weights_6_retrained[2:nrow(weights_6_retrained),'Estimate']

signature_6$retrained_drs <- as.vector(as.matrix(signature_6[,match(names(weights_6_retrained), colnames(signature_6))]) %*% as.matrix(weights_6_retrained))

retrained_6 <- roc(signature_6$db_dv[signature_6$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], 
                   signature_6$retrained_drs[signature_6$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], ci = T)

roc.test(retrained_6, retrained_per)

roc_plot_6 <- roc_plot(six_simple_drs, "D: 6-protein signature", "pink")
roc_plot_simple_5 <- roc_plot(simple_per, "5-protein signature", 'darkgrey')

##### plot all together 
# roc_plot_original
# roc_plot_simple_5
# roc_plot_6
# roc_plot_retrained

ci_original <- ci.se(original_weights_per, specificities=seq(0, 1, l=25))
ci_original <- data.frame(x = as.numeric(rownames(ci_original)),
                          lower = ci_original[, 1],
                          upper = ci_original[, 3])

ci_simple <- ci.se(simple_per, specificities=seq(0, 1, l=25))
ci_simple <- data.frame(x = as.numeric(rownames(ci_simple)),
                        lower = ci_simple[, 1],
                        upper = ci_simple[, 3])

ci_retrained <- ci.se(retrained_per, specificities=seq(0, 1, l=25))
ci_retrained <- data.frame(x = as.numeric(rownames(ci_retrained)),
                           lower = ci_retrained[, 1],
                           upper = ci_retrained[, 3])

ci_6 <- ci.se(six_simple_drs, specificities=seq(0, 1, l=25))
ci_6 <- data.frame(x = as.numeric(rownames(ci_6)),
                   lower = ci_6[, 1],
                   upper = ci_6[, 3])

viridis_pal <- c('#3E92CC', '#29BF12', '#FF715B', '#211A1D')

overlaid_roc_all <- ggroc(list("Original weights (5 proteins)" = original_weights_per,
                               "Simple DRS (5 proteins)" = simple_per, 
                               "Retrained weights (5 proteins)" = retrained_per, 
                               "Simple DRS (6 proteins)" = six_simple_drs))+
  theme_bw()+
  scale_color_manual(values = viridis_pal)+
  geom_ribbon(
    data = ci_original,
    aes(x = x, ymin = lower, ymax = upper),
    alpha = 0.2,
    fill = viridis_pal[1],
    inherit.aes = F)+
  geom_ribbon(
    data = ci_simple,
    aes(x = x, ymin = lower, ymax = upper),
    alpha = 0.3,
    fill = viridis_pal[2],
    inherit.aes = F) +
  geom_ribbon(
    data = ci_retrained,
    aes(x = x, ymin = lower, ymax = upper),
    alpha = 0.2,
    fill = viridis_pal[3],
    inherit.aes = F)+
  geom_ribbon(
    data = ci_6,
    aes(x = x, ymin = lower, ymax = upper),
    alpha = 0.2,
    fill = viridis_pal[4],
    inherit.aes = F)+
  labs(x = '1-Specificity', y = "Sensitivity", color = 'Signature')

png(file = 'ROC2_overlaid.png', 
    height = 5, width = 7.5, units = 'in', res = 1500)
overlaid_roc_all
dev.off()

roc_plot_original <- roc_plot(original_weights_per, title = "A: Original weights", 'darkgrey')
roc_plot_retrained <- roc_plot(retrained_per, title = "B: Retrained weights", "darkgrey")
roc_plot_simple <- roc_plot(simple_per, title = "C: Simple DRS", "darkgrey")

roc_plots_all <- ggarrange(roc_plot_original, roc_plot_retrained, roc_plot_simple,roc_plot_6,
                           ncol=2, nrow=2, common.legend = TRUE,legend="bottom")

png(file = 'ROC2_all.png', 
    height = 10, width = 10, units = 'in', res = 1500)
roc_plots_all
dev.off()

######## test performance of 6 protein signature in other groups 

# use new weights to estimate DRS for other groups 
signature_6 <- signature_6[signature_6$phenotype %in% c( 'DEFINITE BACTERIAL', 'Non-Sterile DB',
                                                         'PROBABLE BACTERIAL', 'BACTERIAL SYNDROME', 
                                                         'DEFINITE VIRAL', 'PROBABLE VIRAL', 'VIRAL SYNDROME'), ]

signature_6$phenotype <- as.character(signature_6$phenotype)
signature_6$phenotype[signature_6$phenotype=="Non-Sterile DB"] <- 'NON-STERILE DB'

signature_6$phenotype <- fct_relevel(as.factor(signature_6$phenotype), 'DEFINITE BACTERIAL', 'NON-STERILE DB',
                                     'PROBABLE BACTERIAL', 'BACTERIAL SYNDROME', 
                                     'VIRAL SYNDROME', 'PROBABLE VIRAL', 
                                     'DEFINITE VIRAL')

signature_6$phenotype <- droplevels(signature_6$phenotype)

coords <- coords(roc=retrained_6, x = "all", transpose = FALSE)
coords <- coords[coords$sensitivity >= .90, ]
thresh <- coords[which(coords$sensitivity==min(coords$sensitivity)),]
thresh <- thresh$threshold[thresh$specificity == max(thresh$specificity)]


pdf(file='6_proteins_rest_retrained.pdf', height = 7, width = 7)
ggplot(signature_6, aes(x = phenotype, 
                        y = retrained_drs, 
                        color = phenotype))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  theme_bw()+
  geom_hline(yintercept = thresh, linetype = 'dashed')+
  labs(x = '', y = "Retrained Disease Risk Score")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_color_manual(values = c('#a50026','#d73027','#f46d43','#fdae61',
                                '#74add1','#4575b4', '#313695'))+
  theme(legend.position = "none") 
dev.off() 

signature_6$classification <- 0
signature_6$classification[signature_6$retrained_drs > thresh] <- 1

table(signature_6$classification, signature_6$phenotype)

####### oved et al 

############
# how does TRAIL and IP10 perform
ip10_trail <-   protein_res_out[,c('TRAIL..pg.ml.', 'CXCL10.IP.10..pg.ml.')]
ip10_trail <- ip10_trail[rownames(ip10_trail) %in% rownames(signature_5),]
ip10_trail$phenotype <- signature_5$phenotype[match(rownames(ip10_trail), rownames(signature_5))]
ip10_trail$db_dv <- signature_5$db_dv[match(rownames(ip10_trail), rownames(signature_5))]

ip10_trail[is.na(ip10_trail)] <- 0

# calculate score by adding up values
ip10_trail$drs <- apply(ip10_trail, 1, function(x){sum(as.numeric(x[1:2]))})

ip_trail_per <- roc(ip10_trail$db_dv[ip10_trail$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], 
                    ip10_trail$drs[ip10_trail$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], ci = T)

# youdens for optimal sens and spec 
coords_ip10_trail <- coords(roc=ip_trail_per, input = 'threshold', transpose = FALSE, x = 'best')
coords_ip10_trail <- coords_ip10_trail$threshold[coords_ip10_trail$specificity == max(coords_ip10_trail$specificity)]

ip10_trail$class <- 0
ip10_trail$class[ip10_trail$drs > coords_trail_ip] <- 1

confmat_trail_ip <- confusionMatrix(as.factor(ip10_trail$db_dv[ip10_trail$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                           'DEFINITE VIRAL')]), 
                                    as.factor(ip10_trail$class[ip10_trail$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                           'DEFINITE VIRAL')]))
confmat_trail_ip$byClass

# add in CRP to TRAIL and IP10 
crp_ip10_trail <-   protein_res_out[,c('CRP.mg.L', 'TRAIL..pg.ml.', 'CXCL10.IP.10..pg.ml.')]
crp_ip10_trail <- crp_ip10_trail[rownames(crp_ip10_trail) %in% rownames(signature_5),]
crp_ip10_trail$phenotype <- signature_5$phenotype[match(rownames(crp_ip10_trail), rownames(signature_5))]
crp_ip10_trail$db_dv <- signature_5$db_dv[match(rownames(crp_ip10_trail), rownames(signature_5))]

crp_ip10_trail[is.na(crp_ip10_trail)] <- 0
#crp_ip10_trail[is.na(crp_ip10_trail)] <- 0

crp_ip10_trail$drs <- apply(crp_ip10_trail, 1, function(x){sum(
  as.numeric(x[1]) - sum(as.numeric(x[2:3]))
)
})

crp_ip_trail_per <- roc(crp_ip10_trail$db_dv[crp_ip10_trail$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], 
                        crp_ip10_trail$drs[crp_ip10_trail$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], ci = T)

coords_crp_ip <- coords(roc=crp_ip_trail_per, input = 'threshold', transpose = FALSE, x = 'best')
coords_crp_ip <- coords_crp_ip$threshold[coords_crp_ip$specificity == max(coords_crp_ip$specificity)]

crp_ip10_trail$class <- 0
crp_ip10_trail$class[crp_ip10_trail$drs > coords_crp_ip] <- 1

confmat_crp_ip <- confusionMatrix(as.factor(crp_ip10_trail$db_dv[crp_ip10_trail$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                                 'DEFINITE VIRAL')]), 
                                  as.factor(crp_ip10_trail$class[crp_ip10_trail$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                                 'DEFINITE VIRAL')]))

confmat_crp_ip$byClass

######## calculate CRP+IP10+TRAIL using their model 
library(nnet)
crp_ip10_trail <- crp_ip10_trail[crp_ip10_trail$phenotype %in% c("DEFINITE VIRAL", 'DEFINITE BACTERIAL'),]

test <- multinom(as.factor(class) ~ CRP.mg.L + TRAIL..pg.ml. + CXCL10.IP.10..pg.ml., data = crp_ip10_trail)

#crp_ip10_trail$multinom_residuals <- test$residuals
crp_ip10_trail$fitted_values <- test$fitted.values[,1]

ggplot(crp_ip10_trail, aes(x = as.factor(class), y = fitted_values, color = class))+
  geom_boxplot()+
  theme_bw()

roc(crp_ip10_trail$db_dv, crp_ip10_trail$fitted_values, ci = T)

#add CRP into the 6 proteins 
signature_6_crp <- cbind(signature_6[,c(1:5,8)], 
                         CRP = protein_res_out$CRP.mg.L[match(rownames(signature_6), rownames(protein_res_out))], 
                         signature_6[,c(6:7)])

signature_6_crp[is.na(signature_6_crp)] <- 0

weights_6_crp <- c(weights_6, 0.5)
names(weights_6_crp)[length(weights_6_crp)] <- 'CRP'

signature_6_crp$drs <- (rowSums(signature_6_crp[,names(weights_6_crp)[weights_6_crp > 0]])-rowSums(signature_6_crp[,names(weights_6_crp)[weights_6_crp < 0]]))

crp_6_per <- roc(signature_6_crp$db_dv[signature_6_crp$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], 
                 signature_6_crp$drs[signature_6_crp$phenotype %in% c("DEFINITE BACTERIAL", "DEFINITE VIRAL")], ci = T)


# sens and spec for simple DRS
coords_crp_6 <- coords(roc=crp_6_per, input = 'threshold', transpose = FALSE, x = 'best')
coords_crp_6 <- coords_crp_6$threshold[coords_crp_6$specificity == max(coords_crp_6$specificity)]

signature_6_crp$crp_class <- 0
signature_6_crp$crp_class[signature_6_crp$drs > coords_crp_6] <- 1

confmat_crp_6 <- confusionMatrix(as.factor(signature_6_crp$db_dv[signature_6_crp$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                                  'DEFINITE VIRAL')]), 
                                 as.factor(signature_6_crp$crp_class[signature_6_crp$phenotype %in% c("DEFINITE BACTERIAL", 
                                                                                                      'DEFINITE VIRAL')]))

confmat_crp_6$byClass

# with CRP
roc.test(crp_ip_trail_per, crp_6_per)

# without CRP 
roc.test(ip_trail_per, six_simple_drs)


#### plot the ROC curves
plot_crp_ip10_trail_simple <- roc_plot(crp_ip_trail_per, title = "A: CRP+IP10+TRAIL", 'darkgrey')
plot_ip10_trail_simple <- roc_plot(ip_trail_per, title = "B: IP10+TRAIL", 'darkgrey')
plot_6_crp_simple <- roc_plot(crp_6_per, title = "C: 6-protein signature + CRP", 'pink')
plot_6_simple <- roc_plot(six_simple_drs, "D: 6-protein signature", "pink")

png(file = 'Oved_et_al_comparison.png', 
    height = 10, width = 10, res = 1500, units = 'in')
grid.arrange(plot_crp_ip10_trail_simple, plot_ip10_trail_simple,
             plot_6_crp_simple, plot_6_simple, ncol = 2)
dev.off()


########### test 6 protein signature using simple DRS by age 
final_phenos$age <- NA

# read in complete metadata 
perform <- read.csv(file = 'PERFORM_BIVA_ED_EXPORT_SEPT22_2021.csv', header = T)
patients_match <- patients_match[patients_match$UNIQUE_PATIENT_ID %in% rownames(protein_res_out),]
perform <- perform[perform$UNIQUE_PATIENT_ID %in% patients_match$UNIQUE_PATIENT_ID,]

perform$age <- interval(as.Date(dmy(perform$DATE_OF_BIRTH)), 
                        as.Date(dmy_hm(perform$DATE_TIME_FIRST_RESEARCH_BLOOD))) %/% months(1)

signature_6$age <- perform$age[match(rownames(signature_6), 
                                     perform$UNIQUE_PATIENT_ID)]

signature_6_db_dv <- signature_6[signature_6$phenotype  %in% c('DEFINITE BACTERIAL', 'DEFINITE VIRAL'), ]

hist(signature_6_db_dv$age)
range(signature_6_db_dv$age, na.rm = T)

signature_6_db_dv$age_group <- NA
signature_6_db_dv$age_group[signature_6_db_dv$age > 120] <- "10-18"
signature_6_db_dv$age_group[signature_6_db_dv$age < 120] <- "5-10"
signature_6_db_dv$age_group[signature_6_db_dv$age < 60] <- "2-5"
signature_6_db_dv$age_group[signature_6_db_dv$age < 24] <- "Under 2"

# Under 2 performance
roc(signature_6_db_dv$db_dv[signature_6_db_dv$age_group == "Under 2"], 
    signature_6_db_dv$simple_drs[signature_6_db_dv$age_group == "Under 2"], ci = T)

roc.test(six_simple_drs, roc(signature_6_db_dv$db_dv[signature_6_db_dv$age_group == "Under 2"], 
                             signature_6_db_dv$simple_drs[signature_6_db_dv$age_group == "Under 2"], ci = T))

# 2-5 performance
roc(signature_6_db_dv$db_dv[signature_6_db_dv$age_group == "2-5"], 
    signature_6_db_dv$simple_drs[signature_6_db_dv$age_group == "2-5"], ci = T)

roc.test(six_simple_drs, roc(signature_6_db_dv$db_dv[signature_6_db_dv$age_group == "2-5"], 
                             signature_6_db_dv$simple_drs[signature_6_db_dv$age_group == "2-5"], ci = T))

# 5-10 performance
roc(signature_6_db_dv$db_dv[signature_6_db_dv$age_group == "5-10"], 
    signature_6_db_dv$simple_drs[signature_6_db_dv$age_group == "5-10"], ci = T)

roc.test(six_simple_drs, roc(signature_6_db_dv$db_dv[signature_6_db_dv$age_group == "5-10"], 
                             signature_6_db_dv$simple_drs[signature_6_db_dv$age_group == "5-10"], ci = T))

# 10-18 performance
roc(signature_6_db_dv$db_dv[signature_6_db_dv$age_group == "10-18"], 
    signature_6_db_dv$simple_drs[signature_6_db_dv$age_group == "10-18"], ci = T)

roc.test(six_simple_drs, roc(signature_6_db_dv$db_dv[signature_6_db_dv$age_group == "10-18"], 
                             signature_6_db_dv$simple_drs[signature_6_db_dv$age_group == "10-18"], ci = T))