############################################################################################################################################################
# ------------- SOMASCAN 
# ------------- This script is for QC, normalisation, metafile data addition, ID matching for SomaScan
############################################################################################################################################################

source("https://bioconductor.org/biocLite.R")

library(pcaMethods)
library(COCONUT)
library(readat)
library(dplyr)
require(gdata)
library(stringr)
library(viridis)

############################################################################################################################################################
# ------------- Read in the files 
############################################################################################################################################################

# keep only commands ensure that normalisation is not done automatically 
rawfile <- readAdat(file = "raw.adat", keepOnlyPasses = FALSE, keepOnlySamples = FALSE)

############################################################################################################################################################
# ------------- NORMALISATION OF DATA
# ------------- Hybridisation normalisation has already been done
# ------------- calculate the place scalar scale factor 
############################################################################################################################################################

# split DF into plate sub-sections
plate1 <- rawfile[rawfile$PlateId == "Set 001",]
plate2 <- rawfile[rawfile$PlateId == "Set 002",]
plate3 <- rawfile[rawfile$PlateId == "Set 003",]

# get sequence data
seqData <- getSequenceData(normalisedfile)
seqData$SeqId <- paste("SeqId", seqData$SeqId, sep = ".")

# apply median scales to all expression measurements in each plate data frame 
# the median is obtained from the calibration scores 
# calibration samples are added, median taken from protein measuremtns across these calibration scores 
# ratio of median to plate scale and then median taken of ratios per plate.
plate1_scaled <- as.data.frame(apply(plate1[,22:1326], c(1,2), function(x){
  y <- x*1.24
  y <- round(y,1)
  return(y)
}))

plate2_scaled <- as.data.frame(apply(plate2[,22:1326], c(1,2), function(x){
  y <- x*0.802  
  y <- round(y,1)
  return(y)
}))

plate3_scaled <- as.data.frame(apply(plate3[,22:1326], c(1,2), function(x){
  y <- x*1.005
  y <- round(y,1)
}))

# rejoin dataframes and add metadata
merged <- rbind(plate1_scaled, plate2_scaled, plate3_scaled)
preExp <- rbind(plate1[,1:21], plate2[,1:21], plate3[,1:21])
plateScaleFile <- cbind(preExp, merged)
plateScaleFile <- plateScaleFile[order(plateScaleFile$SampleDescription),]

############################################################################################################################################################
# ------------- NORMALISATION OF DATA
# ------------- Hybridisation normalisation and plate scaling has already been done
# ------------- Median normalisation now necessary 
############################################################################################################################################################

# proteins have been diluted either 0.05, 1, or 40%
# must be normalised according to their dilution factor
hybFile <- plateScaleFile

# split sequence data frame according to dilutions 
proteins1 <- seqData[seqData$Dilution == "1", ]
proteins05 <- seqData[seqData$Dilution == "0.005",]
proteins40 <- seqData[seqData$Dilution == "40", ]

# add onto the hybridisation file the columns corresponding to the scales for each dilution 
# first the outlier sample (identified by SomaLogic) will have to be removed 
normalisedFileOutlier <- normalisedfile[(normalisedfile$SampleDescription %in% hybFile$SampleDescription), ]
# then add columns onto hybFile 
hybFile <- cbind(normalisedFileOutlier$NormScale_0_005, normalisedFileOutlier$NormScale_1, normalisedFileOutlier$NormScale_40, hybFile)
colnames(hybFile)[1] <- "NormScale_0_005"
colnames(hybFile)[2] <- "NormScale_1"
colnames(hybFile)[3] <- "NormScale_40"

# medianFunction is a function for median normalisation that takes the dilution scale as a variable so can be used for all dilution
# the sequence data file has been split so will be fed into this function in 3 parts and the median normalistion performed for proteins in each dilution class seperately
# median factor is specific per sample (choice of 3) and per protein because each protein is one of 3 dilutions 
medianFunction <- function(expSet, proteinSet, dilution){
  id <- proteinSet[[1]]
  # pull out the normalisation scale factor for all samples 
  normScale <- as.data.frame(expSet[,..dilution])
  proteinValues <- as.data.frame(expSet[,..id])
  # create sub dataframe of the normalisation scale for the sample and the expression values 
  subdf <- cbind(normScale, proteinValues)
  medScaled <- apply(subdf, 1, function(x){
    # multiple the expression by the scale factor 
    normalised <- x[1]*x[2]
    round(normalised, 1)
  })
  medScaled <- as.data.frame(medScaled)
  colnames(medScaled)[1] <- id
  return(medScaled)
}

# call function for median normalisation for proteins in 1% dilution 
median_1 <- apply(proteins1, 1, function(x){
  df <- medianFunction(hybFile, x, "NormScale_1")
})
median_1 <- as.data.frame(median_1)

# call function for median normalisation for proteins in 0.005% dilution 
median_05 <- apply(proteins05, 1, function(x){
  df <- medianFunction(hybFile, x, "NormScale_0_005")
})
median_05 <- as.data.frame(median_05)

# call function for median normalisation for proteins in 40% dilution 
median_40 <- apply(proteins40, 1, function(x){
  df <- medianFunction(hybFile, x, "NormScale_40")
})
median_40 <- as.data.frame(median_40)

# bind the data frames together after having median normalised and add on the initial columns
# hybMedFile has now gone through: HYBRIDISATION NORMALISATION AND MEDIAN NORMALISATION 
hybMedFile <- cbind(hybFile[,1:24], median_1, median_05, median_40)

############################################################################################################################################################
# ------------- NORMALISATION OF DATA
# ------------- Hybridisation normalisation, plate scaling and median normalisation has already been done
# ------------- Plate calibration  
############################################################################################################################################################

# each sample is on one of 3 plates
# each protein has a scale factor for EACH plate 
# therefore depending on which plate the sample is on, the protein expression value is multiplied by one of 3 scale factors 

# split expression matrix according to plate 
hybMedFile <- hybMedFile[order(hybMedFile$PlateId),] 
plateOne <- hybMedFile[(hybMedFile$PlateId == "Set 001"), ]
plateTwo <- hybMedFile[(hybMedFile$PlateId == "Set 002"), ]
plateThr <- hybMedFile[(hybMedFile$PlateId == "Set 003"), ]

# scalingColumnFunction is a function that pulls out the calibration scale factor according to the calibration set being tested i.e. the plate 
scalingColumnFunction <- function(columnX, seqData, calSet){
  seqId <- colnames(columnX)
  # this takes out the scale factor found in the sequence data set 
  scaleFac <- as.numeric(seqData[,..calSet][grep(seqId, seqData$SeqId)])
  colScaled <- lapply(columnX, function(x){
    # multiple the expression value by the scale factor 
    scaled <- x*scaleFac
    round(scaled, 1)
    return(scaled)
  })
  colScaled <- data.frame(matrix(unlist(colScaled)))
}

# apply the function to every column in the data frames split according to plate 
# first plate 1
plateOneCalibration <- lapply(colnames(plateOne)[25:1329], function(x){
  column <- plateOne[,..x]
  plateScaled <- scalingColumnFunction(column, seqData, "Cal_Set_001")
})

# transform into data frame matching orignal expression set 
plateOneCalibration <- as.data.frame(plateOneCalibration)
plateOneCalibration <- round(plateOneCalibration, 1)

# next plate 2
plateTwoCalibration <- lapply(colnames(plateTwo)[25:1329], function(x){
  column <- plateTwo[,..x]
  plateScaled <- scalingColumnFunction(column, seqData, "Cal_Set_002")
})

# transform into data frame matching orignal expression set 
plateTwoCalibration <- as.data.frame(plateTwoCalibration)
plateTwoCalibration <- round(plateTwoCalibration, 1)

# last plate 3
plateThrCalibration <- lapply(colnames(plateThr)[25:1329], function(x){
  column <- plateThr[,..x]
  plateScaled <- scalingColumnFunction(column, seqData, "Cal_Set_003")
})
# transform into data frame matching orignal expression set 
plateThrCalibration <- as.data.frame(plateThrCalibration)
plateThrCalibration <- round(plateThrCalibration, 1)

# replace column names prior to binding dataframes together to reform original expression set   
colnames(plateOneCalibration) <- colnames(hybMedFile[,25:1329])
colnames(plateTwoCalibration) <- colnames(hybMedFile[,25:1329])
colnames(plateThrCalibration) <- colnames(hybMedFile[,25:1329])

# join dataframes together to reform original expression set 
hybPlateFile <- rbind(plateOneCalibration, plateTwoCalibration, plateThrCalibration)
# add columns to beginning 
hybPlateFile <- cbind(hybMedFile[,1:24], hybPlateFile)
hybPlateFile <- hybPlateFile[order(hybPlateFile$SampleDescription), ]

# hybPlateFile has undergone: HYBRIDISATION, MEDIAN AND CALIBRATION SCALE FACTOR NORMALISATION 
# remove outlier
hybPlateFile <- hybPlateFile[hybPlateFile[[23]] < 2.5 & hybPlateFile[[23]] > 0.4, ]

############################################################################################################################################################
# ------------- COCONUT batch effect removal 
############################################################################################################################################################
# log transform normalised data 
hybPlateFile[,25:1329] <- log2(hybPlateFile[,25:1329])

# split into plates 
mynorm.plate1 <- hybPlateFile[hybPlateFile$PlateId == "Set 001",]
mynorm.plate2 <- hybPlateFile[hybPlateFile$PlateId == "Set 002",]
mynorm.plate3 <- hybPlateFile[hybPlateFile$PlateId == "Set 003",]

# Phenofile is a function that for each plate, creates a dataframe containing the phenotype information for each sample 
phenoFile <- function(file){
  phenodf <- data.frame(file[[14]], file[[18]])
  colnames(phenodf) <- c("SampleDescription", "SampleGroup")
  phenodf$DB <- NA
  phenodf$DV <- NA
  phenodf$HC <- NA
  phenodf$OD <- NA
  phenodf$KD <- NA
  phenodf$TB <- NA
  phenodf$DB[phenodf[[2]] == "DB"] <- 0
  phenodf$DV[phenodf[[2]] == "DV"] <- 0
  phenodf$HC[phenodf[[2]] == "HC"] <- 0
  phenodf$OD[phenodf[[2]] == "OD"] <- 0
  phenodf$KD[phenodf[[2]] == "KD"] <- 0
  phenodf$TB[phenodf[[2]] == "TB"] <- 0
  phenodf[is.na(phenodf)] <- 1
  rownames(phenodf) <- phenodf$SampleDescription
  return(phenodf)
}

# call function for each plate 
plate1_pheno <- phenoFile(mynorm.plate1)
plate2_pheno <- phenoFile(mynorm.plate2)
plate3_pheno <- phenoFile(mynorm.plate3)

# create object for each plate
# this is required as input for the coconut batch effect method 
plate1_obj <- vector()
# must contain phenotype information 
plate1_obj$pheno <- plate1_pheno
# must contain gene information in the transformed format 
plate1_obj$genes <- as.data.frame(t(mynorm.plate1[,25:1329]))
# column names must be the same as the row names for the phenotypes 
colnames(plate1_obj$genes) <- plate1_obj$pheno$SampleDescription

# create object for each plate
# this is required as input for the coconut batch effect method 
plate2_obj <- vector()
plate2_obj$pheno <- plate2_pheno
# must contain protein information in the transformed format 
plate2_obj$genes <- as.data.frame(t(mynorm.plate2[,25:1329]))
# column names must be the same as the row names for the phenotypes 
colnames(plate2_obj$genes) <- plate2_obj$pheno$SampleDescription

# create object for each plate
# this is required as input for the coconut batch effect method 
plate3_obj <- vector()
plate3_obj$pheno <- plate3_pheno
# must contain protein information in the transformed format 
plate3_obj$genes <- as.data.frame(t(mynorm.plate3[,25:1329]))
# column names must be the same as the row names for the phenotypes 
colnames(plate3_obj$genes) <- plate3_obj$pheno$SampleDescription

# concatenate into list and rename objects 
gses <- list(plate1_obj, plate2_obj, plate3_obj)
names(gses) <- c("Plate1", "Plate2", "Plate3")
# call coconut function 
# must state the healthy control column 
cocoRes <- COCONUT(gses, control.0.col = "HC",
                   byPlatform = F)

# combine function simplifies output 
cocoCombine <- combineCOCOoutput(cocoRes)
# pull out normalised protein values 
mynorm_coco <- data.frame(t(cocoCombine$genes))
# replace full sample names
rownames(mynorm_coco) <- cocoCombine$pheno$SampleDescription
mynorm_coco <- mynorm_coco[order(rownames(mynorm_coco)),]

# add on the other columns 
mynorm_coco <- data.frame(hybPlateFile[,1:24], mynorm_coco)
# mynorm_coco has now  undergone: HYBRIDISATION, MEDIAN AND PLATE CALIBRATION NORMALISATION, AND COCONUT BATCH EFFECT REMOVAL


############################################################################################################################################################
# ------------- Post QC process
# ------------- Remove the human-virus proteins and CRP protein 
# ------------- Add ensembl transcript IDs to the dataframe
############################################################################################################################################################

# change the column names from mynorm_coco so that can match up with sequence data file 
colnames(mynorm_coco) <- str_replace_all(colnames(mynorm_coco),"\\.", "-")
colnames(mynorm_coco) <- str_replace_all(colnames(mynorm_coco), "SeqId-", "SeqId.")

# remove the CRP protein and human virus proteins
# pull out the sequence IDs corresponding to the CRP protein and human-virus proteins 
virusID <- as.list(seqData$SeqId[grep("Human-virus", seqData$EntrezGeneSymbol)])
CRP <- seqData$SeqId[grep("P02741", seqData$UniProt)]
proteinsRemove <- list.append(virusID, CRP)

# remove these proteins from the data frames
# qcExp is the expression matrix after QC processes and after CRP/human viruses have been removed 
# qcSeq is the corresponding protein information after those proteins have been removed 
qcExp <- mynorm_coco[, !colnames(mynorm_coco) %in% proteinsRemove]
qcSeq <- seqData[!seqData$SeqId %in% proteinsRemove,]


############# PCA plots of SomaScan dataset
pc_soma_f <- qcExp[qcExp$SampleGroup %in% c("DB", "DV"),]

pc_soma <- prcomp(pc_soma_f[,33:ncol(pc_soma_f)], rank = 6, scale. = T)
s_soma <- summary(pc_soma)

pc_dataset_soma <- data.frame(pc_soma$x,
                              Group = pc_soma_f$SampleGroup,
                              Sex = pc_soma_f$Gender,
                              Age = pc_soma_f$ageMonths)

# disease PCA plots
disease.1.2_soma <- ggplot(pc_dataset_soma, aes(x = PC1, y = PC2, color = Group))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("red", "blue"),
                     labels = c("Bacterial", "Viral"))+
  stat_ellipse()+
  theme(text = element_text(size=12, family="Palatino"), 
        legend.position = 'bottom')+
  guides(color=guide_legend(nrow=2,byrow=TRUE))

disease.3.4_soma <- ggplot(pc_dataset_soma, aes(x = PC3, y = PC4, color = Group))+
  geom_point()+
  theme_bw()+
  scale_color_manual(values = c("red", "blue"),
                     labels = c("Bacterial", "Viral"))+
  stat_ellipse()+
  theme(text = element_text(size=12, family="Palatino"), 
        legend.position = 'bottom')+
  guides(color=guide_legend(nrow=2,byrow=TRUE))

# sex PCA plots
sex.1.2_soma <- ggplot(pc_dataset_soma, aes(x = PC1, y = PC2, color = Sex))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_d(option = "C", end = 0.9)+
  stat_ellipse()+
  theme(text = element_text(size=12, family="Palatino"), 
        legend.position = 'bottom')+
  guides(color=guide_legend(nrow=2,byrow=TRUE))

sex.3.4_soma <- ggplot(pc_dataset_soma, aes(x = PC3, y = PC4, color = Sex))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_d(option = "C", end = 0.9)+
  stat_ellipse()+
  theme(text = element_text(size=12, family="Palatino"), 
        legend.position = 'bottom')+
  guides(color=guide_legend(nrow=2,byrow=TRUE))


# age PCA plots
age.1.2_soma <- ggplot(pc_dataset_soma, aes(x = PC1, y = PC2, color = Age))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_c()+
  # stat_ellipse()+
  theme(text = element_text(size=12, family="Palatino"), 
        legend.position = 'bottom')#+
#guides(color=guide_legend(nrow=2,byrow=TRUE))

age.3.4_soma <- ggplot(pc_dataset_soma, aes(x = PC3, y = PC4, color = Age))+
  geom_point()+
  theme_bw()+
  scale_color_viridis_c()+
  # stat_ellipse()+
  theme(text = element_text(size=12, family="Palatino"), 
        legend.position = 'bottom')#+
# guides(color=guide_legend(nrow=2,byrow=TRUE))