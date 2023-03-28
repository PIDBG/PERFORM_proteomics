# functions for package - omicsHJ
load.packages <- function(){
  library('ggplot2')
  library(viridis)
  library(ggpubr)
  library(glmnet)
  library(pROC)
  library(limma)
  library(ROCR)
  library(dplyr)
  library(grid)
  library(gridExtra)
  library(stringr)
  library(clipr)
}


# Limma function
limma.fun <- function(data, p.thresh, comparisons, start.col, end.col, n.features, model.mat){
  # construct linear models
  l.m <- lmFit(t(data[,start.col:end.col]), model.mat)
  # for pairwise comparisons between groups, make a contrast matrix
  constrasts.design <- paste(unlist(comparisons), collapse = "-")
  constrasts.mat <- makeContrasts(contrasts = constrasts.design, levels = colnames(model.mat))
  # linear model for age, gender, disease + plate
  cb.fit <- eBayes(contrasts.fit(l.m, constrasts.mat))
  # list of DE features
  top.features <- topTable(cb.fit, coef = 1, adjust = "BH", number = n.features, p.value = p.thresh)
  print(top.features)
  return(top.features)
}


# fspls packages
fspls.packages <- function(){
  library(nnet)
  library(dplyr)
  library(pROC)
  library(ROCR)
  library(glmnet)
}



fspls.binomial <- function(path.fspls, expression, response, p.thresh, beam, split, split.ratio, max){
  source(path.fspls)
  
  fspls.packages()
  ##FOLLOWING JUST SETS DEFAULT VALUES, can change with options if you want to
  # params =readJSONParams(getOption("fspls.default_params"))
  options("fspls.family"= "binomial") ## binomial method 
  options("method" = "fspls")
  options("fspls.lambda" = 0)       ## if 0, no shrinkage. IF specified, then uses that shrinkage, if NULL, then uses CV
  options("fspls.lambda1" = 1)      ## for pvalue adjustment
  options("fspls.debug" = "off")     ## a lot of debugging information printed
  options("fspls.log" = NULL)
  options("fspls.beam" = beam)
  options("fspls.pv_thresh" = p.thresh)
  options("fspls.max" = max)
  elapse_start = proc.time()
  # create input for fspls
  # check if test training split is wanted 
  if(isTRUE(split)){
    print(paste("Splitting by: ", split.ratio, sep = ""))
    ind <- sample.split(response, SplitRatio = split.ratio)
    train.exp <- expression[ind == TRUE,]
    test.exp <- expression[ind == FALSE,]
    train.res <- response[ind == TRUE]
    test.res <- response[ind == FALSE]
    print(train.res)
    print(test.res)
    # create test and training data items 
    train.obj <- list(data = train.exp, y = data.frame(train.res), 
                      weights = data.frame(rep(1, length(train.res))))
    test.obj <- list(data = test.exp, y = data.frame(test.res), 
                     weights = data.frame(rep(1, length(test.res))))
    train.obj <- list("train" = train.obj)
    test.obj <- list("test" = test.obj)
  }
  else if(isFALSE(split)){
    weights = rep(1, length(response))
    obj <- list(data = expression, y = data.frame(response), weights = data.frame(response = weights))
    train.obj <- list("train" = obj)
    test.obj <- list("test" = obj)
  }
  #FOLLOWING RESTRICTS TO VARIABLES WITH NA
  options("fspls.pheno" = dimnames(train.obj[[1]]$y)[[2]][1])
  print("LOADING DATA ELPASED TIME:")
  print(proc.time()-elapse_start)
  ## p-value threshold will determine the number of selected features
  elapse_start = proc.time()
  # run fs-pls
  
  model = trainModel(trainOriginal_l1 = train.obj,testOriginal_l1 = test.obj)#, pv_thresh = p.thresh, max = max)
}

#### fspls with multiple iterations returning top signaure, auc, acc for each iteration 
fspls.iterate <- function(n.iterations, path.fspls, expression, response, max, p.thresh, beam, split, split.ratio, seed){
  accuracy <- list()
  signature <- list()
  aucs <- list()
  aucs_train <- list()
  models <- list()
  df.names <- data.frame()
  for(i in 1:n.iterations){
    seed.i <- seed+i
    set.seed(seed.i)
    print(paste("Iteration",i, sep = " "))
    model <- fspls.binomial(path.fspls = path.fspls,
                            expression = expression, 
                            response = response, p.thresh = p.thresh, 
                            beam = beam, max = max, 
                            split = split, split.ratio = split.ratio)
    # extract highest AUC 
    idx <- length(model)
    ind <- which(model[[idx]][[1]][['testeval']][,'auc'] == max (model[[idx]][[1]][['testeval']][,'auc']))
    #auc.model <- data.frame(test = model[[idx]][[1]][['testeval']][,'auc'], 
    #                         train = model[[idx]][[1]][['traineval']][,'auc'])
    # auc.model$mean <- apply(auc.model, 1, mean)
    #ind <- which(auc.model$mean == max(auc.model$mean))
    print(ind)
    
    idx_genes <- rownames(model[[idx]][[1]][['testeval']])[ind]
    #idx_genes <- str_replace_all(idx_genes, "\\.test", "\\.train")
    print(paste("ID:", idx_genes, sep = " "))
    
    if (length(idx_genes) > 1) {
      idx_sig <- which(model[[idx]][[1]][["testeval"]][idx_genes,'acc']==max(model[[idx]][[1]][["testeval"]][idx_genes,'acc']))[1]
      idx_sig <- idx_genes[idx_sig]
    } else idx_sig <- idx_genes
    print(idx_sig)
    acc <- model[[idx]][[1]]$testeval[rownames(model[[idx]][[1]]$testeval) == idx_sig,"acc"]
    auc <- model[[idx]][[1]]$testeval[rownames(model[[idx]][[1]]$testeval) == idx_sig, 'auc']
    auc_train <- model[[idx]][[1]]$traineval[rownames(model[[idx]][[1]]$testeval) == idx_sig, 'auc']
    names.frame <- model[[idx]][[1]]$variablesl[rownames(model[[idx]][[1]]$testeval) == idx_sig]
    df <- c(table(names.frame), auc = auc)
    pros <- str_split(idx_sig, "\\_\\_")
    pros <- c(unlist(pros), "AUC")
    pros <- str_remove_all(pros, "\\.test")
    df <- t(data.frame(df))
    colnames(df) <- c(pros)
    accuracy[[i]] <- acc
    signature[[i]] <- idx_sig    
    aucs[[i]] <- auc
    models[[i]] <- model
    aucs_train[[i]] <- auc_train
    #print(df.names)
    df <- data.frame(df)
    df.names <- bind_rows(df.names, df)
  }
  return(list(acc = accuracy, sig = signature, auc = aucs, model = models, table = df.names, training_aucs = aucs_train))
}


