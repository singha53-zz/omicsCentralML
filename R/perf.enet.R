#' cross-validation function for elastic net panel
#'
#' Estimate test error of elastic net panel
#' @param X nxp matrix - training dataset
#' @param Y categorical variables
#' @param alpha = 1 (lasso), alpha = 0 (ridge), 0 < alpha < 1 (elastic net penalty)
#' @param family "binomial" or "multinomial"
#' @param lambda = strength of elastic net penalty
#' @param M = # of folds
#' @param iter = number of times to repeat cross-valiation
#' @param threads - number of cpus (each cross-valiation scheme performed on a separate node)
#' @param progressBar - show progressbar (TRUE/FALE)
#' @export
perf <- function(x, ...){
  UseMethod("perf")
}


perf.enet = function (object, validation = c("Mfold", "loo"), M = 5, iter = 10,
  threads = 4, progressBar = TRUE) {
  library(dplyr)
  library(tidyr)
  X = object$X
  Y = object$Y
  n = nrow(X)
  alpha = object$alpha
  family = object$family
  lambda = object$lambda
  filter = object$filter
  topranked = object$topranked
  keepVar = object$keepVar
  if (validation == "Mfold") {
    folds <- lapply(1:iter, function(i) caret::createFolds(Y, k = M))
    cl <- parallel::makeCluster(mc <- getOption("cl.cores",
      threads))
    parallel::clusterExport(cl, varlist = c("enetCV", "enet",
      "X", "Y", "alpha", "lambda", "M", "folds", "progressBar",
      "family", "filter", "topranked", "keepVar"), envir = environment())
    cv <- parallel::parLapply(cl, folds, function(foldsi,
      X, Y, alpha, lambda, M, progressBar, family, filter,
      topranked, keepVar) {
      omicsCentralML::enetCV(X = X, Y = Y, alpha = alpha, lambda = lambda,
        M = M, folds = foldsi, progressBar = progressBar,
        family = family, filter = filter, topranked = topranked,
        keepVar=keepVar)
    }, X, Y, alpha, lambda, M, progressBar, family, filter,
      topranked, keepVar) %>% amritr::zip_nPure()
    parallel::stopCluster(cl)
    perf <- do.call(rbind, cv$perf) %>% as.data.frame %>%
      gather(ErrName, Err) %>% dplyr::group_by(ErrName) %>%
      dplyr::summarise(Mean = mean(Err), SD = sd(Err))
  }
  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
    cv <- omicsCentralML::enetCV(X, Y, alpha, lambda, M, folds, progressBar,
      family, filter, topranked, keepVar)
    perf <- data.frame(Mean = cv$perf) %>% mutate(ErrName = rownames(.))
    perf$SD <- NA
  }
  result = list()
  result$folds = folds
  result$probs = cv$probs
  result$trueLabels = cv$trueLabels
  result$panels = cv$enet.panel
  result$perf = perf
  return(invisible(result))
}

enetCV = function (X, Y, alpha, lambda, M, folds, progressBar, family,
  filter, topranked, keepVar) {
  probs <- predictResponseList <- enet.panel <- list()
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    X.train = X[-omit, , drop = FALSE]
    Y.train = Y[-omit]
    if (filter == "none") {
      X.train1 <- X.train
    }
    if (filter == "p.value") {
      design <- model.matrix(~Y.train)
      fit <- limma::eBayes(limma::lmFit(t(X.train), design))
      top <- limma::topTable(fit, coef = 2, adjust.method = "BH",
        n = nrow(fit))
      X.train1 <- X.train[, rownames(top)[1:topranked],
        drop = FALSE]
    }
    if(is.null(keepVar)){
      penalty.factor <- rep(1, ncol(X.train1))
      X.train2 <- X.train1
    } else {
      X.train1 <- X.train1[, setdiff(colnames(X.train1), keepVar)]
      X.train2 <- as.matrix(cbind(X.train1, X.train[, keepVar]))
      colnames(X.train2) <- c(colnames(X.train1), keepVar)
      penalty.factor <- c(rep(1, ncol(X.train1)), rep(0, length(keepVar)))
    }


    X.test1 = X[omit, colnames(X.train2), drop = FALSE]
    if (family == "binomial") {
      fit <- glmnet::glmnet(X.train2, Y.train, family = "binomial", alpha = alpha, penalty.factor=penalty.factor)
      cv.fit <- glmnet::cv.glmnet(X.train2, Y.train, family = "binomial")
      if (is.null(lambda)) {
        lambda = cv.fit$lambda.min
      }
      else {
        lambda = lambda
      }
      Coefficients <- coef(fit, s = lambda)
      Active.Index <- which(Coefficients[, 1] != 0)
      Active.Coefficients <- Coefficients[Active.Index,
        ]
      enet.panel[[i]] <- names(Active.Coefficients)[-1]
    }
    if (family == "multinomial") {
      fit <- glmnet::glmnet(X.train2, Y.train, family = "multinomial",
        alpha = alpha, type.multinomial = "grouped", penalty.factor=penalty.factor)
      cv.fit <- glmnet::cv.glmnet(X.train2, Y.train, family = "multinomial")
      if (is.null(lambda)) {
        lambda = cv.fit$lambda.min
      }
      else {
        lambda = lambda
      }
      Coefficients <- coef(fit, s = lambda)
      Active.Index <- which(Coefficients[[1]][, 1] != 0)
      Active.Coefficients <- Coefficients[[1]][Active.Index,
        ]
      enet.panel[[i]] <- names(Active.Coefficients)[-1]
    }
    probs[[i]] <- predict(fit, newx = X.test1, s = lambda,
      type = "response")
    predictResponseList[[i]] <- predict(fit, newx = X.test1,
      s = lambda, type = "class")
  }
  predictResponse <- unlist(predictResponseList)
  if (family == "binomial") {
    probs <- unlist(probs)
    trueLabels = Y[unlist(folds)]
    perf <- omicsCentralML::tperformance(weights = probs, trueLabels = trueLabels)
  }
  else {
    trueLabels = Y[unlist(folds)]
    mat <- table(factor(trueLabels, levels(Y)), factor(predictResponse,
      levels(Y)))
    mat2 <- mat
    diag(mat2) <- 0
    classError <- colSums(mat2)/colSums(mat)
    er <- sum(mat2)/sum(mat)
    ber <- mean(classError)
    perf <- c(classError, er, ber)
    names(perf) <- c(names(classError), "ER", "BER")
  }
  return(list(probs = probs, trueLabels = trueLabels, perf = perf,
    enet.panel = enet.panel, predictResponse = predictResponse))
}
