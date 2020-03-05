#' Elastic net classification panel
#'
#' build an elastic net classification panel
#' @param X nxp matrix - training dataset
#' @param Y categorical variables
#' @param alpha = 1 (lasso), alpha = 0 (ridge), 0 < alpha < 1 (elastic net penalty)
#' @param lambda = strength of elastic net penalty
#' @param family "binomial" or "multinomial"
#' @param X.test nxp matrx  - test dataset
#' @param filter = "none" or "p.value"
#' @param topranked = 50 (top number of features to select and build a classifier)
#' @export
enet = function (X, Y, alpha, lambda = NULL, family, X.test = NULL,
  Y.test = NULL, filter = "p.value", topranked = 50, keepVar = NULL){
  if (filter == "none") {
    X1 <- X
  }
  if (filter == "p.value") {
    design <- model.matrix(~Y)
    fit <- limma::eBayes(limma::lmFit(t(X), design))
    top <- limma::topTable(fit, coef = 2, adjust.method = "BH",
      n = nrow(fit))
    X1 <- X[, rownames(top)[1:topranked]]
  }
  if(is.null(keepVar)){
    penalty.factor <- rep(1, ncol(X1))
    X2 <- X1
  } else {
    X1 <- X1[, setdiff(colnames(X1), keepVar)]
    X2 <- as.matrix(cbind(X1, X[, keepVar]))
    colnames(X2) <- c(colnames(X1), keepVar)
    penalty.factor <- c(rep(1, ncol(X1)), rep(0, length(keepVar)))
  }

  if (family == "binomial") {
    fit <- glmnet::glmnet(X2, Y, family = "binomial", alpha = alpha, penalty.factor=penalty.factor)
    cv.fit <- glmnet::cv.glmnet(X2, Y, family = "binomial")
    if (is.null(lambda)) {
      lambda = cv.fit$lambda.min
    } else {
      lambda = lambda
    }
    Coefficients <- coef(fit, s = lambda)
    Active.Index <- which(Coefficients[, 1] != 0)
    Active.Coefficients <- Coefficients[Active.Index, ]
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
  }
  if (family == "multinomial") {
    fit <- glmnet::glmnet(X2, Y, family = "multinomial", alpha = alpha,
      type.multinomial = "grouped", penalty.factor=penalty.factor)
    cv.fit <- glmnet::cv.glmnet(X2, Y, family = "multinomial")
    if (is.null(lambda)) {
      lambda = cv.fit$lambda.min
    } else {
      lambda = lambda
    }
    Coefficients <- coef(fit, s = lambda)
    Active.Index <- which(Coefficients[[1]][, 1] != 0)
    Active.Coefficients <- Coefficients[[1]][Active.Index,
      ]
    enet.panel <- names(Active.Coefficients)[-1]
    enet.panel.length <- length(enet.panel)
  }
  if (!is.null(X.test)) {
    probs <- predict(fit, newx = X.test[, colnames(X2)], s = lambda, type = "response")
    predictResponse <- unlist(predict(fit, newx = X.test,
      s = lambda, type = "class"))
    if (family == "binomial") {
      perfTest <- tperformance(weights = as.numeric(as.matrix(probs)),
        trueLabels = Y.test)
    } else {
      mat <- table(factor(as.character(predictResponse),
        levels = levels(Y.test)), Y.test)
      mat2 <- mat
      diag(mat2) <- 0
      classError <- colSums(mat2)/colSums(mat)
      er <- sum(mat2)/sum(mat)
      ber <- mean(classError)
      perfTest <- c(classError, er, ber)
      names(perfTest) <- c(names(classError), "ER", "BER")
    }
  } else {
    perfTest <- predictResponse <- probs <- NA
  }
  result <- list(X = X, Y = Y, fit = fit, enet.panel = enet.panel,
    lambda = lambda, alpha = alpha, family = family, probs = probs,
    Active.Coefficients = Active.Coefficients, perfTest = perfTest,
    predictResponse = predictResponse, filter = filter, topranked = topranked, keepVar=keepVar)
  class(result) <- "enet"
  return(result)
}
