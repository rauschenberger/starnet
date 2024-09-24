
#--- Main function -------------------------------------------------------------

#' @export
#' @title
#' Stacked Elastic Net Regression
#' 
#' @description
#' Implements stacked elastic net regression.
#'  
#' @param y
#' response\strong{:}
#' numeric vector of length \eqn{n}
#' 
#' @param X
#' covariates\strong{:}
#' numeric matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#'
#' @param family
#' character "gaussian", "binomial" or "poisson"
#'
#' @param nalpha
#' number of \code{alpha} values
#' 
#' @param alpha
#' elastic net mixing parameters\strong{:}
#' vector of length \code{nalpha} with entries
#' between \eqn{0} (ridge) and \eqn{1} (lasso);
#' or \code{NULL} (equidistance)
#'
#' @param nfolds
#' number of folds
#'
#' @param foldid
#' fold identifiers\strong{:}
#' vector of length \eqn{n} with entries between \eqn{1} and \code{nfolds};
#' or \code{NULL} (balance)
#' 
#' @param type.measure
#' loss function\strong{:}
#' character "deviance", "class", "mse" or "mae"
#' (see \code{\link[glmnet]{cv.glmnet}})
#'
#' @param alpha.meta
#' meta-learner\strong{:}
#' value between \eqn{0} (ridge) and \eqn{1} (lasso)
#' for elastic net regularisation; 
#' \code{NA} for convex combination
#' 
#' @param intercept,upper.limit,unit.sum
#' settings for meta-learner\strong{:} logical,
#' or \code{NULL}
#' (\code{intercept=!is.na(alpha.meta)},
#' \code{upper.limit=TRUE},
#' \code{unit.sum=is.na(alpha.meta)})
#' 
#' @param penalty.factor
#' differential shrinkage\strong{:}
#' vector of length \eqn{n} with entries
#' between \eqn{0} (include) and \eqn{Inf} (exclude), 
#' or \code{NULL} (all \eqn{1})
#' 
#' @param ...
#' further arguments passed to \code{\link[glmnet]{glmnet}}
#' 
#' @inherit starnet-package references
#' 
#' @details
#' Post hoc feature selection\strong{:} consider
#' argument \code{nzero} in functions
#' \code{\link{coef}} and \code{\link{predict}}.
#' 
#' @return
#' Object of class \code{\link[starnet]{starnet}}.
#' The slots \code{base} and \code{meta}
#' contain \code{\link[glmnet]{cv.glmnet}}-like objects,
#' for the base and meta learners, respectively.
#' 
#' @examples
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' set.seed(1)
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X,family="gaussian")}}
#' \donttest{
#' set.seed(1)
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X,family="gaussian")}
#' 
starnet <- function(y,X,family="gaussian",nalpha=21,alpha=NULL,nfolds=10,foldid=NULL,type.measure="deviance",alpha.meta=1,penalty.factor=NULL,intercept=NULL,upper.limit=NULL,unit.sum=NULL,...){
  
  if(is.na(alpha.meta) && (!"CVXR" %in% .packages(all.available=TRUE))){
    stop("Install CVXR from CRAN for alpha.meta=NA.",call.=FALSE)
  }
  
  #--- temporary ---
  # family <- "gaussian"; nalpha <- 21; alpha <- NULL; nfolds <- 10; foldid <- NULL; type.measure <- "deviance"; alpha.meta <- 0; penalty.factor <- NULL; intercept <- TRUE; upper.limit=TRUE; unit.sum=FALSE
  
  #--- default ---
  #if(all(is.null(intercept),is.null(upper.limit),is.null(unit.sum))){
  #  if(is.na(alpha.meta)){
  #    intercept <- FALSE
  #    upper.limit <- unit.sum <- TRUE
  #  } else if(alpha.meta==1){
  #    intercept <- TRUE
  #    upper.limit <- unit.sum <- FALSE
  #  } else if(alpha.meta==0){
  #    intercept <- upper.limit <- TRUE
  #    unit.sum <- FALSE
  #  }
  #}
  
  if(is.null(intercept)){intercept <- !is.na(alpha.meta)}
  if(is.null(upper.limit)){upper.limit <- TRUE}  
  if(is.null(unit.sum)){unit.sum <- is.na(alpha.meta)}
    
  #--- checks ---
  cornet:::.check(x=y,type="vector")
  cornet:::.check(x=X,type="matrix")
  cornet:::.check(x=family,type="string",values=c("gaussian","binomial","poisson"))
  cornet:::.check(x=nalpha,type="scalar",min=0)
  cornet:::.check(x=alpha,type="vector",min=0,max=1,null=TRUE)
  if(length(y)!=nrow(X)){stop("Contradictory sample size.",call.=FALSE)}
  cornet:::.check(x=nfolds,type="scalar",min=3)
  cornet:::.check(x=foldid,type="vector",values=seq_len(nfolds),null=TRUE)
  cornet:::.check(x=type.measure,type="string",values=c("deviance","class","mse","mae")) # not auc (min/max confusion)
  if(any(!is.na(alpha.meta))){cornet:::.check(x=alpha.meta[!is.na(alpha.meta)],type="vector",min=0,max=1)}
  if(!is.null(c(list(...)$lower.limits,list(...)$upper.limits))){
    stop("Reserved arguments \"lower.limits\" and \"upper.limits\".",call.=FALSE)
  }
  if(!is.null(list(...)$intercept)){
    stop("Reserved argument \"intercept\".",call.=FALSE)
  }
  if(is.null(penalty.factor)){
    penalty.factor <- rep(x=1,times=ncol(X))
    # never pass NULL to glmnet/cv.glmnet!
  }
  
  #--- fold identifiers ---
  if(is.null(foldid)){
    foldid <- .folds(y=y,nfolds=nfolds)
  }
  nfolds <- max(foldid)
    
  #--- alpha sequence ---
  if(is.null(alpha)){
    alpha <- seq(from=0,to=1,length.out=nalpha)
  } else {
    alpha <- alpha
    nalpha <- length(alpha)
  }
  
  #--- full fit ---
  base <- lapply(alpha,function(x) list())
  names(base) <- paste0("alpha",alpha)
  nlambda <- numeric()
  for(i in seq_len(nalpha)){
    base[[i]]$glmnet.fit <- glmnet::glmnet(y=y,x=X,family=family,alpha=alpha[i],penalty.factor=penalty.factor,...)
    base[[i]]$lambda <- base[[i]]$glmnet.fit$lambda
    nlambda[i] <- length(base[[i]]$glmnet.fit$lambda)
  }
  
  #--- predictions ---
  n <- length(y)
  link <- list()
  for(i in seq_len(nalpha)){
    link[[i]] <- matrix(data=NA,nrow=n,ncol=nlambda[i])
  }
  
  #--- base cross-validation ---
  for(k in seq_len(nfolds)){
    y0 <- y[foldid!=k]
    y1 <- y[foldid==k]
    X0 <- X[foldid!=k,,drop=FALSE]
    X1 <- X[foldid==k,,drop=FALSE]
    for(i in seq_len(nalpha)){
      object <- glmnet::glmnet(y=y0,x=X0,family=family,alpha=alpha[i],penalty.factor=penalty.factor,...)
      temp <- stats::predict(object=object,newx=X1,type="link",
                          s=base[[i]]$glmnet.fit$lambda)
      link[[i]][foldid==k,seq_len(ncol(temp))] <- temp
    }
  } 
  
  #--- tune base lambdas ---
  for(i in seq_len(nalpha)){
    fit <- .mean.function(link[[i]],family=family)
    base[[i]]$cvm <- apply(fit,2,function(x) .loss(y=y,x=x,family=family,type.measure=type.measure,foldid=foldid))
    base[[i]]$id.min <- which.min(base[[i]]$cvm)
    base[[i]]$lambda.min <- base[[i]]$lambda[base[[i]]$id.min]
  }
  
  #--- predictions ---
  hat <- matrix(NA,nrow=nrow(X),ncol=nalpha)
  for(i in seq_len(nalpha)){
    hat[,i] <- link[[i]][,base[[i]]$id.min]
  }
  
  meta <- list()
  
  #--- meta cross-validation ---
  for(i in seq_along(alpha.meta)){
    if(is.na(alpha.meta[i])){next}
    meta[[i]] <- glmnet::cv.glmnet(y=y,x=hat,
                                     lower.limits=0,
                                     upper.limits=ifelse(upper.limit,1,Inf), # was 1
                                     foldid=foldid,
                                     family=family,
                                     type.measure=type.measure,
                                     intercept=intercept,
                                     alpha=alpha.meta[i],...)
    if(unit.sum){
      warning("inequality constraint",call.=FALSE)
      cond <- Matrix::colSums(meta[[i]]$glmnet.fit$beta)>1
      meta[[i]]$cvm[cond] <- Inf
      meta[[i]]$lambda.min <- meta[[i]]$lambda[which.min(meta[[i]]$cvm)]
    }
  }
  
  #--- convex combination ---  
  for(i in seq_along(alpha.meta)){
    if(!is.na(alpha.meta[i])){next}
    glm <- .glm(y=y,X=hat,family=family,intercept=intercept,lower.limit=TRUE,upper.limit=upper.limit,unit.sum=unit.sum)
    meta[[i]] <- list()
    meta[[i]]$glmnet.fit <- list()
    meta[[i]]$glmnet.fit$a0 <- rep(x=glm$alpha,times=2)
    meta[[i]]$glmnet.fit$beta <- cbind(0,glm$beta) # matrix(data=rep(glm$beta,times=2),ncol=2)
    meta[[i]]$glmnet.fit$lambda <- meta[[i]]$lambda <- c(99e99,0) # was c(1,0)
    meta[[i]]$glmnet.fit$offset <- FALSE
    if(length(alpha.meta)>1){
      link <- rep(NA,times=length(y))
      for(k in seq_len(nfolds)){
        glm <- .glm(y=y[foldid!=k],X=hat[foldid!=k,],family=family,intercept=intercept,lower.limit=TRUE,upper.limit=upper.limit,unit.sum=unit.sum)
        link[foldid==k] <- glm$alpha + hat[foldid==k,] %*% glm$beta
      }
      y_hat <- .mean.function(x=link,family=family)
      cvm <- .loss(y=y,x=y_hat,family=family,type.measure=type.measure,foldid=foldid)
    } else {
      cvm <- 0
    }
    meta[[i]]$cvm <- c(Inf,cvm)
    meta[[i]]$lambda.min <- 0
    class(meta[[i]]) <- "cv.glmnet"
    class(meta[[i]]$glmnet.fit) <- "glmnet"
  }
  
  names(meta) <- paste0("alpha",alpha.meta)
  
  #--- tune meta alpha ---
  cvm_meta <- sapply(meta,function(x) min(x$cvm))
  #message(paste0(names(cvm_meta)," ",round(cvm_meta,digits=2)," "))
  id_meta <- which.min(cvm_meta)
  
  #--- message ---
  meta <- meta[[names(id_meta)]]
  #message(paste(paste(round(stats::coef(meta,s="lambda.min")[1],digits=3)),"_",
  #              paste(round(stats::coef(meta,s="lambda.min")[-1],digits=3),collapse=" ")))
  
  if(FALSE){
    ## debugging: stacking returns same output as tuning
    # cvm <- sapply(base,function(x) x$cvm[x$id.min])
    # id <- which.min(cvm)
    # meta$glmnet.fit$a0[] <- 0
    # meta$glmnet.fit$beta[,] <- 0
    # meta$glmnet.fit$beta[id,] <- 1
  }
  
  #--- return ---
  info <- list(id_meta=id_meta,
               type.measure=type.measure,
               family=family,
               mean=mean(y),
               foldid=foldid,X=X)
  
  list <- list(base=base,meta=meta,info=info)
  class(list) <- "starnet"
  return(list)
}

#--- Methods for class "starnet" ------------------------------------------------

#' @export
#' @title
#' Makes Predictions
#'
#' @description
#' Predicts outcome from features with stacked model.
#' 
#' @param object
#' \link[starnet]{starnet} object
#' 
#' @param newx
#' covariates\strong{:}
#' numeric matrix with \eqn{n} rows (samples)
#' and \eqn{p} columns (variables)
#' 
#' @param type
#' character "link" or "response"
#' 
#' @param nzero
#' maximum number of non-zero coefficients\strong{:}
#' positive integer, or \code{NULL}
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @return
#' Matrix of predicted values, with samples in the rows,
#' and models in the columns. Included models are
#' \code{alpha} (fixed elastic net),
#' \code{ridge} (i.e. \code{alpha0}),
#' \code{lasso} (i.e. \code{alpha1}),
#' \code{tune} (tuned elastic net),
#' \code{stack} (stacked elastic net),
#' and \code{none} (intercept-only model).
#' 
#' @examples
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' set.seed(1)
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' y_hat <- predict(object,newx=X[c(1),,drop=FALSE])}}
#' \donttest{
#' set.seed(1)
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' y_hat <- predict(object,newx=X[c(1),,drop=FALSE])}
#' 
predict.starnet <- function(object,newx,type="response",nzero=NULL,...){
  if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
  
  x <- object
  
  cornet:::.check(x=newx,type="matrix")
  
  # all alphas
  nalpha <- length(x$base)
  link <- matrix(NA,nrow=nrow(newx),ncol=nalpha)
  colnames(link) <- names(x$base)
  for(i in seq_len(nalpha)){
    link[,i] <- as.numeric(stats::predict(object=x$base[[i]]$glmnet.fit,
                            newx=newx,s=x$base[[i]]$lambda.min,type="link"))
  }
  
  # elastic net
  cvm <- sapply(x$base[seq_len(nalpha)],function(x) x$cvm[x$id.min])
  id <- which.min(cvm)
  
  frame <- as.data.frame(link)
  frame$ridge <- frame$alpha0
  frame$lasso <- frame$alpha1
  frame$tune <- frame[[names(id)]]
  if(any(frame$tune!=frame[[id]])){stop("Invalid.")}
  if(is.null(nzero)){
    frame$stack <- as.numeric(stats::predict(object=x$meta,newx=link,type="link",s="lambda.min"))
  } else {
    coef <- coef.starnet(object=x,nzero=nzero)
    frame$stack <- as.numeric(coef$alpha + newx %*% coef$beta)
  }
  frame$none <- as.numeric(stats::predict(object=x$base$alpha1$glmnet.fit,
                  newx=newx,s=Inf,type="link"))
  
  if(type=="link"){
    return(frame)
  } else  if(type=="response"){
    frame <- lapply(frame,function(y) .mean.function(y,family=x$info$family))
    return(as.data.frame(frame))
  } else {
    stop("Invalid type.",call.=FALSE)
  }
  
}

#' @export
#' @title
#' Extract Coefficients
#'
#' @description
#' Extracts pooled coefficients.
#' (The meta learners weights the coefficients from the base learners.)
#' 
#' @inheritParams predict.starnet
#' 
#' @return
#' List of scalar \code{alpha} and vector \code{beta},
#' containing the pooled intercept and the pooled slopes,
#' respectively.
#' 
#' @examples
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' set.seed(1)
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' coef <- coef(object)}}
#' \donttest{
#' set.seed(1)
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' coef <- coef(object)}
#' 
coef.starnet <- function(object,nzero=NULL,...){
  if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
  
  # base coefficients
  base <- list()
  coef <- sapply(object$base,function(x) stats::coef(object=x$glmnet.fit,s=x$lambda.min))
  base$alpha <- sapply(coef,function(x) x[1,])
  base$beta <- sapply(coef,function(x) x[-1,])
  
  # meta coefficients
  meta <- list()
  weights <- weights.starnet(object)
  meta$alpha <- weights[1]
  meta$beta <- weights[-1]
  
  # pooled coefficients
  pool <- list()
  pool$alpha <- meta$alpha + sum(meta$beta * base$alpha)
  pool$beta <- base$beta %*% meta$beta
  
  # post hoc selection
  if(!is.null(nzero)){
     eta <- as.numeric(pool$alpha + object$info$X %*% pool$beta)
     if(stats::sd(eta)==0){return(list(alpha=pool$alpha,beta=0*pool$beta))}
     model <- .cv.glmnet(y=eta,x=object$info$X,family="gaussian",
                 alpha=1,penalty.factor=1/abs(pool$beta),
                 foldid=object$info$foldid,nzero=nzero)
     coef <- stats::coef(model,s="lambda.min")
     return(list(alpha=coef[1],beta=coef[-1]))
     # alternatives: lasso-like elastic net (alpha=0.95);
     # logistic regression of predicted probabilities on X
  }
  
  return(pool)
}

#' @export
#' @importFrom stats weights
#' @title
#' Extract Weights
#'
#' @description
#' Extracts coefficients from the meta learner,
#' i.e. the weights for the base learners.
#' 
#' @param object
#' \link[starnet]{starnet} object
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @return
#' Vector containing intercept and slopes from the meta learner.
#' 
#' @examples
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' set.seed(1)
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' weights(object)}}
#' \donttest{
#' set.seed(1)
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' weights(object)}
#' 
weights.starnet <- function(object,...){
  if(length(list(...))!=0){warning("Ignoring argument.",call.=FALSE)}
  x <- object$meta
  coef <- stats::coef(object=x,s=x$lambda.min)
  names <- rownames(coef)
  coef <- as.numeric(coef)
  names(coef) <- names
  return(coef)
}

#' @export
#' @title
#' Print Values
#'
#' @description
#' Prints object of class \link[starnet]{starnet}.
#' 
#' @param x
#' \link[starnet]{starnet} object
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @return
#' Prints "stacked gaussian/binomial/poisson elastic net".
#' 
#' @examples
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' set.seed(1)
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' print(object)}}
#' \donttest{
#' set.seed(1)
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' print(object)}
#' 
print.starnet <- function(x,...){
  cat(paste0("stacked \"",x$info$family,"\" elastic net"),"\n")
}

#--- Manuscript functions ------------------------------------------------------

#' @export
#' @title
#' Model comparison
#'
#' @description
#' Compares stacked elastic net, tuned elastic net, ridge and lasso.
#' 
#' @inheritParams starnet
#' 
#' @param nzero
#' number of non-zero coefficients\strong{:}
#' scalar/vector including positive integer(s) or \code{NA};
#' or \code{NULL} (no post hoc feature selection)
#' 
#' @param nfolds.ext,nfolds.int,foldid.ext,foldid.int
#' number of folds (\code{nfolds})\strong{:}
#' positive integer;
#' fold identifiers (\code{foldid})\strong{:} 
#' vector of length \eqn{n} with entries between \eqn{1} and \code{nfolds},
#' or \code{NULL},
#' for hold-out (single split) instead of cross-validation (multiple splits)\strong{:}
#' set \code{foldid.ext} to \eqn{0} for training and to \eqn{1} for testing samples
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @return
#' List containing the cross-validated loss
#' (or out-of sample loss if \code{nfolds.ext} equals two,
#' and \code{foldid.ext} contains zeros and ones).
#' The slot \code{meta} contains the loss from the stacked elastic net
#' (\code{stack}), the tuned elastic net (\code{tune}), ridge, lasso,
#' and the intercept-only model (\code{none}).
#' The slot \code{base} contains the loss from the base learners.
#' And the slot \code{extra} contains the loss from the restricted
#' stacked elastic net (\code{stack}), lasso, and lasso-like elastic net
#' (\code{enet}),
#' with the maximum number of non-zero coefficients shown in the column name.
#' 
#' @examples
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' set.seed(1)
#' n <- 50; p <- 20
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' y <- rnorm(n=n,mean=rowSums(X[,1:20]))
#' loss <- cv.starnet(y=y,X=X,nfolds.ext=2,nfolds.int=3)}}
#' \donttest{
#' loss <- cv.starnet(y=y,X=X)}
#' 
cv.starnet <- function(y,X,family="gaussian",nalpha=21,alpha=NULL,nfolds.ext=10,nfolds.int=10,foldid.ext=NULL,foldid.int=NULL,type.measure="deviance",alpha.meta=1,nzero=NULL,intercept=NULL,upper.limit=NULL,unit.sum=NULL,...){
  
  # family <- "gaussian"; nfolds.ext <- nfolds.int <- 10; foldid.ext <- foldid.int <- NULL; type.measure <- "deviance"; alpha.meta <- 0; alpha = NULL; nalpha <- 21; nzero <- NULL; intercept <- upper.limit <- unit.sum <- NULL
  
  #--- fold identifiers ---
  if(is.null(foldid.ext)){
      fold <- .folds(y=y,nfolds=nfolds.ext)
  } else {
    fold <- foldid.ext
  }
  nfolds.ext <- max(fold)
  # start trial
  if(is.null(foldid.int)){
    foldid.int <- .folds(y=y,nfolds=nfolds.int)
  }
  nfolds.int <- max(foldid.int)
  # end trial

  #--- alpha sequence ---
  if(is.null(alpha)){
    alpha <- seq(from=0,to=1,length.out=nalpha)
  } else {
    nalpha <- length(alpha)
  }
  
  #--- cross-validated loss ---
  
    meta <- c("stack","tune","ridge","lasso","none")
    base <- paste0("alpha",alpha)
    cols <- c(meta,base)
    if(!is.null(nzero)){cols <- c(cols,paste0(rep(c("stack","lasso","enet"),times=length(nzero)),rep(nzero,each=3)))}
    pred <- matrix(data=NA,nrow=length(y),ncol=length(cols),
                   dimnames=list(NULL,cols))
    for(i in seq_len(nfolds.ext)){
      # dense models
      fit <- starnet(y=y[fold!=i],X=X[fold!=i,],alpha=alpha,nalpha=nalpha,nfolds=nfolds.int,foldid=foldid.int[fold!=i],family=family,type.measure=type.measure,alpha.meta=alpha.meta,intercept=intercept,upper.limit=upper.limit,unit.sum=unit.sum,...) # INSERT ,...
      temp <- predict.starnet(fit,newx=X[fold==i,,drop=FALSE])
      for(j in c(meta,base)){
        pred[fold==i,j] <- unlist(temp[[j]])
      }
      # sparse models
      if(!is.null(nzero)){
        for(j in seq_along(nzero)){
          pred[fold==i,paste0("stack",nzero[j])] <- predict.starnet(fit,newx=X[fold==i,,drop=FALSE],nzero=nzero[j])[,"stack"] # added drop=FALSE 2020-04-02
          temp <- .cv.glmnet(y=y[fold!=i],x=X[fold!=i,],alpha=1,family=family,type.measure=type.measure,foldid=fit$info$foldid,nzero=nzero[j],...)
          pred[fold==i,paste0("lasso",nzero[j])] <- stats::predict(temp,
                  newx=X[fold==i,,drop=FALSE],type="response",s="lambda.min") # added drop=FALSE 2020-04-02
          temp <- .cv.glmnet(y=y[fold!=i],x=X[fold!=i,],alpha=0.95,family=family,type.measure=type.measure,foldid=fit$info$foldid,nzero=nzero[j],...)
          pred[fold==i,paste0("enet",nzero[j])] <- stats::predict(temp,
                  newx=X[fold==i,,drop=FALSE],type="response",s="lambda.min") # added drop=FALSE 2020-04-02
        }
      }
    }
  
  if(length(type.measure)!=1){stop("Implement multiple type measures!")}
  loss <- apply(pred[fold!=0,],2,function(x) .loss(y=y[fold!=0],x=x,family=family,type.measure=type.measure,foldid=fold))
  
  ylim <- range(c(loss[base],loss[c("tune","stack")]))
  tryCatch(graphics::plot(y=loss[base],x=as.numeric(substring(text=base,first=6)),xlab="alpha",ylab="loss",ylim=ylim),error=function(x) NULL)
  tryCatch(graphics::abline(h=loss["tune"],lty=2,col="grey"),error=function(x) NULL)
  tryCatch(graphics::abline(h=loss["stack"],lty=2,col="black"),error=function(x) NULL)
  
  meta <- loss[names(loss) %in% meta]
  base <- loss[names(loss) %in% base]
  
  type <- c("stack","lasso","enet")
  extra <- matrix(NA,nrow=3,ncol=length(nzero),dimnames=list(type,nzero))
  for(i in seq_len(3)){
    for(j in seq_along(nzero)){
      extra[i,j] <- loss[paste0(type[i],nzero[j])]
    }
  }
  
  return(list(meta=meta,base=base,extra=extra))
}

#' @name .simulate
#' 
#' @title
#' Simulation
#'
#' @description
#' Functions for simulating data
#' 
#' @param n
#' sample size\strong{:}
#' positive integer
#' 
#' @param p
#' dimensionality\strong{:}
#' positive integer
#' 
#' @param mode
#' character \code{"sparse"}, \code{"dense"} or \code{"mixed"}
#' 
#' @param family
#' character \code{"gaussian"}, \code{"binomial"} or \code{"poisson"}
#' 
#' @return
#' List of vector \code{y} and matrix \code{X}.
#' 
#' @examples
#' NA
#' 
.simulate.block <- function(n,p,mode,family="gaussian"){
  Z <- matrix(data=stats::rnorm(n*3),nrow=n,ncol=3)
  y <- rowSums(Z)
  #if(family=="binomial"){y <- round(1/(1+exp(-y)))}
  #if(family=="poisson"){y <- round(exp(y))}
  y <- switch(family,gaussian=y,binomial=round(1/(1+exp(-y))),
              poisson=round(exp(y)),stop("Invalid."))
  X <- matrix(data=stats::rnorm(n*p),nrow=n,ncol=p)
  if(mode=="sparse"){
    X[,1] <- sqrt(0.1)*X[,1] + sqrt(0.9)*Z[,1]
    X[,p] <- sqrt(0.1)*X[,p] + sqrt(0.9)*Z[,2]
  } else if(mode=="dense"){
    X[,1:250] <- sqrt(0.9)*X[,1:250] + sqrt(0.1)*Z[,1]
    X[,(p-250+1):p] <- sqrt(0.9)*X[,(p-250+1):p] + sqrt(0.1)*Z[,2]
  } else if(mode=="mixed"){
    X[,1:25] <- sqrt(0.5)*X[,1:25] + sqrt(0.5)*Z[,1]
    X[,(p-25+1):p] <- sqrt(0.5)*X[,(p-25+1):p] + sqrt(0.5)*Z[,2]
  } else {
    stop("Invalid.",.call=FALSE)
  }
  return(list(y=y,X=X))
}
#' @rdname .simulate
#' @param rho
#' correlation\strong{:}
#' numeric between \eqn{0} and \eqn{1}
#' @param pi
#' effects\strong{:}
#' numeric between \eqn{0} (sparse) and \eqn{1} (dense)
#' 
.simulate.grid <- function(n,p,rho,pi,family="gaussian"){
  if(!"mvtnorm" %in% .packages(all.available=TRUE)){
    stop("Install mvtnorm from CRAN for simulation.",call.=FALSE)
  }
  mean <- rep(x=0,times=p)
  sigma <- matrix(data=NA,nrow=p,ncol=p)
  sigma <- rho^abs(row(sigma)-col(sigma))
  diag(sigma) <- 1
  X <- mvtnorm::rmvnorm(n=n,mean=mean,sigma=sigma)
  nzero <- round(pi*p)
  beta <- abs(stats::rnorm(p))*sample(x=rep(x=c(0,1),times=c(p-nzero,nzero)))
  eta <- as.numeric(X %*% beta)
  y <- eta + stats::rnorm(n=n,sd=0.5*stats::sd(eta))
  y <- switch(family,gaussian=y,binomial=round(1/(1+exp(-y))),stop("Invalid."))
  return(list(y=y,X=X,beta=beta))
}
#' @rdname .simulate
#' 
.simulate.mode <- function(n,p,mode,family="gaussian"){
  if(!"mvtnorm" %in% .packages(all.available=TRUE)){
    stop("Install mvtnorm from CRAN for simulation.",call.=FALSE)
  }
  mean <- rep(x=0,times=p)
  sigma <- matrix(data=0.1,nrow=p,ncol=p)
  diag(sigma) <- 1
  X <- mvtnorm::rmvnorm(n=n,mean=mean,sigma=sigma)
  q <- switch(mode,sparse=5,dense=50,mixed=20,stop("Invalid.",.call=FALSE))
  # mixed: 15 is too lasso-friendly
  beta <- sample(rep(c(0,1),times=c(p-q,q)))
  eta <- as.numeric(X %*% beta)
  y <- eta + stats::rnorm(n=n,sd=0.5*stats::sd(eta))
  y <- switch(family,gaussian=y,binomial=round(1/(1+exp(-y))),stop("Invalid."))
  return(list(y=y,X=X,beta=beta))
}

#' @title 
#' Loss
#' 
#' @description
#' Calculate loss from predicted and observed values
#' 
#' @param y
#' observed values\strong{:}
#' numeric vector of length \eqn{n}
#' 
#' @param x
#' predicted values\strong{:}
#' numeric vector of length \eqn{n}
#' 
#' @param family
#' character \code{"gaussian"}, \code{"binomial"}, \code{"poisson"},
#' \code{"mgaussian"}, or \code{"multinomial"} (to implement: \code{"cox"})
#' 
#' @param type.measure
#' character \code{"deviance"}, \code{"mse"}, \code{"mae"}, \code{"class"},
#' or \code{"auc"}
#' 
#' @param foldid
#' fold identifiers\strong{:}
#' integer vector of length \eqn{n},
#' or \code{NULL}
#' 
#' @param grouped
#' logical (for \code{"cox"} only)
#'
#' @examples
#' NA
#' 
.loss <- function(y,x,family,type.measure,foldid=NULL,grouped=TRUE){
  
  if(family=="cox" & grouped){stop("Implement \"grouped Cox\"! See unit tests.",call.=FALSE)}
  
  if(is.null(foldid)&(family=="cox"|type.measure=="auc")){
    stop("Missing foldid.",call.=FALSE)
  }
  
  if(family=="multinomial"){
    Y <- matrix(0,nrow=length(y),ncol=length(unique(y)))
    Y[cbind(seq_along(y),as.numeric(as.factor(y)))] <- 1
  }
  
  if(family=="mgaussian"){
    Y <- y
  }
  
  if(type.measure=="deviance"){
    
    if(family=="gaussian"){
      
      type.measure <- "mse"
      
    } else if(family=="binomial"){
      
      limit <- 1e-05
      x[x<limit] <- limit
      x[x>1-limit] <- 1-limit
      return(mean(-2*(y*log(x)+(1-y)*log(1-x))))
      
    } else if(family=="poisson"){
      
      return(mean(2*(ifelse(y==0,0,y*log(y/x))-y+x)))
             
    } else if(family=="cox"){
      
      cvraw <- numeric()
      for(i in seq_along(unique(foldid))){
        if(grouped){
          stop("Invalid \"grouped Cox\"! See unit tests!",call.=FALSE)
          full <- glmnet::coxnet.deviance(pred=x,y=y)
          mink <- glmnet::coxnet.deviance(pred=x[foldid!=i],y=y[foldid!=i])
          cvraw[i] <- full-mink
        } else {
          cvraw[i] <- glmnet::coxnet.deviance(pred=x[foldid==i],y=y[foldid==i])
        }
      }
      weights <- tapply(X = y[, "status"], INDEX = foldid, FUN = sum)
      return(stats::weighted.mean(x = cvraw/weights, w = weights, na.rm = TRUE))
      
    } else if(family=="multinomial"){
      
      limit <- 1e-05
      x[x<limit] <- limit
      x[x>1-limit] <- limit
      lp <- Y * log(x)
      ly <- Y * log(Y)
      ly[Y==0] <- 0
      return(mean(rowSums(2*(ly-lp))))
      
    }
  }
  
  if(type.measure=="mse" & family %in% c("gaussian","binomial","poisson")){
    return((1 + (family=="binomial"))*mean((x-y)^2))
  }
  
  if(type.measure %in% c("deviance","mse") & family %in% c("mgaussian","multinomial")){
    return(mean(rowSums((Y-x)^2)))
  } 
  
  if(type.measure=="mae" & family %in% c("mgaussian","multinomial")){
    return(mean(rowSums(abs(Y-x))))
  }
  
  if(type.measure=="mae" & family %in%  c("gaussian","binomial","poisson")){
    return((1 + (family=="binomial"))*mean(abs(x-y)))
  }
  
  if(type.measure=="class" & family=="binomial"){
    return(mean(y!=round(x)))
  }
  
  if(type.measure=="class" & family=="multinomial"){
    temp <- colnames(x)[apply(x,1,which.max)]
    return(mean(y!=temp))
  }
  
  if(type.measure=="auc" & family=="binomial"){
    weights <- table(foldid)
    cvraw <- rep(x=NA,times=length(weights))
    for(k in seq_along(weights)){
      cvraw[k] <- glmnet.auc(y=y[foldid==k],prob=x[foldid==k])
    }
    return(stats::weighted.mean(x=cvraw,w=weights,na.rm=TRUE))
  }
  
  stop(paste0("type measure \"",type.measure,"\" not implemented for family \"",family,"\"."))
  
}

#' @title 
#' glmnet:::auc
#' 
#' @description
#' Import of \code{\link[glmnet]{auc}} (internal function)
#' 
#' @param y
#' observed classes
#' 
#' @param prob
#' predicted probabilities
#' 
#' @param w
#' (ignored here)
#' 
#' @return
#' area under the ROC curve
#' 
#' @examples
#' NA
#' 
glmnet.auc <- get("auc",envir=asNamespace("glmnet"))

# same as cv.glmnet but with nzero (similar to dfmax and pmax)

#' @title 
#' glmnet::cv.glmnet
#' 
#' @description
#' Wrapper for \code{\link[glmnet]{cv.glmnet}},
#' with different handling of sparsity constraints.
#' 
#' @param ... 
#' see \code{\link[glmnet]{cv.glmnet}}
#' 
#' @param nzero
#' maximum number of non-zero coefficients\strong{:}
#' positive integer
#' 
#' @return
#' Object of class \code{\link[glmnet]{cv.glmnet}}.
#' 
#' @examples
#' NA
#' 
.cv.glmnet <- function(...,nzero){
  
  model <- glmnet::cv.glmnet(...)
  #cat(unique(model$nzero),"\n")
  
  for(i in seq_len(3)){ # original: 2
    # change from 1.05*min to min(1.05*)
    from <- 1.05*min(model$lambda[model$nzero==0],max(model$lambda),na.rm=TRUE) # original: 1.05*min
    if(is.na(nzero)||is.infinite(nzero)){ # unlimited model
      if(model$lambda.min==rev(model$lambda)[1]){
        to <- 0.01*min(model$lambda)
      } else {
        to <- max(model$lambda[model$lambda<model$lambda.min])
      }
    } else { # sparsity constraint
      if(all(model$nzero>=nzero)){
        to <- 0.01*min(model$lambda)
      } else {
        to <- min(model$lambda[model$nzero<=(nzero+1)],na.rm=TRUE)
        to <- max(model$lambda[model$lambda<model$lambda.min],to)
      }
    }
    
    lambda <- exp(seq(from=log(from),to=log(to),length.out=100))
    model <- glmnet::cv.glmnet(...,lambda=lambda)
    #cat(unique(model$nzero),"\n")
  }
  
  if(is.na(nzero)){return(model)}

  cond <- model$nzero<=nzero
  model$lambda <- model$lambda[cond]
  model$cvm <- model$cvm[cond]
  model$cvsd <- model$cvsd[cond]
  model$cvup <- model$cvup[cond]
  model$cvlo <- model$cvlo[cond]
  model$nzero <- model$nzero[cond]
  model$lambda.min <- model$lambda[which.min(model$cvm)]
  model$lambda.1se <- max(model$lambda[model$cvm<min(model$cvup)])
  model$glmnet.fit$a0 <- model$glmnet.fit$a0[cond]
  model$glmnet.fit$beta <- model$glmnet.fit$beta[,cond,drop=FALSE]
  model$glmnet.fit$df <- model$glmnet.fit$df[cond]
  model$glmnet.fit$lambda <- model$glmnet.fit$lambda[cond]
  model$glmnet.fit$dev.ratio <- model$glmnet.fit$dev.ratio[cond]
  model$glmnet.fit$call <- NULL
  model$glmnet.fit$dim[2] <- sum(cond)
  #cat(unique(model$nzero),"\n")
  
  return(model)
}


.glm <- function(y,X,family,intercept=TRUE,lower.limit=FALSE,upper.limit=FALSE,unit.sum=FALSE){
  
  alpha <- CVXR::Variable(1)
  beta <- CVXR::Variable(ncol(X))
  
  if(family=="gaussian"){
    objective <- CVXR::Minimize(mean((y-alpha-X%*%beta)^2))
  } else if(family=="binomial"){
    objective <- CVXR::Minimize(sum(alpha+X[y<=0,]%*%beta)+sum(CVXR::logistic(-alpha-X%*%beta)))
  } else if(family=="poisson"){
    stop("Family not implemented.")
    #objective <- CVXR::Minimize(mean(2*(y[y>=1]*log(y[y>=1]/(alpha+X[y>=1,]%*%beta)))-y[y>=1]+(alpha-X[y>=1,]%*%beta)))
  } else {
    stop("Family not implemented.")
  }
  
  constraints <- list(alpha==0,beta>=0,beta<=1,sum(beta)==1)[c(!intercept,lower.limit,upper.limit,unit.sum)]
  # either sum(beta)==1 or sum(beta)<=1
  problem <- CVXR::Problem(objective,constraints=constraints)
  result <- CVXR::solve(problem)

  alpha <- round(result$getValue(alpha),digits=6)
  beta <- round(result$getValue(beta),digits=6)
  return(list(alpha=alpha,beta=beta))
}


# This function should replace the one in palasso!
.folds <- function (y, nfolds, foldid = NULL) {
  if(nfolds > length(y)){nfolds <- length(y); warning("too many folds!")} # added
  if (!is.null(foldid)) {
    return(foldid)
  }
  if (survival::is.Surv(y)) {
    y <- y[,"status"]
  }
  if (all(y %in% c(0,1))) {
    foldid <- rep(x=NA,times=length(y))
    foldid[y==0] <- sample(x=rep(x=seq_len(nfolds),length.out=sum(y==0)))
    #foldid[y==1] <- sample(x=rep(x=seq(from=nfolds,to=1),length.out=sum(y==1))) # balanced size
    
    # balanced proportions
    max <- max(foldid,na.rm=TRUE)
    if(max==nfolds){
      seq <- seq_len(nfolds)
    } else {
      seq <- c(seq(from=max+1,to=nfolds,by=1),seq(from=1,to=max,by=1))
    }
    foldid[y==1] <- sample(x=rep(x=seq,length.out=sum(y==1)))
    
  } else {
    foldid <- sample(x=rep(x=seq_len(nfolds),length.out=length(y)))
  }
  return(foldid)
}

# This function is a copy from joinet!
.mean.function <- function(x,family){
  if(family %in% c("gaussian","cox")){
    return(x)
  } else if(family=="binomial"){
    return(1/(1+exp(-x)))
  } else if(family=="poisson"){
    return(exp(x))
  } else {
    stop("Family not implemented.",call.=FALSE)
  }
}



#y <- c(0,0,0,0,0,0,0,0,1,1)
#for(nfolds in 1:10){
#  folds <- .folds(y=y,nfolds=nfolds)
#  table <- table(folds,y)
#  print(table)
#}

#y <- c(0,0,0,0,0,0,0,0,1,1)
#fold <- rep(NA,times=length(y))
#fold[y==0] <- rep(1:min(sum(y==0),nfolds),length.out=sum(y==0))


### plot matrix
# #x <- matrix(c(-(12:1),1:12),nrow=4,ncol=6)
# plot.matrix <- function(x,cutoff=0,main="",xlab="col",ylab="row",xscale=NULL,yscale=NULL,xlas=1,ylas=3,...){
#   
#   if(is.null(xscale)){xscale <- seq_len(ncol(x))}
#   if(is.null(yscale)){yscale <- rev(seq_len(nrow(x)))}
#   
#   if(length(list(...))!=0){warning("Ignoring arguments.")}
#   
#   #k <- 100
#   #breaks <- stats::quantile(x,probs=seq(from=0,to=1,length.out=k+1))
#   
#   # start original 
#   #step <- 0.005
#   #min <- pmin(-0.05,min(x,na.rm=TRUE))
#   #max <- pmax(+0.05,max(x,na.rm=TRUE))
#   #breaks <- c(min,rev(seq(from=0,to=min,by=-step)),seq(from=step,max,by=step),max)
#   #k <- length(breaks)
#   
#   # start trial
#   step <- 0.01
#   min <- pmin(-1,min(x,na.rm=TRUE)) # was -1
#   max <- pmax(+1,max(x,na.rm=TRUE)) # was +1
#   breaks <- c(min,rev(seq(from=0,to=-0.5,by=-step)),seq(from=step,+0.5,by=step),max) # was -0.2 and +0.2
#   #k <- length(breaks)
#   
#   if("RColorBrewer" %in% .packages(all.available=TRUE)){
#     blue <- RColorBrewer::brewer.pal(n=9,name="Blues")
#     blue <- grDevices::colorRampPalette(colors=blue)(sum(breaks<cutoff)) # was k
#     #white <- ifelse(any(breaks==cutoff),"white",NULL)
#     red <- RColorBrewer::brewer.pal(n=9,name="Reds")
#     red <- grDevices::colorRampPalette(colors=red)(sum(breaks>cutoff)) # was k
#     col <- c(rev(blue[seq_len(sum(breaks<cutoff)-1)]),
#              rep("white",times=sum(breaks==cutoff)+1),
#              red[seq_len(sum(breaks>cutoff)-1)])
#   } else {
#     stop("install.packages(\"RColorBrewer\")")
#   }
#   
#   nrow <- nrow(x)
#   ncol <- ncol(x)
#   
#   graphics::plot.new()
#   graphics::par(xaxs="i",yaxs="i")
#   graphics::plot.window(xlim=c(1-0.5,ncol+0.5),ylim=c(1-0.5,nrow+0.5))
#   
#   graphics::title(main=main,line=0.5)
#   graphics::title(xlab=xlab,ylab=ylab,cex.lab=1,line=2.5)
#   graphics::image(x=seq_len(ncol),y=seq_len(nrow),z=t(x),breaks=breaks,col=col,add=TRUE) # t(x[rev(seq_len(nrow)),])
#   graphics::box()
#   
#   at <- unique(round(seq(from=1,to=ncol,by=ceiling(ncol/15))))
#   graphics::axis(side=1,at=seq_len(ncol)[at],labels=signif(xscale,digits=2)[at],las=xlas)
#   at <- unique(round(seq(from=1,to=nrow,by=ceiling(nrow/15))))
#   graphics::axis(side=2,at=seq_len(nrow)[at],labels=signif(yscale,digits=2)[at],las=ylas)
#   
# }

### simulate data
# .simulate <- function(n,p,rho,pi,family,factor){
#   if(pi==base::pi){stop("Invalid pi!")}
#   mu <- rep(x=0,times=p)
#   Sigma <- matrix(data=NA,nrow=p,ncol=p)
#   Sigma <- rho^abs(row(Sigma)-col(Sigma))
#   diag(Sigma) <- 1
#   X <- MASS::mvrnorm(n=n,mu=mu,Sigma=Sigma)
#   mean <- rowSums(X[,seq_len(pi*p),drop=FALSE])
#   # consider "rowMeans" or "scale(mean)"
#   if(family=="gaussian"){
#     y <- stats::rnorm(n=n,mean=mean,sd=factor*stats::sd(mean))
#   } else if(family=="binomial"){
#     y <- stats::rbinom(n=n,size=1,prob=1/(1+exp(-mean)))
#   } else if(family=="poisson"){
#     y <- stats::rpois(n=n,lambda=exp(mean))
#     y[is.na(y)] <- exp(20)
#   }
#   list <- list(y=y,X=X)
#   return(list)
# }

### test significant difference
# .test <- function (y,X){
#   fold <- palasso:::.folds(y=y,nfolds=5)
#   fold <- ifelse(fold==1,1,0)
#   fit <- starnet(y=y[fold==0],X=X[fold==0,],family="binomial")
#   pred <- predict.starnet(fit,newx=X[fold==1,])
#   if(any(pred < 0 | pred > 1)) {
#     stop("Outside unit interval.",call.=FALSE)
#   }
#   limit <- 1e-05
#   pred[pred < limit] <- limit
#   pred[pred > 1 - limit] <- 1 - limit
#   res <- -2*(y[fold==1]*log(pred)+(1-y[fold==1])*log(1 - pred))
#   old <- as.data.frame(res[,c("ridge","lasso","tune")])
#   new <- res[,"stack"]
#   pvalue <- apply(old,2,function(x) stats::wilcox.test(x=x,y=new,paired=TRUE,alternative="greater")$p.value)
#   return(pvalue)
# }

# Adapt code to family="cox", which has no intercept,
# i.e. coef.glmnet and coef.starnet do not yet work.
