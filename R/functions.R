
#--- Main function -------------------------------------------------------------

#' @export
#' @title
#' Stacked Elastic Net Regression
#' 
#' @description
#' Implements stacked elastic net regression.
#'  
#' @param y
#' numeric response\strong{:}
#' vector of length \eqn{n}
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
#' vector of values between \eqn{0} (ridge) and \eqn{1} (lasso);
#' or \code{NULL} (equidistance)
#'
#' @param nfolds
#' number of folds
#'
#' @param foldid
#' fold identifiers\strong{:}
#' vector with entries between \eqn{1} and \code{nfolds};
#' or \code{NULL} (balance)
#' 
#' @param type.measure
#' loss function\strong{:}
#' character "deviance", "class", "mse" or "mae"
#' (see \code{\link[glmnet]{cv.glmnet}})
#'
#' @param alpha.meta
#' elastic net mixing parameters for stacking\strong{:}
#' vector of values between \eqn{0} (ridge) and \eqn{1} (lasso),
#' \emph{see details}
#' 
#' @param intercept,upper.limit,unit.sum
#' settings for meta-learner\strong{:} logical
#' 
#' @param penalty.factor
#' penalty factors for base-learners
#' 
#' @param grouped
#' logical
#' 
#' @param ...
#' further arguments passed to \code{\link[glmnet]{glmnet}}
#' 
#' @references 
#' A Rauschenberger, E Glaab, and MA van de Wiel (2020)
#' "Improved elastic net regression through stacked generalisation"
#' \emph{Manuscript in preparation.}
#' 
#' @details
#' Combine predictions from \emph{some} \code{alpha} with \code{alpha.meta}\eqn{=1},
#' or from \emph{all} \code{alpha} with \code{alpha.meta}\eqn{=0}.
#' We recommend to use \code{alpha.meta}\eqn{=0} (default) for stability.
#' 
#' @examples
#' set.seed(1)
#' n <- 30; p <- 50
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet::starnet(y=y,X=X,family="gaussian")
#' 
starnet <- function(y,X,family="gaussian",nalpha=21,alpha=NULL,nfolds=10,foldid=NULL,type.measure="deviance",alpha.meta=1,grouped=TRUE,penalty.factor=NULL,intercept=NULL,upper.limit=NULL,unit.sum=NULL,...){
  
  #--- temporary ---
  # family <- "gaussian"; nalpha <- 21; alpha <- NULL; nfolds <- 10; foldid <- NULL; type.measure <- "deviance"; alpha.meta <- 0; grouped <- TRUE; penalty.factor <- NULL; intercept <- TRUE; upper.limit=TRUE; unit.sum=FALSE
  
  #--- default ---
  if(all(is.null(intercept),is.null(upper.limit),is.null(unit.sum))){
    if(is.na(alpha.meta)){
      intercept <- FALSE
      upper.limit <- unit.sum <-TRUE
    } else if(alpha.meta==1){
      intercept <- TRUE
      upper.limit <- unit.sum <- FALSE
    } else if(alpha.meta==0){
      intercept <- upper.limit <- TRUE
      unit.sum <- FALSE
    }
  }
  
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
  
  #if(nfolds > length(y)){nfolds <- length(y); warning("Using LOOCV!")}
  
  if(family=="cox" & grouped){stop("Implement grouped!",call.=FALSE)}
  
  #--- fold identifiers ---
  if(is.null(foldid)){
    foldid <- .folds(y=y,nfolds=nfolds) # was palasso:::.folds 2020-02-19
  } # else {
    nfolds <- length(unique(foldid)) # outside "if" since 2020-02-17
  #}
    
  #message("trial - unbalanced folds")
  #foldid <- palasso:::.folds(y=y,nfolds=nfolds) # remove this
  #nfolds <- length(unique(foldid)) # remove this
  
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
    fit <- joinet:::.mean.function(link[[i]],family=family)
    base[[i]]$cvm <- apply(fit,2,function(x) .loss(y=y,x=x,family=family,type.measure=type.measure,foldid=foldid,grouped=grouped))
    base[[i]]$id.min <- which.min(base[[i]]$cvm)
    base[[i]]$lambda.min <- base[[i]]$lambda[base[[i]]$id.min]
  }
  
  #--- predictions ---
  hat <- matrix(NA,nrow=nrow(X),ncol=nalpha)
  for(i in seq_len(nalpha)){
    hat[,i] <- link[[i]][,base[[i]]$id.min]
  }
  
  meta <- list()
  
  #message("reshuffle folds for meta-learner")
  #foldid <- .folds(y=y,nfolds=nfolds) # delete this
  
  #--- meta cross-validation ---
  for(i in seq_along(alpha.meta)){
    #message("trial - force parameters")
    #intercept <- TRUE; upper.limit <- FALSE; unit.sum <- FALSE # remove this
    
    if(is.na(alpha.meta[i])){next}
    
    #message("trial: without base-lasso")
    #penalty.factor=ifelse(alpha==1,Inf,1) # REMOVE THIS
    
    meta[[i]] <- glmnet::cv.glmnet(y=y,x=hat,
                                     lower.limits=0,
                                     upper.limits=ifelse(upper.limit,1,Inf), # was 1
                                     foldid=foldid,
                                     family=family,
                                     type.measure=type.measure,
                                     intercept=intercept,
                                     alpha=alpha.meta[i],
                                     grouped=grouped,...)
    
    if(unit.sum){
      cond <- Matrix::colSums(meta[[i]]$glmnet.fit$beta)>1
      meta[[i]]$cvm[cond] <- Inf
      meta[[i]]$lambda.min <- meta[[i]]$lambda[which.min(meta[[i]]$cvm)]
      #sum <- sum(coef(meta[[i]],s="lambda.min")[-1])
      #message("sum = ",round(sum,2))
    }

  }
  
  #browser()
  #plot(sapply(base,function(x) min(x$cvm)))
  #as.numeric(coef(meta[[1]],s="lambda.min")[-1])

  #--- convex combination ---  
  for(i in seq_along(alpha.meta)){
    #message("trial - force parameters")
    #intercept <- TRUE; upper.limit <- TRUE; unit.sum <- TRUE # remove this
    
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
      #link <- glm$alpha + hat %*% glm$beta
      y_hat <- joinet:::.mean.function(x=link,family=family)
      cvm <- .loss(y=y,x=y_hat,family=family,type.measure=type.measure,foldid=foldid,grouped=grouped)
    } else {
      cvm <- 0
    }
    
    #if(round(cvm,digits=2)==0){browser()}
    
    meta[[i]]$cvm <- c(Inf,cvm)
    meta[[i]]$lambda.min <- 0
    class(meta[[i]]) <- "cv.glmnet"
    class(meta[[i]]$glmnet.fit) <- "glmnet"
  }
  
  names(meta) <- paste0("alpha",alpha.meta)
  
  #--- tune meta alpha ---
  cvm_meta <- sapply(meta,function(x) min(x$cvm))
  message(paste0(names(cvm_meta)," ",round(cvm_meta,digits=2)," "))
  id_meta <- which.min(cvm_meta)
  #message(names(id_meta))
  meta <- meta[[names(id_meta)]]
  message(paste(paste(round(stats::coef(meta,s="lambda.min")[1],digits=3)),"_",
                paste(round(stats::coef(meta,s="lambda.min")[-1],digits=3),collapse=" ")))
  #plot(coef(meta,s="lambda.min")[-1]) 
  #message(paste(round(glm$alpha,digits=3),"_",paste(round(glm$beta,digits=3),
  #      collapse=" ")))
  
  if(FALSE){
    # message("trial - weights from cvm")
    # cvm <- sapply(base,function(x) min(x$cvm))
    # prior <- 1-(cvm-min(cvm))/(max(cvm)-min(cvm))
    # weight <- prior/sum(prior)
    # meta$glmnet.fit$a0 <- rep(x=0,times=2)
    # meta$glmnet.fit$beta <- matrix(data=rep(weight,times=2),ncol=2)
  }
  
  if(FALSE){
    # # penalty factors from cvm
    # cvm <- sapply(base,function(x) min(x$cvm))
    # prior <- 1-(cvm-min(cvm))/(max(cvm)-min(cvm))
    # weight <- prior/sum(prior)
    # #penalty.factor=1/weight, # DELETE THIS!
    # ##penalty.factor <- ifelse(alpha==1,Inf,1) # exclude lasso
    # # ALSO DELETE PENALTY.FACTOR BELOW !
  }
  
  if(FALSE){
    ## fused lasso
    #if(family=="gaussian"){
    #  model <- "linear"
    #} else if(family=="binomial"){
    #  model <- "logistic"
    #} else if(family=="poisson"){
    #  model <- "poisson"
    #} else {
    #  stop("Not implemented.")
    #}
    #optL1 <- penalized::optL1(response=y,penalized=hat,positive=TRUE,
    #               fusedl=TRUE,fold=foldid,lambda2=0,model=model)
    #optL1$lambda
    #optL2 <- penalized::optL2(response=y,penalized=hat,positive=TRUE,
    #              fusedl=TRUE,fold=foldid,lambda1=optL1$lambda,model=model)
    #optL2$lambda
  }
  
  if(FALSE){
    #message("trial - convex combination with CVXR")
    # x <- hat
    # p <- nalpha
    # # different code for linear and logistic regression!
    # 
    # alpha <- CVXR::Variable(1)
    # beta <- CVXR::Variable(p)
    # if(family=="gaussian"){
    #   objective <- CVXR::Minimize(mean((y-alpha-x%*%beta)^2))
    # } else if(family=="binomial"){
    #   objective <- CVXR::Minimize(sum(alpha+x[y<=0,]%*%beta)+sum(CVXR::logistic(-alpha-x%*%beta)))
    # } else {
    #   stop("Not implemented!")
    # }
    # #constraints <- list() # for testing purposes
    # constraints <- list(alpha==0,beta>=0,sum(beta)==1)
    # problem <- CVXR::Problem(objective,constraints=constraints)
    # result <- CVXR::solve(problem)
    # #result$getValue(alpha)
    # #result$getValue(beta)
    # 
    # if(family=="gaussian"){
    #   glm <- stats::glm(y~x,family=gaussian)
    # } else if(family=="binomial"){
    #   glm <- stats::glm(y~x,family=binomial)
    # } else {
    #   stop("Not implemented!")
    # }
    # #glm$coefficients
    # 
    # beta <- round(result$getValue(beta),digits=6)
    
    # Consider different coding (-1/1 vs 0/1)! See
    # https://cvxr.rbind.io/cvxr_examples/cvxr_logistic-regression/
    
    # temporary solution
    #glm <- .glm(y=y,X=hat,family=family)# trial
    #meta$glmnet.fit$a0[] <- glm$alpha
    #meta$glmnet.fit$beta[,] <- glm$beta
    #message(paste(round(glm$beta,3),collapse=" "))
    
  }
  
  ## trial start 
  ## Plot loss and coefs.
  #plot(sapply(base,function(x) x$cvm[x$id.min]))
  #plot(as.numeric(cor(y,hat)))
  #plot(coef(meta[[1]],s=0))
  ## trial start
  
  ## trial start ## added on 2020-02-18
  #for(i in seq_along(alpha.meta)){
  #  meta[[i]]$glmnet.fit$beta <- apply(meta[[i]]$glmnet.fit$beta,2,function(x) x/sum(x))
  #}
  ## Note: this should be inside cv.glmnet, not outside!
  
  if(FALSE){ # 2020-02-18
    # message("trial - convex combination")
    # if(length(alpha.meta)!=1){stop("Implement multiple alpha!")}
    # # full fit
    # meta <- list()
    # meta$glmnet.fit <- glmnet::glmnet(y=y,x=hat,lower.limits=0,
    #                                   family=family,alpha=alpha.meta,intercept=FALSE,...)
    # meta$glmnet.fit$beta <- apply(meta$glmnet.fit$beta,2,function(x) x/sum(x))
    # meta$glmnet.fit$beta[is.na(meta$glmnet.fit$beta)] <- 0
    # meta$lambda <- meta$glmnet.fit$lambda
    # link <- matrix(data=NA,nrow=n,ncol=length(meta$lambda))
    # # meta cv
    # for(k in seq_len(nfolds)){
    #   y0 <- y[foldid!=k]
    #   y1 <- y[foldid==k]
    #   hat0 <- hat[foldid!=k,,drop=FALSE]
    #   hat1 <- hat[foldid==k,,drop=FALSE]
    #   object <- glmnet::glmnet(y=y0,x=hat0,family=family,alpha=alpha.meta,intercept=FALSE,...)
    #   object$beta <- apply(object$beta,2,function(x) x/sum(x))
    #   object$beta[is.na(object$beta)] <- 0
    #   temp <- stats::predict(object=object,newx=hat1,type="link",
    #                            s=meta$glmnet.fit$lambda)
    #   link[foldid==k,seq_len(ncol(temp))] <- temp
    # }
    # # tune meta lambda
    # fit <- joinet:::.mean.function(link,family=family)
    # meta$cvm <- apply(fit,2,function(x) .loss(y=y,x=x,family=family,type.measure=type.measure,foldid=foldid,grouped=grouped))
    # meta$id.min <- which.min(meta$cvm)
    # meta$lambda.min <- meta$lambda[meta$id.min]
    # class(meta) <- "cv.glmnet"
  }
  
  if(FALSE){
    # # meta-learner without intercept might work better
    # meta <- list()
    # for(i in seq_along(alpha.meta)){
    #   message("temporary - no upper limits")
    #   meta[[i]] <- glmnet::cv.glmnet(y=y,x=hat,
    #                                  lower.limits=0,
    #                                  upper.limits=1,
    #                                  foldid=foldid,
    #                                  family=family,
    #                                  type.measure=type.measure,
    #                                  alpha=alpha.meta[i],
    #                                  intercept=FALSE,
    #                                  grouped=grouped,...)
    # }
    # names(meta) <- paste0("alpha",alpha.meta)
  }
  
  if(FALSE){ # 2020-02-19
    # message("trial - simple convex combination")
    # if(length(alpha.meta)!=1){stop("Implement multiple alpha!")}
    # lambda <- seq(from=0.1,to=0,by=-0.01)
    # meta <- list()
    # meta$glmnet.fit <- glmnet::glmnet(y=y,x=hat,lower.limits=0,
    #                                 family=family,alpha=alpha.meta,
    #                                 lambda=lambda,intercept=FALSE,...)
    # meta$glmnet.fit$beta <- apply(meta$glmnet.fit$beta,2,function(x) x/sum(x))
    # meta$glmnet.fit$beta[is.na(meta$glmnet.fit$beta)] <- 0
    # meta$lambda <- meta$glmnet.fit$lambda
    # meta$id.min <- length(lambda)
    # meta$lambda.min <- 0
    # class(meta) <- "cv.glmnet"
  }
  
  if(FALSE){ # 2020-02-19
    # message("debugging - to be removed")
    # # Here stacking returns same output as tuning.
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
#' number of non-zero coefficients\strong{:}
#' positive integer, or \code{NULL}
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @examples
#' set.seed(1)
#' n <- 30; p <- 50
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' y_hat <- predict(object,newx=X[c(1),,drop=FALSE])
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
  #frame$none <- as.numeric(stats::predict(object=x$meta,newx=link,s=Inf,type="link")) # original
  frame$none <- as.numeric(stats::predict(object=x$base$alpha1$glmnet.fit,
                  newx=newx,s=Inf,type="link")) # trial 2020-02-20
  
  if(type=="link"){
    return(frame)
  } else  if(type=="response"){
    #frame <- sapply(frame,function(y) joinet:::.mean.function(y,family=x$info$family))
    frame <- lapply(frame,function(y) joinet:::.mean.function(y,family=x$info$family))
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
#' @param object
#' \link[starnet]{starnet} object
#' 
#' @param nzero
#' number of non-zero coefficients\strong{:}
#' positive integer, or \code{NULL}
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @examples
#' set.seed(1)
#' n <- 30; p <- 50
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' coef <- coef(object)
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
  
  # posthoc selection
  if(!is.null(nzero)){
     eta <- as.numeric(pool$alpha + object$info$X %*% pool$beta)
     if(stats::sd(eta)==0){return(list(alpha=pool$alpha,beta=0*pool$beta))}
     model <- .cv.glmnet(y=eta,x=object$info$X,family="gaussian",
                 alpha=1,penalty.factor=1/abs(pool$beta),
                 foldid=object$info$foldid,nzero=nzero)
     coef <- stats::coef(model,s="lambda.min")
     return(list(alpha=coef[1],beta=coef[-1]))
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
#' @examples
#' set.seed(1)
#' n <- 30; p <- 50
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' weights(object)
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
#' @examples
#' set.seed(1)
#' n <- 30; p <- 20
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' object <- starnet(y=y,X=X)
#' print(object)
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
#' or \code{NULL} (no posthoc feature selection)
#' 
#' @param grouped
#' logical
#' 
#' @param nfolds.ext,nfolds.int,foldid.ext,foldid.int
#' (number of) external/internal folds
#' 
#' @param ...
#' further arguments (not applicable)
#' 
#' @examples
#' set.seed(1)
#' n <- 50; p <- 20
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' y <- rnorm(n=n,mean=rowSums(X[,1:20]))
#' loss <- cv.starnet(y=y,X=X)
#' 
cv.starnet <- function(y,X,family="gaussian",nalpha=21,alpha=NULL,nfolds.ext=10,nfolds.int=10,foldid.ext=NULL,foldid.int=NULL,type.measure="deviance",alpha.meta=1,grouped=TRUE,nzero=NULL,intercept=NULL,upper.limit=NULL,unit.sum=NULL,...){
  
  # family <- "gaussian"; nfolds.ext <- 5; nfolds.int <- 10; foldid.ext <- NULL; foldid.int <- NULL; type.measure <- "deviance"; alpha.meta <- 0; alpha = NULL; nalpha <- 21; grouped <- TRUE; nzero <- NULL; intercept <- upper.limit <- unit.sum <- NULL
  
  #--- fold identifiers ---
  if(is.null(foldid.ext)){
    #if(nfolds>=length(y)){ # add this ? 2020-02-17
    #  fold <- sample(x = rep(x = seq_len(nfolds), length.out = length(y)))
    #} else {
      fold <- .folds(y=y,nfolds=nfolds.ext) # was palasso:::.folds 2020-02-19
    #}
  } else {
    fold <- foldid.ext
    #nfolds <- length(unique(foldid)) # removed 2020-02-17
  }
  #nfolds.ext <- length(unique(fold)) # added 2020-02-17
  nfolds.ext <- max(fold) # added 2020-02-20
  # Maybe use max(fold) or length(unique(fold,na.rm=TRUE)),
  # such that we can use giant testing sets.
  
  #--- alpha sequence ---
  if(is.null(alpha)){
    alpha <- seq(from=0,to=1,length.out=nalpha)
  } else {
    alpha <- alpha
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
      #if(i<6){next} # remove this!
      # dense models
      fit <- starnet(y=y[fold!=i],X=X[fold!=i,],alpha=alpha,nalpha=nalpha,nfolds=nfolds.int,foldid=foldid.int[fold!=i],family=family,type.measure=type.measure,alpha.meta=alpha.meta,grouped=grouped,intercept=intercept,upper.limit=upper.limit,unit.sum=unit.sum,...) # INSERT ,...
      temp <- predict.starnet(fit,newx=X[fold==i,,drop=FALSE]) # added drop=FALSE 2020-02-17
      #if(sum(fold==i)==1){temp <- as.data.frame(t(temp))} # trial 2020-02-17
      for(j in c(meta,base)){
        pred[fold==i,j] <- unlist(temp[[j]])
      }
      # sparse models
      if(!is.null(nzero)){
        for(j in seq_along(nzero)){
          #if(i==6){browser()} # remove this!
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
  
  loss <- apply(pred[fold!=0,],2,function(x) .loss(y=y[fold!=0],x=x,family=family,type.measure=type.measure,foldid=fold,grouped=grouped))
 # inserted "fold!=0" (twice)
  
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


#' @title
#' Simulation
#'
#' @description
#' Prints object of class \link[starnet]{starnet}.
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
#' character \code{"sparse"}, \code{"dense"} or \code{"blocks"}
#' 
#' @param family
#' \code{"gaussian"}, \code{"binomial"} or \code{"poisson"}
#' 
#' @examples
#' NA
#' 
.simulate <- function(n,p,mode,family="gaussian"){
  q <- 3
  Z <- matrix(data=stats::rnorm(n*q),nrow=n,ncol=q)
  y <- rowSums(Z)
  if(family=="binomial"){y <- round(1/(1+exp(-y)))}
  if(family=="poisson"){y <- round(exp(y))}
  X <- matrix(data=stats::rnorm(n*p),nrow=n,ncol=p)
  if(mode=="sparse"){
    X[,1] <- sqrt(0.1)*X[,1] + sqrt(0.9)*Z[,1]
    X[,p] <- sqrt(0.1)*X[,p] + sqrt(0.9)*Z[,2]
  } else if(mode=="dense"){
    X[,1:250] <- sqrt(0.9)*X[,1:250] + sqrt(0.1)*Z[,1]
    X[,(p-250):p] <- sqrt(0.9)*X[,(p-250):p] + sqrt(0.1)*Z[,2]
  } else if(mode=="blocks"){
    X[,1:25] <- sqrt(0.5)*X[,1:25] + sqrt(0.5)*Z[,1]
    X[,(p-25):p] <- sqrt(0.5)*X[,(p-25):p] + sqrt(0.5)*Z[,2]
  }
  return(list(y=y,X=X))
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
#' @examples
#' NA
#' 
.cv.glmnet <- function(...,nzero){
  
  model <- glmnet::cv.glmnet(...)
  #model <- glmnet::cv.glmnet(y=y,x=X,family="binomial")
  #model <- glmnet::cv.glmnet(y=y,x=X,family="binomial",lambda=lambda)
  
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
        to <- max(model$lambda[model$lambda<model$lambda.min],to) # trial
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
  
  # maybe use intercept=TRUE|FALSE|
  
  #if(intercept){
  #  message("trial - force intercept")
  #  if(family=="gaussian"){alpha <- mean(y)}
  #  if(family=="binomial"){alpha <- log(mean(y)/(1-mean(y)));message(coef(glm(y~1,family="binomial"))[1])}
  #} else {
    alpha <- CVXR::Variable(1)
  #}
  
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
  
  #if(convex.combination){
  #  if(intercept){
  #    constraints <- list(beta>=0,sum(beta)<=1) # note: "<=" not "=="
  #    #message("temporary: no convex combination!")
  #    #constraints <- list(beta>=0,beta<=1)
  #  } else {
  #    constraints <- list(alpha==0,beta>=0,sum(beta)<=1)
  #    #message("temporary: no convex combination!")
  #    #constraints <- list(alpha==0,beta>=0,beta<=1)
  #  }
  #} else {
  #  constraints <- list()
  #}
  
  constraints <- list(alpha==0,beta>=0,beta<=1,sum(beta)==1)[c(!intercept,lower.limit,upper.limit,unit.sum)]
  
  #message("trial - force sum to one") # alternative: sum(beta)<=1
  #message("trial - less or equal to one") # alternative: sum(beta)==1
  
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
  if (all(y %in% c(0, 1))) {
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
    
    foldid <- sample(x = rep(x = seq_len(nfolds), length.out = length(y)))
  }
  return(foldid)
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
