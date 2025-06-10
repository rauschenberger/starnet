
if(!grepl('SunOS',Sys.info()['sysname'])){
for(family in c("gaussian","binomial","poisson")){
  n <- 50; p <- 100 # was 30 and 50
  if(family=="cox"){
    list <- list()
    list$y <- survival::Surv(time=abs(rnorm(n)),
                             event=rbinom(n=n,size=1,prob=0.2))
    list$X <- matrix(rnorm(n*p),nrow=n,ncol=p)
  } else {
      list <- .simulate.block(n=n,p=p,mode="sparse",family=family)
  }
  
  foldid <- .folds(y=list$y,nfolds=5)
  
  glmnet <- glmnet::cv.glmnet(y=list$y,x=list$X,family=family,foldid=foldid,alpha=0.5)
  object <- starnet(y=list$y,X=list$X,family=family,foldid=foldid)
  
  #--- Equality glmnet and starnet ---
  
  reprod <- object$base$alpha0.5

  testthat::test_that("lambda: glmnet = starnet",{
    a <- glmnet$lambda
    b <- reprod$lambda
    seq <- seq_len(min(length(a),length(b)))
    cond <- all(a[seq]==b[seq])
    testthat::expect_true(cond)
  })
  
  testthat::test_that("lambda.min: glmnet = starnet",{
    a <- glmnet$lambda.min
    b <- reprod$lambda.min
    cond <- a==b
    testthat::expect_true(cond)
  })
  
  testthat::test_that("cvm: glmnet = starnet",{
    a <- glmnet$cvm
    b <- reprod$cvm
    seq <- seq_len(min(length(a),length(b)))
    cond <- all(abs(a[seq]-b[seq]) < 1e-06)
    testthat::expect_true(cond)
  })
  
  testthat::test_that("predict: glmnet = starnet",{
    if(family=="cox"){
      a <- glmnet::predict.coxnet(glmnet$glmnet.fit,newx=list$X,type="response")
      b <- glmnet::predict.coxnet(reprod$glmnet.fit,newx=list$X,type="response")
    } else {
      a <- glmnet::predict.glmnet(glmnet$glmnet.fit,newx=list$X,type="response")
      b <- glmnet::predict.glmnet(reprod$glmnet.fit,newx=list$X,type="response")
    }
    cond <- all(a==b)
    testthat::expect_true(cond)
  })
  
  #--- Equivalence stacking and pooling ---
  
  testthat::test_that("stacking = pooling",{
      pred <- predict.starnet(object=object,newx=list$X,type="response")
      a <- pred$stack
  
      coef <- coef.starnet(object=object)
      b <- .mean.function(coef$alpha + list$X %*% coef$beta,family=family)
  
      pred <- predict.starnet(object=object,newx=list$X,type="link")
      weights <- weights.starnet(object=object)
      sub <- pred[,grepl(pattern="alpha",x=colnames(pred))]
      c <- .mean.function(weights[1] + as.matrix(sub) %*% weights[-1],family=family)
  
      cond <- all(abs((a-b)<1e-06) & abs((a-c)<1e-06) & abs((b-c)<1e-06))
      
      testthat::expect_true(cond)
  })
  
}

#--- testing the convex combination ---

if("CVXR" %in% .packages(all.available=TRUE)){
  for(family in c("gaussian","binomial")){
    
    n <- 100; p <- 5
    list <- .simulate.block(n=n,p=p,mode="sparse")
    if(family=="binomial"){
      list$y <- stats::rbinom(n=n,size=1,prob=1/(1+exp(-list$y)))
    } else if(family=="poisson"){
      list$y <- stats::rpois(n=n,lambda=exp(list$y))
    }
    
    glm0 <- .glm(y=list$y,X=list$X,family=family)
    
    if(family=="gaussian"){
      glm1 <- stats::glm(y~X,family=gaussian,data=list)
    } else if(family=="binomial"){
      glm1 <- stats::glm(y~X,family=binomial,data=list)
    } else if(family=="poisson") {
      glm1 <- stats::glm(y~X,family=poisson,data=list)
    }
    
    testthat::test_that("same coefficients",{
      change <- (c(glm0$alpha,glm0$beta)-glm1$coefficients)/glm1$coefficients
      testthat::expect_true(all(change<0.01))
    })
    
  }
}

#--- testing the loss function ---

set.seed(1)
for(family in c("gaussian","binomial","poisson","mgaussian","multinomial","cox")){
  for(type.measure in c("deviance","mse","mae","class","auc")){
    
    if(type.measure=="class" & !family %in% c("binomial","multinomial")){next}
    if(type.measure=="auc" & family!="binomial"){next}
    if(family=="cox" & type.measure!="deviance"){next}
    
    if(family=="binomial" & type.measure=="auc"){warning("AUC requires signal!"); next}

    # data simulation
    n <- 200; p <- 150
    alpha <- 0 # runif(1)
    
    if(family=="gaussian"){
      y <- stats::rnorm(n=n)
    } else if(family=="binomial"){
      y <- stats::rbinom(n=n,size=1,prob=0.2)
    } else if(family=="poisson"){
      y <- stats::rpois(n=n,lambda=4)
    } else if(family=="multinomial"){
      y <- sample(c("a","b",1,2),replace=TRUE,size=n)
    } else if(family=="mgaussian"){
      y <- cbind(stats::rnorm(n=n),stats::rnorm(n=n),stats::rnorm(n=n))
    } else if(family=="cox"){
      y <- survival::Surv(time=rpois(n=n,lambda=4)+1,
                          event=rbinom(n=n,size=1,prob=0.8))
    }
    
    X <- matrix(data=stats::rnorm(n=n*p),nrow=n,ncol=p)
    nfolds <- 10
    foldid <- sample(seq_len(nfolds),replace=TRUE,size=n)
    
    if(family %in% c("multinomial","mgaussian")){
      if(family=="multinomial"){dim <- length(unique(y))}
      if(family=="mgaussian"){dim <- ncol(y)}
      y_hat <- matrix(NA,nrow=n,ncol=dim)
      if(family=="multinomial"){colnames(y_hat) <- sort(unique(y))}
    } else {
      y_hat <- rep(NA,times=n)
    }
    
    lambda <- min(glmnet::glmnet(y=y,x=X,alpha=alpha,family=family)$lambda)
    
    beta <- pred <- list()
    for(k in seq_len(nfolds)){
      if(family=="mgaussian"){
        y0 <- y[foldid!=k,]
      } else {
        y0 <- y[foldid!=k]
      }
      X0 <- X[foldid!=k,,drop=FALSE]
      X1 <- X[foldid==k,,drop=FALSE]
      glmnet <- glmnet::glmnet(y=y0,x=X0,alpha=alpha,lambda=lambda,family=family)
      if(family %in% c("multinomial","mgaussian")){
          y_hat[foldid==k,] <- stats::predict(object=glmnet,newx=X1,type="response")
      } else {
          temp <- stats::predict(object=glmnet,newx=X1,type="link")
          y_hat[foldid==k] <- .mean.function(temp,family=family)
      }
      beta[[k]] <- coef(glmnet)
      pred[[k]] <- predict(glmnet,newx=X,type="link")
    }
    
    grouped <- FALSE # family!="cox"
    # Implement loss function for "cox-grouped"!
    # See notes below! See palasso implementation!
    a <- glmnet::cv.glmnet(y=y,x=X,alpha=alpha,lambda=c(lambda,0.5*lambda),foldid=foldid,family=family,type.measure=type.measure,grouped=grouped)$cvm[1]
    b <- .loss(y=y,x=y_hat,family=family,type.measure=type.measure,foldid=foldid,grouped=grouped)
    
    if(family=="cox"){next} # due to update in glmnet
    
    testthat::test_that("same loss",{
      testthat::expect_true(abs(a-b)<1e-03)
    })
    
  }
}

if(FALSE){
  # Run code from above before running this code!
  
  ### Cox: ungrouped deviance
  old <- new <- cvraw <- rep(NA,times=10)
  for(i in 1:10){
    old[i] <- glmnet::coxnet.deviance(x=X[foldid==i,],y=y[foldid==i],beta=beta[[i]])
    new[i] <- glmnet::coxnet.deviance(pred=pred[[i]][foldid==i],y=y[foldid==i])
  }
  all(abs(old-new)<1e-06)
  
  ### Cox: grouped deviance
  old <- new <- rep(NA,times=10)
  for(i in 1:10){
    full <- glmnet::coxnet.deviance(x=X,y=y,beta=beta[[i]])
    mink <- glmnet::coxnet.deviance(x=X[foldid!=i,],y=y[foldid!=i],beta=beta[[i]])
    old[i] <- full-mink
    full <- glmnet::coxnet.deviance(pred=pred[[i]],y=y)
    mink <- glmnet::coxnet.deviance(pred=pred[[i]][foldid!=i],y=y[foldid!=i])
    new[i] <- full-mink
  }
  all(abs(old-new)<1e-06)
  
  ### Cox: averaging
   # choose one from above
  # old
  cvraw <- old
  status = y[, "status"]
  weights = as.vector(tapply(status, foldid, sum))
  temp = as.matrix(cvraw/weights,ncol=1)
  apply(temp, 2, weighted.mean, w = weights, na.rm = TRUE)
  # new
  cvraw <- new
  weights <- tapply(X = y[, "status"], INDEX = foldid, FUN = sum)
  stats::weighted.mean(x = cvraw/weights, w = weights, na.rm = TRUE)
  
  grouped <- TRUE # FALSE or TRUE
  glmnet::cv.glmnet(y=y,x=X,alpha=alpha,lambda=c(lambda,0.5*lambda),foldid=foldid,family="cox",grouped=grouped)$cvm[1]
  
}
}
