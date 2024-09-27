
#' @name starnet-package
#' @keywords documentation
#' @docType package
#' 
#' @aliases starnet-package
#' 
#' @title
#' Stacked Elastic Net Regression
#' 
#' @description
#' The R package \code{starnet} implements stacked elastic net regression.
#' The elastic net generalises ridge and lasso regularisation.
#' Instead of fixing or tuning the mixing parameter alpha,
#' we combine multiple alphas by stacked generalisation.
#' 
#' @details
#' Use function \code{\link{starnet}} for model fitting.
#' Type \code{library(starnet)} and then \code{?starnet} or
#' \code{help("starnet)"} to open its help file.
#' 
#' See the vignette for further examples.
#' Type \code{vignette("starnet")} or \code{browseVignettes("starnet")}
#' to open the vignette.
#' 
#' @references
#' \href{https://orcid.org/0000-0001-6498-4801}{Armin Rauschenberger},
#' \href{https://orcid.org/0000-0003-3977-7469}{Enrico Glaab},
#' and \href{https://orcid.org/0000-0003-4780-8472}{Mark A. van de Wiel}
#' (2021).
#' "Predictive and interpretable models via the stacked elastic net".
#' \emph{Bioinformatics} 37(14):2012-2016.
#' \doi{10.1093/bioinformatics/btaa535}.
#' (Click
#' \href{https://orbilu.uni.lu/bitstream/10993/43221/1/starnet.pdf}{here}
#' to access PDF.)
#' 
#' @examples
#' \dontshow{
#' if(!grepl('SunOS',Sys.info()['sysname'])){
#' #--- data simulation ---
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' # n samples, p features
#' 
#' #--- model fitting ---
#' object <- starnet(y=y,X=X)
#' # "base": one model for each alpha
#' # "meta": model for stacking them 
#' 
#' #--- make predictions ---
#' y_hat <- predict(object,newx=X)
#' # one column for each alpha,
#' # and for tuning and stacking
#' 
#' #--- extract coefficients ---
#' coef <- coef(object)
#' # scalar "alpha": intercept
#' # vector "beta": slopes}}
#' 
#' \donttest{
#' #--- data simulation ---
#' n <- 50; p <- 100
#' y <- rnorm(n=n)
#' X <- matrix(rnorm(n*p),nrow=n,ncol=p)
#' # n samples, p features
#' 
#' #--- model fitting ---
#' object <- starnet(y=y,X=X)
#' # "base": one model for each alpha
#' # "meta": model for stacking them 
#' 
#' #--- make predictions ---
#' y_hat <- predict(object,newx=X)
#' # one column for each alpha,
#' # and for tuning and stacking
#' 
#' #--- extract coefficients ---
#' coef <- coef(object)
#' # scalar "alpha": intercept
#' # vector "beta": slopes
#' 
#' #--- model comparison ---
#' loss <- cv.starnet(y=y,X=X)
#' # cross-validated loss for different alpha,
#' # and for tuning and stacking}}
#' 
"_PACKAGE"

