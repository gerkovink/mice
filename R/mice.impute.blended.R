#' Imputation by blended predictive mean matching
#'
#' @aliases mice.impute.blended blended
#' @param y Vector to be imputed
#' @param ry Logical vector of length \code{length(y)} indicating the
#' the subset \code{y[ry]} of elements in \code{y} to which the imputation
#' model is fitted. The \code{ry} generally distinguishes the observed
#' (\code{TRUE}) and missing values (\code{FALSE}) in \code{y}.
#' @param x Numeric design matrix with \code{length(y)} rows with predictors for
#' \code{y}. Matrix \code{x} may have no missing values.
#' @param wy Logical vector of length \code{length(y)}. A \code{TRUE} value
#' indicates locations in \code{y} for which imputations are created.
#' @param blend The blending factor where a value 1 would result in PMM and a 
#' value 0 would only use the Mahalanobis distance
#' @param donors The size of the donor pool among which a draw is made.
#' The default is \code{donors = 5L}. Setting \code{donors = 1L} always selects
#' the closest match, but is not recommended. Values between 3L and 10L
#' provide the best results in most cases (Morris et al, 2015).
#' @param matchtype Type of matching distance. The default choice
#' (\code{matchtype = 1L}) calculates the distance between
#' the \emph{predicted} value of \code{yobs} and
#' the \emph{drawn} values of \code{ymis} (called type-1 matching).
#' Other choices are \code{matchtype = 0L}
#' (distance between predicted values) and \code{matchtype = 2L}
#' (distance between drawn values).
#' @param rank Logical flag if the distance must be calculated by rank
#' @param ridge The ridge penalty used in \code{.norm.draw()} to prevent
#' problems with multicollinearity. The default is \code{ridge = 1e-05},
#' which means that 0.01 percent of the diagonal is added to the cross-product.
#' Larger ridges may result in more biased estimates. For highly noisy data
#' (e.g. many junk variables), set \code{ridge = 1e-06} or even lower to
#' reduce bias. For highly collinear data, set \code{ridge = 1e-04} or higher.
#' \code{2.22} (June 2014) to \code{3.11.7} (Oct 2020). Since version \code{3.12.0}
#' \code{mice()} uses the much faster \code{matchindex} C function. Use
#' the deprecated \code{matcher} function only for exact reproduction.
#' @param \dots Other named arguments.
#' @return Vector with imputed data, same type as \code{y}, and of length
#' \code{sum(wy)}
#' @author Anaïs Fopma, Mingyang Cai, Gerko Vink
#' @details
#' Imputation by blended pmm
#' 
#' @family univariate imputation functions
#' @keywords datagen
#' @export
mice.impute.blended <- function(y, ry, x, wy = NULL,
                                blend, 
                                donors = 5L,
                                matchtype = 1L, 
                                ridge = 1e-05, 
                                rank = TRUE,
                                ...) {
  
  if (is.null(wy)) {
    wy <- !ry
  }
  x <- cbind(1, as.matrix(x))
  ynum <- y
  if (is.factor(y)) {
    ynum <- as.integer(y)
  }
  parm <- .norm.draw(ynum, ry, x, ridge = ridge)
  
  if (matchtype == 0L) {
    yhatobs <- x[ry, , drop = FALSE] %*% parm$coef
    yhatmis <- x[wy, , drop = FALSE] %*% parm$coef
  }
  if (matchtype == 1L) {
    yhatobs <- x[ry, , drop = FALSE] %*% parm$coef
    yhatmis <- x[wy, , drop = FALSE] %*% parm$beta
  }
  if (matchtype == 2L) {
    yhatobs <- x[ry, , drop = FALSE] %*% parm$beta
    yhatmis <- x[wy, , drop = FALSE] %*% parm$beta
  }
  
  Sx <- cov(x[ry, -1 , drop = FALSE])
  maha.dis <- apply(x[wy, -1 , drop = FALSE], 1, 
                    function(z) mahalanobis(x[ry, -1 , drop = FALSE], z, Sx))
  
  predict.match <- function(z){
    d <- abs(yhatobs - z)
    f <- d > 0
    a1 <- ifelse(any(f), min(d[f]), 1)
    d <- d + runif(length(d), 0, a1 / 10^10)
    return(d)
  }
  
  predict.dis <- map(yhatmis, predict.match) %>% 
    unlist %>% 
    matrix(nrow = length(yhatobs), 
           ncol = length(yhatmis))
  
  if (rank){ 
    pd.rank <- apply(predict.dis, 2, function(x) rank(x, ties.method  = "random"))
    mh.rank <- apply(maha.dis, 2, function(x) rank(x, ties.method  = "random"))
    blend.dis <- blend * pd.rank + (1 - blend) * mh.rank
  } else {
    blend.dis <- blend * scale(predict.dis) + (1 - blend) * scale(maha.dis)
  }
  
  donor.select <- function(z){
    if (donors == 1) {
      return(y[which.min(z)])
    }
    donors <- min(donors, length(z))
    donors <- max(donors, 1)
    y.obs<-y[ry]
    ds <- sort.int(z, partial = donors) 
    m <- sample(y.obs[z <= ds[donors]], 1)
  }
  return(apply(blend.dis, 2, donor.select))
}  

