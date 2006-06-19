#-->> BEGIN kendallATS aka cenken code

## Generics

setGeneric("cenken", 
           function(x, x.cen, y, y.cen) standardGeneric("cenken"))

## Classes 


## Private supporting routines

cenken0 =
function(x, x.cen, y, y.cen)
{
    kendallATS(x, x.cen, y, y.cen)
}

cenken1 =
function(x, x.cen, y, y.cen)
{
    xc = rep(FALSE, length(x))
    kendallATS(x, xc, x.cen, y)
}

# routine for cenken("numeric", "Cen")
cenken.unicen =
function(x, x.cen, y, y.cen)
{
    cenken1 (x, x.cen@Surv[,1], !x.cen@Surv[,2])
}

# routine for cenken("Cen", "Cen")
cenken.bicen =
function(x, x.cen, y, y.cen)
{
    kendallATS(x@Surv[,1], !x@Surv[,2], x.cen@Surv[,1], !x.cen@Surv[,2])
}

## Methods

setMethod("cenken",
          signature(x="numeric", x.cen="logical", y="numeric", y.cen="logical"),
          cenken0)

setMethod("cenken",
          signature(x="numeric", x.cen="missing", y="numeric", y.cen="logical"),
          cenken1)

setMethod("cenken",
          signature(x="numeric", x.cen="numeric", y="logical", y.cen="missing"),
          cenken1)

setMethod("cenken",
          signature(x="numeric", x.cen="Cen", y="missing", y.cen="missing"),
          cenken.unicen)

setMethod("cenken",
          signature(x="Cen", x.cen="Cen", y="missing", y.cen="missing"),
          cenken.bicen)

## kendallATS function -- the heart of the cenken routines
#
## Original code written by D Lorenz (USGS) for S-Plus
## Port to R by Lopaka Lee,  June 2006
## Requires ktau.f built as a native shared object

# Kendall's tau used as an estimate of the relation between x and y
#    with the ATS slope estimator (x and y left-censored).
#
# Coding History:
#    2005Jun30 DLLorenz original code modifed from Dennis Helsel Minitab macro
#    2005Jul14 DLLorenz date fix
#    2006Apr11 DLLorenz Modifed to print slope and medians
#    2006Jun05 DLLorenz Modifed bounds on slope estimates
#    2006Jun05     This version.
#

kendallATS = 
function(x, x.cen, y, y.cen) 
{

#kendallATS <- function(x, y, na.rm = T) {
  ## y must be of class lcens.
  ## x can be class lcens or a vector.
  ## Error checking.
  ## Make sure we have enough data.
#  if(!(class(x) == 'lcens'))
#    x <- lcens(x, rep(0, length(x)), interval=T)
#  ## Remove NAs if desired.
#  if(na.rm) {
#    OKs <- !is.na(y[,1]) & !is.na(x[,1])
#    y <- y[OKs,]
#    x <- x[OKs,]
#  }
#  else
  if(any(is.na(x) | is.na(y))) stop("Missing values not allowed")

  n <- length(x)
  if(length(x) < 3) stop("y and x should effectively be longer than 2")

  ## Median of the y values.
  median.y <- median(cenros(y, y.cen))
  if(is.na(median.y)) median.y <- min(y)

  ## Median of the x values.
  median.x <- median(cenros(x, x.cen))
  if(is.na(median.x)) median.x <- min(x)

  ## compute s, vars, tau, tau-b
  xx <- x
  cx <- x.cen
  yy <- y
  cy <- y.cen
  delx <- min(diff(sort(unique(xx)))) / 1000.
  dely <- min(diff(sort(unique(yy)))) / 1000.
  dupx <- xx - delx * cx
  diffx <- outer(dupx, dupx, "-")
  diffcx <- outer(cx, cx, "-")
  xplus <- outer(cx, -cx, "-")
  dupy <- yy - dely * cy
  diffy <- outer(dupy, dupy, "-")
  diffcy <- outer(cy, cy, "-")
  yplus <- outer(cy, -cy, "-")
  signyx <- sign(diffy * diffx)
  tt <- (sum(1-abs(sign(diffx)))-n)/2   # no. of pairwise ties in x
  uu <- (sum(1-abs(sign(diffy)))-n)/2   # no. of pairwise ties in y
  cix <- sign(diffcx)*sign(diffx)
  cix <- ifelse(cix <= 0, 0, 1)
  tt <- tt + sum(cix) / 2
  signyx <- signyx * (1 - cix)
  ciy <- sign(diffcy)*sign(diffy)
  ciy <- ifelse(ciy <= 0, 0, 1)
  uu <- uu + sum(ciy) / 2
  signyx <- signyx * (1 - ciy)
  xplus <- ifelse(xplus <= 1, 0, 1)
  yplus <- ifelse(yplus <= 1, 0, 1)
  diffx <- abs(sign(diffx))
  diffy <- abs(sign(diffy))
  tplus <- xplus*diffx
  uplus <- yplus*diffy
  tt <- tt + sum(tplus) / 2
  uu <- uu + sum(uplus) / 2
  itot <- sum(signyx * (1- xplus) * (1 - yplus))
  kenS <- itot/2 # kendall's S for original data
  tau <- (itot)/(n*(n-1))    # Kendall's tau-a for original data
  J <- n * (n - 1) / 2
  taub <- kenS/(sqrt(J - tt) * sqrt(J - uu))
  varS <- n*(n-1)*(2*n+5)/18   # var without tie correction
  ##  adjust varS for ties
  ##  taken from tiesboth.mac by Ed Gilroy
  ##  and modified for censored data by Dennis Helsel.
  ##  extracted from ckend.mac
  ## variance adjusted for ties between censored observations
  intg <- 1:n
  dupx <- xx - delx * cx # offset censored values by a small amount
  dupy <- yy - dely * cy
  dorder <- order(dupx)
  dxx <- dupx[dorder]
  dcx <- cx[dorder]
  dorder <- order(dupy)
  dyy <- dupy[dorder]
  dcy <- cy[dorder]
  tmpx <- dxx - intg * (1 - dcx) * delx
  tmpy <- dyy - intg * (1 - dcy) * dely
  rxlng <- rle(rank(tmpx))$lengths
  nrxlng <- table(rxlng)
  rxlng <- as.integer(names(nrxlng))
  x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
  x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
  x3 <- nrxlng * rxlng * (rxlng - 1)
  rylng <- rle(rank(tmpy))$lengths
  nrylng <- table(rylng)
  rylng <- as.integer(names(nrylng))
  y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
  y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
  y3 <- nrylng * rylng * (rylng - 1)
  delc <- (sum(x1) + sum(y1)) / 18 -
    sum(x2) * sum(y2) / (9 * n * (n - 1) * (n - 2)) -
      sum(x3) * sum(y3) / ( 2 * n * (n - 1))
  ## correction for C-C pairs double counted as UC-C pairs
  ## made to deluc
  x4 <- nrxlng * (rxlng - 1)
  y4 <- nrylng * (rylng - 1)
  ## variance adjusted for ties in censored and uncensored obs.
  tmpx <- intg * dcx - 1
  tmpx <- ifelse(tmpx < 0, 0, tmpx)
  nrxlng <- sum(tmpx)
  rxlng <- 2
  x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
  x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
  x3 <- nrxlng * rxlng * (rxlng - 1)
  tmpy <- intg * dcy - 1
  tmpy <- ifelse(tmpy < 0, 0, tmpy)
  nrylng <- sum(tmpy)
  rylng <- 2
  y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
  y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
  y3 <- nrylng * rylng * (rylng - 1)
  deluc <- (sum(x1) + sum(y1)) / 18 -
    sum(x2) * sum(y2) / (9 * n * (n - 1) * (n - 2)) -
      sum(x3) * sum(y3) / ( 2 * n * (n - 1)) - (sum(x4) + sum(y4))
  ## variance adjustment for ties between uncensored obs.
  dxx <- dxx - intg * dcx * delx
  dyy <- dyy - intg * dcy * dely
  rxlng <- rle(rank(dxx))$lengths
  nrxlng <- table(rxlng)
  rxlng <- as.integer(names(nrxlng))
  x1 <- nrxlng * rxlng * (rxlng - 1) * (2 * rxlng + 5)
  x2 <- nrxlng * rxlng * (rxlng - 1) * (rxlng - 2)
  x3 <- nrxlng * rxlng * (rxlng - 1)
  rylng <- rle(rank(dyy))$lengths
  nrylng <- table(rylng)
  rylng <- as.integer(names(nrylng))
  y1 <- nrylng * rylng * (rylng - 1) * (2 * rylng + 5)
  y2 <- nrylng * rylng * (rylng - 1) * (rylng - 2)
  y3 <- nrylng * rylng * (rylng - 1)
  delu <- (sum(x1) + sum(y1)) / 18 -
    sum(x2) * sum(y2) / (9 * n * (n - 1) * (n - 2)) -
      sum(x3) * sum(y3) / ( 2 * n * (n - 1))
  varS <- varS - delc - deluc - delu
  p.val <- 2 * (1 - pnorm((abs(kenS - sign(kenS)))/sqrt(varS)))
  ## prepare for ATS slope estimate
  flipx <- max(xx) - xx + 1
  cx <- 1 - cx
  flipy <- max(yy) - yy + 1
  cy <- 1 - cy
  slopes <- unlist(lapply(1:10, function(i, y, x)
                          ((y[i] - y[1:i]) / (x[i] - x[1:i])),
                          yy * cy, xx * cx)) # substitute zeros for censored values
  bounds <- range(slopes, na.rm=T)
  iter <- 1000
  result <- .Fortran("ktau",
                     x = as.double(flipx), y = as.double(flipy),
                     icx = as.integer(cx), icy = as.integer(cy),
                     num = as.integer(n), slope = double(1),
                     lbound = as.double(bounds[1]),
                     ubound = as.double(bounds[2]),
                     dev = double(1), iter = as.integer(iter))
  sen.slope <- result$slope
  if(result$iter >= iter)
    warning("Slope estimation did not converge")
  ## A line representing the trend of the data then is given by
                                        #
  ##    y(t) = sen.slope*(t-med.time)+med.data
                                        #
  ##    This line has the slope of the trend and passes through
  ##       the point (t=med.time, y=med.data)
  ## Compute the coefficients for the line
  coef <- c(sen.slope*(-median.x)+median.y, sen.slope) # intercept and slope
  ## Return the statistics.
  names(taub) <- "tau"
  zero <- 0
  names(zero) <- "slope"
  est <- c(sen.slope, median.y, median.x)
  names(est) <- c("slope", "median.y", "median.x")
  method <- "Kendall's tau with the ATS slope estimator"
  z <- list(method = method, statistic = taub, p.value = p.val,
            estimate = est, alternative="two.sided",
            coef=coef, slope.se = result$dev)
  z$data.name <-  paste(deparse(substitute(x)), "and",deparse(substitute(y)))
  oldClass(z) <- "htest"
  return(z)
}

#-->> END kendallATS aka cenken code
