\name{cenreg}
\alias{cenreg}
\title{
    Compute regression equations and likelihood correlation coefficient
for censored data.
}
\description{
    Computes regression equations for singly (y only) censored data using
    maximum likelihood estimation. Estimates of slopes and intercept, 
    tests for significance of parameters,and predicted quantiles 
    (Median = points on the line) with confidence intervals can be computed.
}
\synopsis{cenreg(x, y, ycen)}
\usage{
    cenreg(x, y, ycen)
}
\arguments{
\item{x}{
    A numeric vector of observations.
} 
\item{y}{
    A numeric vector of observations.
}
\item{ycen}{
    A logical vector indicating TRUE where an observation in y is
    censored (a less-than value) and FALSE otherwise.  
}
}
\details{
    The likelihood r correlation coefficient measures the linear
    association between y and x, based on the difference in log
    likelihoods between the fitted model and the null model.  Slopes
    and intercepts are fit by maximum likelihood.  A lognormal 
    distribution is fit by default, with a normal distribution being
    an option.  Estimates of predicted values on the line can be
    obtained by specifying the values for all x variables at which
    y is to be predicted.  Requesting the median (p=0.5) will provide
    estimates on the line for a lognormal distribution.  Estimates of
    the mean are also possible, as are estimates of other percentiles.
    Equations for confidence intervals follow those of Meeker and Escobar
    (1098).
}
\value{
    Returns intercept, slope(s), and p-value(s) for the regression.
}
\references{
    Helsel, Dennis R. (2005).  Nondectects and Data Analysis; Statistics for
    censored environmental data.  John Wiley and Sons, USA, NJ.

    Meeker, W.Q. and L. A. Escobar (1998).  Statistical Methods for
    Reliability Data.  John Wiley and Sons, USA, NJ.
}
\examples{
    # Both x and y are censored
    data(Golden)
    with(Golden, cenken(Blood, BloodCen, Kidney, KidneyCen))

    # x is not censored
    # (example on pg 213 of the NADA book)
    data(TCEReg)
    with(TCEReg, cenken(PopDensity, log(TCEConc), TCECen))
}
