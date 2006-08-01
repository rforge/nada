\name{cenreg}
\alias{cenreg}
\title{
    Compute regression equations and likelihood correlation coefficient
for censored data.
}
\description{
    Computes regression equations for singly censored data using maximum
    likelihood estimation. Estimates of slopes and intercept, tests
    for significance of parameters,and predicted quantiles (Median =
    points on the line) with confidence intervals can be computed.
}
\synopsis{cenreg(obs, censored, groups)}
\usage{
    cenreg(obs, censored, groups)
}
\arguments{
\item{obs}{
    Either a numeric vector of observations or a formula.
    See examples below.
} 
\item{censored}{
    If a formula is not specified, this should be
    a logical vector indicating TRUE where an observation in obs is
    censored (a less-than value) and FALSE otherwise.  
}
\item{groups}{
    If a formula is not specified, this should be a numeric
    or factor vector that represents the explanatory variable.
}
}
\details{
    If you are using the formula interface:  The \code{censored} and
    \code{groups} parameters are not specified -- all information is
    provided via a formula as the \code{obs} parameter.  The formula
    must have a \code{Cen} object as the response on the left of the
    \code{~} operator and, if desired, terms separated by + operators
    on the right. See examples below.

    The likelihood r correlation coefficient measures the linear
    association between y (groups) and x (obs), based on the difference
    in log likelihoods between the fitted model and the null model.
    Slopes and intercepts are fit by maximum likelihood.  A lognormal
    distribution is fit by default, with a normal distribution being an
    option.  Estimates of predicted values on the line can be obtained
    by specifying the values for all x variables at which y is to be
    predicted.  Requesting the median (p=0.5) will provide estimates on
    the line for a lognormal distribution.  Estimates of the mean are
    also possible, as are estimates of other percentiles.  Equations for
    confidence intervals follow those of Meeker and Escobar (1098).
}
\value{
    Returns a \code{summary.cenreg} object.
}
\references{
    Helsel, Dennis R. (2005).  Nondectects and Data Analysis; Statistics for
    censored environmental data.  John Wiley and Sons, USA, NJ.

    Meeker, W.Q. and L. A. Escobar (1998).  Statistical Methods for
    Reliability Data.  John Wiley and Sons, USA, NJ.
}
\examples{

    # (example on pg 213 of the NADA book)
    data(TCEReg)
    with(TCEReg, cenreg(log(TCEConc), TCECen, PopDensity))

    # Using the formula interface
    with(TCEReg, cenreg(Cen(log(TCEConc), TCECen)~PopDensity))

    # A model with 2 or more explanatory variables requires the formula interface
    with(TCEReg, cenreg(Cen(log(TCEConc), TCECen)~PopDensity+Depth))
}
