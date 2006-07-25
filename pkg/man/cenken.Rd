\name{cenken}
\alias{cenken}
\title{
    Replace: Compute an ECDF for Censored Data
}
\description{
    Replace: Computes an estimate of an empirical cumulative 
    distribution function (ECDF) for censored data using the 
    Kaplan-Meier method.
}
\synopsis{cenken(x, xcen, y, ycen)}
\usage{
    cenken(x, xcen, y, ycen)
    cenken(x, y, ycen)
}
\arguments{
\item{x}{
    A numeric vector of observations.
} 
\item{xcen}{
    A  logical vector indicating TRUE where an observation in x is
    censored (a less-than value) and FALSE otherwise.  Can be 
    missing/omitted for the case where x is not censored.
}
\item{y}{
    A numeric vector of observations.
}
\item{ycen}{
    A logical vector indicating TRUE where an observation in y is
    censored (a less-than value) and FALSE otherwise.  
}
%\item{formula}{
%    The formula must have a \code{Cen} object as the response on the
%    left of the \code{~} operator and, if desired, terms separated by +
%    operators on the right.  One of the terms may be a \code{strata}
%    object.  For a single survival curve the \code{"~ 1"} part of the
%    formula is not required.
%}

}
\details{
    Computational details and other discussions of the method go here.
}
\value{
    Returns tau (Kendall's tau), slope, and p-value for the regression.
}
\references{
    Helsel, Dennis R. (2005).  Nondectects and Data Analysis; Statistics for
    censored environmental data.  John Wiley and Sons, USA, NJ.
}
%\author{
%    Lopaka(Rob) Lee <rclee@usgs.gov>
%
%    Dennis Helsel <dhelsel@usgs.gov>
%} 
%\seealso{
%    \code{\link{cenken-methods}},
%}
\examples{
    # Both x and y are censored
    data(Golden)
    with(Golden, cenken(Blood, BloodCen, Kidney, KidneyCen))

    # x is not censored
    # (example on pg 213 of the NADA book)
    data(TCEReg)
    with(TCEReg, cenken(PopDensity, log(TCEConc), TCECen))

%    # x (Year) is not censored -- formula interface
%    data(DFe)
%    with(DFe, cenken(Cen(Summer, SummerCen)~Year)
}
\keyword{survival}
