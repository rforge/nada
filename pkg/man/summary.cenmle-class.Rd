\name{summary.cenmle-class}
\docType{class}
\alias{summary.cenmle-class}
\alias{summary.cenmle}

\title{Class "summary.cenmle"}
\description{
    A "summary.cenmle" object is returned from \code{summary}.
}
\section{Objects from the Class}{
    Objects can be created by calls of the form 
    \code{summary(cenmle(obs, censored))}.
}
\section{Slots}{
  \describe{
    No slots, it just extends a list.
    %\item{\code{survreg}:}{Object of class \code{"survreg"}  }
  }
}
\section{Extends}{
    Class \code{"list"}.
    Class \code{"vector"}, by class \code{"list"}.
}
\section{Methods}{
  \describe{
    %\item{summary}{\code{signature(object = "cenmle")}: ... }
  }
}
\author{
    Lopaka(Rob) Lee <rclee@usgs.gov>
} 
\seealso{
  \code{\link{cenmle}}
}
\examples{
    obs      = c(0.5,    0.5,   1.0,  1.5,   5.0,    10,   100)
    censored = c(TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)

    summary(cenmle(obs, censored))
}
\keyword{classes}
