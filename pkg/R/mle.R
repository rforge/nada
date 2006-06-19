#-->> BEGIN Maximum Likelihood Estimation (MLE) Regression section

## Generics

setGeneric("cenmle",
  function(obs, censored, groups, ...) standardGeneric("cenmle"))

## Classes

setOldClass("survreg")

setClass("cenmle", representation(survreg="survreg"))
setClass("cenmle-gaussian", representation("cenmle"))
setClass("cenmle-lognormal", representation("cenmle"))
setClass("summary.cenmle", representation("list"))

## Core Methods

# This keeps things consistent with the survival package
cenreg = cenmle

setMethod("cenmle",
          signature(obs="formula", censored="missing", groups="missing"),
                    function(obs, censored, groups, dist, ...)
{
    dist = ifelse(missing(dist), "lognormal", dist)
    switch(dist,
        gaussian  = new_cenmle_gaussian(obs, dist, ...),
        lognormal = new_cenmle_lognormal(obs, dist, ...),
        survreg(asSurv(obs), dist=dist, ...)
    )
})

# Use the cencen.* framework (in All.R) to dispatch from these cases
setMethod("cenmle", 
          signature(obs="Cen", censored="missing", groups="missing"), 
          cencen.Cen)

setMethod("cenmle",
          signature(obs="numeric", censored="logical", groups="missing"),
                    cencen.vectors)

setMethod("cenmle", 
          signature(obs="numeric", censored="logical", groups="factor"), 
          cencen.vectors.groups)

setMethod("cenmle", 
          signature(obs="numeric", censored="logical", groups="numeric"), 
          cencen.vectors.groups)

setMethod("summary", signature(object="cenmle"), function(object, ...)
{
    s = unclass(summary(object@survreg, ...))
    s$r = rloglik(s)
    return (new("summary.cenmle", s))
})


setMethod("print", signature(x="cenmle"), function(x, ...)
{
    coef    = x@survreg$coefficients

    median  = median(x)
    mean    = mean(x)
    sd      = sd(x)

    ret = matrix(c(median(x), mean(x), sd(x)), ncol=3, nrow=length(coef),
                 dimnames=list(names(coef), c("median", "mean", "sd")))

    print(ret)
    invisible(ret)
})


setMethod("predict", signature(object="cenmle"), summary) 

#setMethod("predict", signature(object="cenmle"), 
#          function(object, newdata, conf.int=FALSE, ...)
#{
#    predict(object@survreg, newdata, ...)
#})

setMethod("residuals", signature(object="cenmle"), function(object, ...)
{
    residuals(object@survreg, ...)
})

setMethod("coef", signature(object="cenmle"), function(object, ...)
{
    coef(object@survreg, ...)
})

setMethod("median", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    as.vector(exp(x@survreg$coef))
})

setMethod("sd", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    ret = exp(2*x@survreg$coef + x@survreg$scale^2)*(exp(x@survreg$scale^2)-1)
    as.vector(sqrt(ret))
})

setMethod("mean", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    as.vector(exp(x@survreg$coef + 0.5*(x@survreg$scale)^2))
})

setMethod("median", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    as.vector(x@survreg$coef)
})

setMethod("sd", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    as.vector(x@survreg$scale)
})

setMethod("mean", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    as.vector(x@survreg$coef)
})

setMethod("cor", signature(x="cenmle"), function(x, y, use, method)
{
    rloglik(summary(x))
})

# This is the hideous summary.survreg from the survival package --
# for now we hack it to do what we want. 
setMethod("print", signature(x="summary.cenmle"), function(x, ...)
{
    digits = max(options()$digits - 4, 3)
    n <- x$n

    print(x$table, digits = digits)
    if (nrow(x$var)==length(x$coefficients)) 
      {
	    cat("\nScale fixed at", format(x$scale, digits=digits),"\n") 
      }
    else if (length(x$scale)==1) 
      {
	    cat ("\nScale=", format(x$scale, digits=digits), "\n")
      }
    else 
      {
	    cat("\nScale:\n")
	    print(x$scale, digits=digits, ...)
	  }

    cat("\n", x$parms, "\n", sep='')
    df  <- sum(x$df) - x$idf   # The sum is for penalized models
    cat("Loglik(model)=", format(round(x$loglik[2],1)),
	"  Loglik(intercept only)=", format(round(x$loglik[1],1)), "\n")

    cat("Loglik-r: ", x$r, "\n");

    if (df <= 0) cat ("\n")
    else
      {
	    cat("\nChisq=", format(round(x$chi,2)), "on", round(df,1),
		"degrees of freedom, p=", 
		format(signif(1-pchisq(x$chi, df),2)), "\n")
      }

    if (x$robust) cat("(Loglikelihood assumes independent observations)\n")

    cat("Number of Newton-Raphson Iterations:", format(trunc(x$iter)), "\n")

    if (!length(x$na.action)) cat("n =", x$n, "\n")
	else cat("n =", x$n, " (", naprint(x$na.action), ")\n", sep="")

    if(!is.null(x$correl)) 
      {
        p <- dim(x$correl)[2]
        if(p > 1) 
          {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(x$correl)
            x$correl[ll] <- format(round(x$correl[ll], digits=digits))
            x$correl[!ll] <- ""
            print(x$correl[-1,  - p, drop = FALSE], quote = FALSE)
          }
      }
    cat("\n")
    invisible(NULL)
})

## Supporting Functions -- private

# Compute the log-likelihood correlation coef (r) from a cenmle summary obj.
# Eq 11.4 pg 187 of Dennis' book
rloglik =
function(x)
{
    n = x$n
    G = -2 * diff(x$loglik)

    sqrt(1 - exp(G/n))
}

# cenmle for lognormal distributions

new_cenmle_lognormal =
function(formula, dist, ...)
{
    new("cenmle-lognormal", survreg=survreg(asSurv(formula), dist=dist, ...))
}

# cenmle for gaussian, or normal, distributions

# If a normal distribution is assumed the input data must be expressed
# as an interval between zero and the DL.  They cannot simply be stated as
# 'left' censored, because that allows some probability of going below 0.
# Since environmental/analytical data are _usually_ not negative,
# estimates will be biased low and wrong.  So with the normal option
# and left censoring, internally we must use interval censoring.  The end of
# the interval are the detected values.  The start of the interval will
# have identical numbers in it for the detects, and a 0 for the 
# nondetects (a simple trick is: start = obs - obs * censored).

new_cenmle_gaussian =
function(formula, dist, ...)
{
    obs      = formula[[2]][[2]]
    censored = formula[[2]][[3]]

    ## This assumes a single-level grouping -- fix it!
    groups   = formula[[3]]

    f = as.formula(substitute(Surv(o - o * c, o, type="interval2")~g, 
                                   list(o=obs, c=censored, g=groups)))
    environment(f) = parent.frame()
    new("cenmle-gaussian", survreg=survreg(f, dist=dist, ...))
}

#-->> END Regression on Maximum Likelihood Estimation (MLE) section

