#-->> BEGIN Regression on Maximum Likelihood Estimation (MLE) section

## Generics

setGeneric("cenmle",
  function(obs, censored, groups, ...) standardGeneric("cenmle"))

## Classes

setOldClass("survreg")

setClass("cenmle", representation(survreg="survreg"))

## Utility functions -- for use with Methods 

## Begin cenmle.* functions
# These routines are allow cenmle
# to be used with Cen objects, vectors, or formulas as input.
# Note they all convert input to a formula and call the formula method.

cenmle.Cen =
function(obs, censored, groups, ...)
{
    cl = match.call()
    f = as.formula(substitute(a~1, list(a=cl[[2]])))
    environment(f) = parent.frame()
    callGeneric(f, ...)
}

cenmle.vectors =
function(obs, censored, groups, ...)
{
    cl = match.call()
    f = as.formula(substitute(Cen(a, b)~1, list(a=cl[[2]], b=cl[[3]])))
    environment(f) = parent.frame()
    callGeneric(f, ...)
}

cenmle.vectors.groups =
function(obs, censored, groups, ...)
{
    cl = match.call()
    f = substitute(Cen(a, b)~g, list(a=cl[[2]], b=cl[[3]], g=cl[[4]]))
    f = as.formula(f)
    environment(f) = parent.frame()
    callGeneric(f, ...)
}

## End cenmle.* routines


## Core Methods

# cenmle for formulas
setMethod("cenmle",
          signature(obs="formula", censored="missing", groups="missing"),
                    function(obs, censored, groups, ...)
{
    obj = survreg(asSurv(obs), ...)
    new("cenmle", survreg=obj)
})

setMethod("cenmle",
          signature(obs="numeric", censored="logical", groups="missing"),
                    cenmle.vectors)

setMethod("print", signature(x="cenmle"), function(x, ...)
{
    ret = c(mean(x), median(x), sd(x))
    names(ret) = c("mean", "median", "sd")
    print(ret)
    invisible(ret)
})

setMethod("median", signature(x="cenmle"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(exp(x@survreg$coef))
})

setMethod("sd", signature(x="cenmle"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    ret = exp(2*x@survreg$coef + x@survreg$scale^2)*(exp(x@survreg$scale^2)-1)
    as.vector(sqrt(ret))
})

setMethod("mean", signature(x="cenmle"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(exp(x@survreg$coef + 0.5*(x@survreg$scale)^2))
})

#-->> END Regression on Maximum Likelihood Estimation (MLE) section

