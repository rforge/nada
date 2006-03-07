#-->> BEGIN Regression on Maximum Likelihood Estimation (MLE) section

## Generics

setGeneric("cenmle",
  function(obs, censored, groups, ...) standardGeneric("cenmle"))

## Classes

setOldClass("survreg")

setClass("cenmle", representation(survreg="survreg"))

## Methods

setMethod("cenmle",
          signature(obs="formula", censored="missing", groups="missing"),
                    function(obs, censored, groups, ...)
{
    obj = survreg(asSurv(obs), ...)
    new("cenmle", survreg=obj)
})

setMethod("cenmle", 
          signature(obs="Cen", censored="missing", groups="missing"), 
          cencen.Cen)

setMethod("cenmle",
          signature(obs="numeric", censored="logical", groups="missing"),
                    cencen.vectors)

setMethod("cenmle", 
          signature(obs="numeric", censored="logical", groups="factor"), 
          cencen.vectors.groups)


setMethod("print", signature(x="cenmle"), function(x, ...)
{
    ret = c(mean(x), median(x), sd(x))
    names(ret) = c("mean", "median", "sd")
    print(ret)
    invisible(ret)
})

setMethod("transform", signature(x="cenmle"), function(x, ...)
{
    switch (x@survreg$dist, 
        lognormal = transform_cenmle_lognormal(x),
        x
    )
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

## Private utility functions

#transform_cenmle_default = function(x) return(x)

transform_cenmle_gaussian =
function(x) 
{
    ret = x;
    return(ret);
}

transform_cenmle_lognormal =
function(x) 
{
    ret = x;
    ret@survreg$coef  = exp(ret@survreg$coef)
    ret@survreg$scale = exp(ret@survreg$scale)
    return(ret);
}

#-->> END Regression on Maximum Likelihood Estimation (MLE) section

