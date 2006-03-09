#-->> BEGIN Regression on Maximum Likelihood Estimation (MLE) section

## Generics

setGeneric("cenmle",
  function(obs, censored, groups, ...) standardGeneric("cenmle"))

## Classes

setOldClass("survreg")

setClass("cenmle", representation(survreg="survreg"))
setClass("cenmle-gaussian", representation("cenmle"))
setClass("cenmle-lognormal", representation("cenmle"))

## Core Methods

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

setMethod("summary", signature(object="cenmle"), function(object, ...)
{
    summary(object@survreg, ...)
})

setMethod("predict", signature(object="cenmle"), function(object, newdata, ...)
{
    predict(object@survreg, newdata, ...)
})

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
    # To do: remove NAs?
    as.vector(exp(x@survreg$coef))
})

setMethod("sd", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    ret = exp(2*x@survreg$coef + x@survreg$scale^2)*(exp(x@survreg$scale^2)-1)
    as.vector(sqrt(ret))
})

setMethod("mean", signature(x="cenmle-lognormal"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(exp(x@survreg$coef + 0.5*(x@survreg$scale)^2))
})

setMethod("median", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(x@survreg$coef)
})

setMethod("sd", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(x@survreg$scale)
})

setMethod("mean", signature(x="cenmle-gaussian"), function(x, na.rm=FALSE)
{
    # To do: remove NAs?
    as.vector(x@survreg$coef)
})


## Supporting Functions 

new_cenmle_lognormal =
function(formula, dist, ...)
{
    new("cenmle-lognormal", survreg=survreg(asSurv(formula), dist=dist, ...))
}

new_cenmle_gaussian =
function(formula, dist, ...)
{
    obs      = eval.parent(formula[[2]][[2]])
    censored = eval.parent(formula[[2]][[3]])
    groups   = eval.parent(formula[[3]])

    start = obs - obs * censored
    end = obs

    f = as.formula(substitute(Surv(s, e, type="interval2")~g, 
                                   list(s=start, e=end, g=groups)))
    environment(f) = parent.frame()
    new("cenmle-gaussian", survreg=survreg(f, dist=dist, ...))
}

#-->> END Regression on Maximum Likelihood Estimation (MLE) section

