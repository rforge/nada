###
# NADA for R by Lopaka Lee.
#
# Version 1.3-0
# Copyright (2004, 2005, 2006) Lopaka Lee
#
# A S-language software module based on 
# methodologies described by Dennis R. Helsel in his book 
# Nondetects and Data Analysis: Statistics for Censored Environmental Data.
#
# NADA is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
# 
# NADA is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.  You should have received a copy of the GNU General
# Public License along with NADA; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
###

### Definitions common to All sections in this package

## Globals

.onLoad = function(lib, pkg) 
{
    require(methods)
    library.dynam("NADA", pkg, lib)
}

NADAprobs = c(0.05,0.10,0.25,0.50,0.75,0.90,0.95)

## Generics

setGeneric("plot", function(x, y, ...) standardGeneric("plot"))

setGeneric("print", function(x, ...) standardGeneric("print"))

setGeneric("summary", function(object, ...) standardGeneric("summary"))

setGeneric("mean", function(x, ...) standardGeneric("mean"))

setGeneric("sd", function(x, na.rm=FALSE) standardGeneric("sd"))

setGeneric("median", function(x, na.rm=FALSE) standardGeneric("median"))

setGeneric("min", function(..., na.rm=FALSE) standardGeneric("min"))
setGeneric("max", function(..., na.rm=FALSE) standardGeneric("max"))

setGeneric("quantile", function(x, ...) standardGeneric("quantile"))

setGeneric("predict", function(object, ...) standardGeneric("predict"))

setGeneric("pexceed", function(object, ...) standardGeneric("pexceed"))

setGeneric("lines", function(x, ...) standardGeneric("lines"))

setGeneric("boxplot", function(x, ...) standardGeneric("boxplot"))

setGeneric("cor", function(x, y = NULL, use = "all.obs",
          method = c("pearson", "kendall", "spearman")) standardGeneric("cor"))

## Broken for the time being
#setGeneric("abline", 
#           function(a, b, h, v, reg, coef, untf, col, lty, lwd, ...) 
#           standardGeneric("abline"))

setGeneric("residuals", function(object, ...) standardGeneric("residuals"))

setGeneric("coef", function(object, ...) standardGeneric("coef"))

#if (as.numeric(version$minor) < 3) {
#    if (!isGeneric("transform"))
#      setGeneric("transform", function(x, ...) standardGeneric("transform"))
#} else {
#    if (!isGeneric("transform"))
#      setGeneric("transform", function(`_data`, ...)
#                 standardGeneric("transform"))
#}


## Classes

setClass("NADAlist", "list")

## Methods

setMethod("print", signature("NADAlist"), function(x, ...)
{
    tag = names(x)
    for (i in 1:length(x))
      {
        cat(tag[i], "\n")
        print(x[[i]])
        cat("\n")
      }
})

#-->> BEGIN general utility functions

##
# split_qual extracts qualifed and unqualifed vectors from a vector
# containing concatenated qualifiying characters and numeric values
# like "<0.5".  Only handles one kind of censoring character/symbol.
split_qual =
function(v, qual.symbol = "<")
{
    v = as.character(v)

    obs = as.numeric(sub(qual.symbol, "", v))
    cen = rep(FALSE, length(obs))
    cen[grep(qual.symbol, v)] = TRUE 

    return(list(obs=obs, cen=cen))
}
splitQual = split_qual

## pct_cen -- Simple function to save some typing
pct_cen =
function(obs, censored)
{
    if (!is.logical(censored)) 
      {
        stop("censored indicator must be logical vector!\n")
      }

    return(100*(length(obs[censored])/length(obs)))
}

#-->> END general utility functions

