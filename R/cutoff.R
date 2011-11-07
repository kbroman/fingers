######################################################################
#
# cutoff.R
#
# copyright (c) 2001-6, Karl W Broman
# last modified Oct, 2006
# first written May, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/fingers package
# Contains: cutoff, cutoff.llr
#           Msib, Mur, avesd.Dur, avesd.Dsib 
#
######################################################################

# Probability of matching phenotype for siblings
Msib <-
function(p)
  p^4+4*p^3*(1-p)+2*p*(1-p)^3+4.5*p^2*(1-p)^2+(1-p)^4

# Prob of matching phe for unrelateds
Mur <-
function(p)
 (1-p)^4+(1-(1-p)^2)^2

avesd.Dur <-
function(p)
{
  x <- 1-Mur(p)
  ave <- mean(x)
  sd <- sqrt(sum(x*(1-x)))/length(x)
  c(ave=ave,sd=sd)
}

avesd.Dsib <-
function(p)
{
  x <- 1-Msib(p)
  ave <- mean(x)
  sd <- sqrt(sum(x*(1-x)))/length(x)
  c(ave=ave,sd=sd)
}


cutoff <-
function(f,method=c("qu","meansib","qs","lr"),value=0.2)
{
  method <- match.arg(method)

  # mean distance among siblings
  if(method=="meansib") return(1-mean(Msib(f)))

  # quantile from unrelateds
  if(method=="qu") 
    return(.C("cutoff",
              as.integer(length(f)), 
              as.double(f),
              as.double(value),
              cutoff=as.double(0),
              PACKAGE="fingers")$cutoff)
    
  if(method=="qs") 
    return(.C("cutoff_s",
              as.integer(length(f)), 
              as.double(f),
              as.double(value),
              cutoff=as.double(0),
              PACKAGE="fingers")$cutoff)

  if(method=="lr")
    return(.C("cutoff_lr",
              as.integer(length(f)), 
              as.double(f),
              as.double(value),
              cutoff=as.double(0),
              PACKAGE="fingers")$cutoff)
}


cutoff.llr <-
function(f,method=c("qu","meansib","qs","lr"),value=0.2)
{
  method <- match.arg(method)

  # mean distance among siblings
  if(method=="meansib")
    return(.C("cutoff_llr_meansib",
              as.integer(length(f)),
              as.double(f),
              cutoff=as.double(0),
              PACKAGE="fingers")$cutoff)

  # quantile from unrelateds
  if(method=="qu") 
    return(.C("cutoff_llr",
              as.integer(length(f)), 
              as.double(f),
              as.double(value),
              cutoff=as.double(0),
              PACKAGE="fingers")$cutoff)
    
  if(method=="qs") 
    return(.C("cutoff_llr_s",
              as.integer(length(f)), 
              as.double(f),
              as.double(value),
              cutoff=as.double(0),
              PACKAGE="fingers")$cutoff)

  if(method=="lr")
    return(.C("cutoff_llr_lr",
              as.integer(length(f)), 
              as.double(f),
              as.double(value),
              cutoff=as.double(0),
              PACKAGE="fingers")$cutoff)
}

# end of cutoff.R
