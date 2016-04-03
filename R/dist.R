######################################################################
#
# dist.R
#
# copyright (c) 2001-2016, Karl W Broman
# Last modified Nov, 2002
# First written Apr, 2016
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/fingers package
# Contains: calc.dist, llrdist
#
######################################################################


# calculate D = 1-M matrix
#    d(i,j) = proportion of non-matches between
#             worms i and j
calc.dist <-
function(dat)
{
  # make sure dat is a matrix,
  # 'cause the stuff below screws up if it's a data.frame
  dat <- as.matrix(dat)

  # data and (1-data) with NA's replaced by 0's
  X1 <- dat; X1[is.na(X1)] <- 0
  X2 <- 1-dat; X2[is.na(X2)] <- 0

  # 1's where there is data; 0's where there is not
  Y <- dat
  Y[!is.na(Y)] <- 1
  Y[is.na(Y)] <- 0

  dist <- (X1 %*% t(X2) + X2 %*% t(X1)) / (Y %*% t(Y))
  dimnames(dist) <- list(rownames(dat),rownames(dat))
  dist
}



######################################################################
#
# llrdist
#
# Function to calculate the distance matrix using the log likelihood
# ratio in place of simply proportion mismatches for the Fingers
# program.
#
# calculations are done in C
#
######################################################################

llrdist <-
function(dat, p=freq(dat))
{
#  if(!is.loaded("llrdist_wrap")) {
#    lib.file <- file.path("./", paste("llrdist", .Platform$dynlib.ext,sep=""))
#    dyn.load(lib.file)
#    cat(paste(" -Loaded", lib.file), "\n")
#  }

  nind <- nrow(dat)
  nmar <- ncol(dat)

  z <- .C("llrdist_wrap",
          as.integer(nind),
          as.integer(nmar),
          as.integer(dat),
          as.double(p),
          dist=as.double(rep(0,nind^2)),
          PACKAGE="fingers")

  z <- matrix(z$dist,ncol=nind)
  dimnames(z) <- list(rownames(dat),rownames(dat))
  z

}


# end of dist.R
