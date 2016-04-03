######################################################################
#
# fingers.R
#
# copyright (c) 2001-2016, Karl W Broman
# May and July, 2001
# Last modified Apr, 2016
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/fingers package
# Contains: fingers, parse.hclust
#           fingers2, parse.hclust.all
#
######################################################################

# reproduce Bill Black's FINGERS program to infer
# groups of full siblings based on RAPD data.

# dist = distance matrix
# cutoff = cutoff value; if NULL, one must include "truefam" instead
# method = hierachical clustering method (see hclust)
# truefam = true family structures (include if cutoff is NULL,
#           in which case we find the cutoff that maximizes the
#           adjusted Rand index
# make.plot = plot the dendrogram
# just.plot = just make the plot; don't return the inferred families

fingers <-
function(dist,cutoff=NULL,method=c("average","complete",
         "mcquitty","single","ward"),truefam=NULL,
         make.plot=FALSE,just.plot=FALSE)
{
  method <- match.arg(method)

  n.ind <- nrow(dist)
  ind.names <- rownames(dist)

  # calculate distance matrix; put in form that hclust needs
  d <- as.dist(dist)

  # run hclust
  ho <- hclust(d,method=method)

  # plot the hclust tree?
  if(make.plot || just.plot) {
    plot(ho)
    if(!is.null(cutoff))
      abline(h=cutoff,lwd=2,lty=2)
  }

  if(!just.plot) {
    if(is.null(cutoff)) { # find best possible cutoff
      if(is.null(truefam))
        stop("Include either cutoff or truefam as an argument.")
      co <- parse.hclust.all(ho,truefam)
    }
    else co <- cutoff

    if(make.plot) abline(h=co,lwd=2,lty=2)

    # figure out the inferred families:
    #     hclust merges that are below the cutoff
    fam <- parse.hclust(ho,co)

    # add individual names
    if(!is.null(rownames(dist)))
      fam <- lapply(fam,function(x,y) { names(x) <- y[x]; x }, ind.names)

    fam <- lapply(fam,sort)
    attr(fam,"cutoff") <- co
    return(fam)
  }
  else invisible()

}


# Parse the output from hclust
#     This was a real pain!
#     (Thus I'm too tired to add comments to the code.)
#
parse.hclust <-
function(hclust.out, cutoff)
{
  n.ind <- nrow(hclust.out$merge)+1
  # hclust.out$merge is a (n.ind-1) x 2 matrix
  #     indicating the successive "merges" in the
  #     hierarchical clustering
  #         negative numbers indicate individual numbers
  #         positive numbers indicate cluster numbers
  #
  # height gives, in essence, the value at which the
  #     clusters merge.
  #
  # We look at the clusters which result when we chop
  #     the tree at "cutoff"
  #
  if(cutoff > max(hclust.out$height)) {
    # all one cluster
    fam <- list(1:n.ind)
  }
  else if(cutoff < min(hclust.out$height)) {
    # all in separate clusters
    fam <- as.list(1:n.ind)
  }
  else if(cutoff < hclust.out$height[2]) {
    # two clusters
    a <- 1:n.ind
    b <- hclust.out$merge[1,]
    fam <- as.list(a[b])
    fam[[length(fam)+1]] <- -b
  }
  else {
    m <- hclust.out$merge[hclust.out$height < cutoff,]

    fam <- vector("list",0)
    n.fam <- 0

    for(i in 1:nrow(m)) {
      if(m[i,1] < 0 && m[i,2] < 0) {   # totally new family
        fam[[n.fam+1]] <- -c(m[i,1],m[i,2])
        n.fam <- n.fam + 1
        names(fam)[[n.fam]] <- paste(":",":",sep=as.character(i))
      }
      else {
        if(m[i,1] < 0) {  # m[i,1] is to be added
          x <- paste(":",":",sep=as.character(m[i,2]))
          p <- grep(x,names(fam))
          fam[[p]] <- c(fam[[p]], -m[i,1])
          names(fam)[p] <- paste(names(fam[p]),
                                 paste(":",":",sep=as.character(i)),
                                 sep="")
        }
        else {
          if(m[i,2] < 0) { # m[i,2] is to be added
            x <- paste(":",":",sep=as.character(m[i,1]))
            p <- grep(x,names(fam))
            fam[[p]] <- c(fam[[p]], -m[i,2])
            names(fam)[p] <- paste(names(fam[p]),
                                   paste(":",":",sep=as.character(i)),
                                   sep="")
          }
          else { # merge groups
            x1 <- paste(":",":",sep=as.character(m[i,1]))
            x2 <- paste(":",":",sep=as.character(m[i,2]))
            p1 <- grep(x1,names(fam))
            p2 <- grep(x2,names(fam))

            pmin <- min(c(p1,p2))
            pmax <- max(c(p1,p2))

            fam[[pmin]] <- c(fam[[pmin]],fam[[pmax]])
            names(fam)[pmin] <- paste(names(fam[pmin]),
                                      names(fam[pmax]),
                                      sep="");
            names(fam)[pmin] <- paste(names(fam[pmin]),
                                      paste(":",":",sep=as.character(i)),
                                      sep="")
            fam <- fam[-pmax]
            n.fam <- n.fam - 1
          }
        }
      }
    }
  }

  names(fam) <- NULL

  # order families according to length
  fam <- fam[rev(order(sapply(fam,length)))]

  # find lone individuals and tack on to end
  loners <- (1:n.ind)[is.na(sapply(1:n.ind,function(x,y) match(x,y),unlist(fam)))]
  fam <- c(fam,as.list(loners))

  fam

}





# This is just like parse.hclust, but
#   runs through *all* possible cutoffs, and calculates
#   the adjusted Rand index for each, using a given set
#   of "true" families, and returns the "ideal" cutoff.
#
parse.hclust.all <-
function(hclust.out, truefam)
{
  n.ind <- nrow(hclust.out$merge)+1

  m <- hclust.out$merge
  fam <- vector("list",0)
  n.fam <- 0

  h <- hclust.out$height
  h <- c(apply(cbind(h[-1],h[-length(h)]),1,mean),max(h)+0.5)
  index <- 1:length(h)

  for(i in 1:nrow(m)) {
    if(m[i,1] < 0 && m[i,2] < 0) {   # totally new family
      fam[[n.fam+1]] <- -c(m[i,1],m[i,2])
      n.fam <- n.fam + 1
      names(fam)[[n.fam]] <- paste(":",":",sep=as.character(i))
    }
    else {
      if(m[i,1] < 0) {  # m[i,1] is to be added
        x <- paste(":",":",sep=as.character(m[i,2]))
        p <- grep(x,names(fam))
        fam[[p]] <- c(fam[[p]], -m[i,1])
        names(fam)[p] <- paste(names(fam[p]),
                               paste(":",":",sep=as.character(i)),
                               sep="")
      }
      else {
        if(m[i,2] < 0) { # m[i,2] is to be added
          x <- paste(":",":",sep=as.character(m[i,1]))
          p <- grep(x,names(fam))
          fam[[p]] <- c(fam[[p]], -m[i,2])
          names(fam)[p] <- paste(names(fam[p]),
                                 paste(":",":",sep=as.character(i)),
                                 sep="")
        }
        else { # merge groups
          x1 <- paste(":",":",sep=as.character(m[i,1]))
          x2 <- paste(":",":",sep=as.character(m[i,2]))
          p1 <- grep(x1,names(fam))
          p2 <- grep(x2,names(fam))

          pmin <- min(c(p1,p2))
          pmax <- max(c(p1,p2))

          fam[[pmin]] <- c(fam[[pmin]],fam[[pmax]])
          names(fam)[pmin] <- paste(names(fam[pmin]),
                                    names(fam[pmax]),
                                    sep="")
          names(fam)[pmin] <- paste(names(fam[pmin]),
                                    paste(":",":",sep=as.character(i)),
                                    sep="")
          fam <- fam[-pmax]
          n.fam <- n.fam - 1
        }
      }
    } # end else

    loners <- (1:n.ind)[is.na(sapply(1:n.ind,
                                     function(x,y) match(x,y), unlist(fam)))]

    index[i] <- cluster.stat(c(fam,as.list(loners)),truefam,method="adj")

  } # end of loop

  h[index==max(index)][1]
}


# end of fingers.R
