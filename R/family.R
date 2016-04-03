######################################################################
#
# family.R
#
# copyright (c) 2001-2016, Karl W Broman
# Last modified Nov, 2002
# First written Apr, 2016
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/fingers package
# Contains: true.fams, cluster.stat, comp.fams, dist.image,
#           cluster.stat.sub
#
######################################################################

# Identify the true family structure using the rownames in
# a RAPD data set
true.fams <-
function(dat)
{
  a <- rownames(dat)
  b <- 1:length(a)
  true <- as.numeric(sapply(strsplit(a,"-"), function(y) y[1]))
  x <- split(a,true)
  y <- split(b,true)
  for(i in 1:length(x)) names(y[[i]]) <- x[[i]]
  y
}

# Statistic for measuring the similarity between true
#   and estimated clustering
cluster.stat <-
function(fam1,fam2,method=c("all","rand","adj","fm","kb"))
{
  method <- match.arg(method)
  result <- NULL
  n <- length(unlist(fam2))
  nc2 <- choose(n,2)

  # rand or kb
  if(method=="rand" || method=="kb" || method=="all") {
    a <- cluster.stat.sub(fam1)
    b <- cluster.stat.sub(fam2)

    X <- outer(a,a,function(a,b) a==b)
    Y <- outer(b,b,function(a,b) a==b)

    if(method=="rand" || method == "all") {
      result <- c(result,
                  sum(((X & Y) | (!X & !Y)) & row(X) > col(X))/
                  choose(length(unlist(fam1)),2))
    }
    if(method=="kb" || method=="all") {
      n.same.tf <- sum(sapply(fam2,function(a) choose.rev(length(a))))
      n.diff.tf <- nc2 - n.same.tf
      n.same.both <- sum((X & Y) & row(X) > col(X))
      n.diff.both <- sum((!X & !Y) & row(X) > col(X))

      result <- c(result,
                  0.5*((n.same.both/n.same.tf) + (n.diff.both/n.diff.tf)))
    }
  }

  # adj-Rand or Fowlkes and Mallows
  if(method=="adj" || method=="fm" || method=="all") {
    f1 <- cbind(unlist(fam1), rep(1:length(fam1), sapply(fam1, length)))
    f2 <- cbind(unlist(fam2), rep(1:length(fam2), sapply(fam2, length)))
    tab <- table(f2[order(f2[, 1]), 2], f1[order(f1[, 1]), 2])

    n.row <- apply(tab,1,sum)
    n.col <- apply(tab,2,sum)

    if(method=="adj" || method=="all") {
      n1 <- sum(choose.rev(n.row))
      n2 <- sum(choose.rev(n.col))
      a <- n1*n2/nc2
      N <- sum(choose.rev(tab))-a
      D <- (n1+n2)/2-a
      result <- c(result, N/D)
    }
    if(method=="fm" || method=="all") {
      P <- sum(n.row^2)-n
      Q <- sum(n.col^2)-n
      R <- sum(tab^2)-n
      if(P*Q == 0) result <- c(result,0)
      else result <- c(result,R/sqrt(P*Q))
    }
  }
  if(method == "all") {
    result <- result[c(1,3,4,2)]
    names(result) <- c("rand","adj","fm","kb")
  }
  result
}

# n choose 2, returning 0 when n < 2
choose.rev <-
function(n)
{
  x <- n
  x[n<2] <- 0
  x[n>=2] <- choose(n[n>=2],2)
  x
}


cluster.stat.sub <-
function(x)
  rep(1:length(x),sapply(x,length))[match(1:length(unlist(x)),unlist(x))]


# Compare two different sets of inferred families for the same data
#     The first component of the output is a two-way table comparing
#         the two family groupings
#     The second component lists the individuals according to the first
#         family grouping and their corresponding family in the second
#         family grouping.
#
comp.fams <-
function(fam1,fam2)
{
  f1 <- cbind(unlist(fam1),rep(1:length(fam1),sapply(fam1,length)))
  f2 <- cbind(unlist(fam2),rep(1:length(fam2),sapply(fam2,length)))

  comp <- lapply(fam1,function(x,y) y[match(x,y[,1]),],f2)

  x <- table(f2[order(f2[,1]),2],f1[order(f1[,1]),2])

  if(all(apply(x,1,function(x) sum(x!=0))==1))
    cat("    All families are the same.\n")
  else cat("    There are some differences.\n")
  list(x,comp)

}

# use "image" to plot the distance matrix
#
# if "fams" is given, reorder so that families
#     are together
#
dist.image <-
function(dist,fams=NULL,col=topo.colors(1+ncol(dist)),...)
{
  n <- nrow(dist)

  if(!is.null(fams)) {
    dist <- dist[unlist(fams),unlist(fams)]
    fam.len <- cumsum(c(0,sapply(fams,length)))+0.5
  }

  image(1:n,1:n,dist,xlab="Worms",ylab="Worms",main="Distance matrix",
        col=col,...)

  a <- par("usr")
  abline(v=a[1:2])
  abline(h=a[3:4])

  if(!is.null(fams)) {
    abline(h=fam.len)
    abline(v=fam.len)
    for(i in 2:length(fam.len)) {
      text(mean(fam.len[i+(-1:0)]),a[4]+diff(a[3:4])*0.02,as.character(i-1),xpd=TRUE)
      text(a[2]+diff(a[1:2])*0.02,mean(fam.len[i+(-1:0)]),as.character(i-1),xpd=TRUE)
    }
  }
}


#
# determine quality of results of the cluster analysis
#
#check.fams <-
#function(fam)
#{
#  true <- lapply(fam,function(x) as.numeric(sapply(strsplit(names(x),"-"),
#                                            function(y) y[1])))
#
#  vote <- sapply(true,function(x) {
#    y <- table(x)
#    as.numeric(names(y)[y==max(y)])[1] })
#
#  o <- order(vote)
#  true <- true[o]
#  vote <- vote[o]
#
#  n.fam <- length(unique(unlist(true)))
#
#  x <- matrix(ncol=length(vote),nrow=n.fam)
#  dimnames(x) <- list(as.character(1:n.fam),as.character(vote))
#  for(j in 1:length(vote))
#    x[,j] <- table(factor(true[[j]],levels=1:n.fam))
#
#  y <- sum(sapply(true,function(x) {
#    z <- table(x)
#    zz <- names(z)[z==max(z)][1]
#    sum(z[names(z) != zz]) }))
#
#  list(table=x,n.err=y)
#
#}

# Add names to a set of families, using the row names from data
#add.names <-
#function(fam,dat)
#{
#  lapply(fam,function(x,y) { names(x) <- y[x]; x }, rownames(dat))
#}

# end of family.R
