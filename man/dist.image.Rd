\name{dist.image}
\alias{dist.image}

\title{Plot distance matrix}

\description{
  Plot the distance matrix for a RAPD data set, with (optionally) lines
  drawn separating clusters of individuals.
}

\usage{
dist.image(dist,fams=NULL,col=topo.colors(1+ncol(dist)),\dots)
}

\arguments{
 \item{dist}{A matrix of size (n.ind x n.ind), containing the distances
   between pairs of individuals.}
 \item{fams}{A list of clusters; each component in the list is one
   inferred family, containing the indices of individuals placed in that
   family.}
 \item{col}{Colors to use in the plot; see \code{\link[graphics]{image}}.}
 \item{\dots}{Other arguments to pass to \code{\link[graphics]{image}}.}
}

\value{
  The function calls \code{\link[graphics]{image}} in order to create an
  image of the distance matrix.
}

\author{Karl W Broman \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{calc.dist}},
  \code{\link{true.fams}} }

\examples{
data(aedes)
f <- freq(aedes)
co <- cutoff(f)
d <- calc.dist(aedes)
fam <- fingers(d,co,make.plot=TRUE)
dist.image(d,fam)
}

\keyword{cluster}
