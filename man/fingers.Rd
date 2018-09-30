\name{fingers}
\alias{fingers}

\title{Infer clusters of related individuals}

\description{
  Perform hierarchical clustering to infer groups of related individuals
  with RAPD data.
}

\usage{
fingers(dist,cutoff=NULL,method=c("average","complete",
         "mcquitty","single","ward"),truefam=NULL,
         make.plot=FALSE,just.plot=FALSE)
}

\arguments{
 \item{dist}{A matrix of size (n.ind x n.ind) containing the distances
   between individuals.}
 \item{cutoff}{A value to use to cut off the dendogram formed by
   hierarchical clustering in order to define a set of
   clusters. (Optional, but if NULL, the argument \code{truefam}
   must be included.)}
 \item{method}{A hierarchical clustering method.  See
   \code{\link[stats]{hclust}}. \emph{Note:} We haven't allowed
   \code{centroid} or \code{median}, because these weren't working for
   us.}
 \item{truefam}{The true family structure; used only if \code{cutoff} is
   NULL, in which case all possible cutoffs are tried, and that
   giving the maximum adjusted Rand index is used.}
 \item{make.plot}{If TRUE, make a plot of the dendogram formed by
   hierarchical clustering.}
 \item{just.plot}{If TRUE, just make the plot; don't return the inferred
   families.  (In this case, the \code{cutoff} argument is not needed.)}
}

\details{
  We use the function \code{\link[stats]{hclust}}
  to do the cluster analysis.
}

\value{
  A list of clusters; each component in the list is one inferred family,
  containing the indices of individuals placed in that family.  The
  cutoff used is included as an attribute.  Use
  \code{attr(result,"cutoff")} to obtain this value.
}

\author{Karl W Broman \email{broman@wisc.edu}}

\seealso{ \code{\link{cutoff}},
  \code{\link{cutoff.llr}},
  \code{\link{calc.dist}},
  \code{\link{llrdist}},
  \code{\link{cluster.stat}},
  \code{\link{true.fams}},
  \code{\link{freq}},
  \code{\link{pull.markers}} }

\references{
  BL Apostol, WC Black IV, BR Miller, P Reiter, BJ Beaty (1993)
  Estimation of the number of full sibling families at an oviposition
  site using RAPD-PCR markers: applications to the mosquito \emph{Aedes
  aegypti}.  Theor Appl Genet 86:991-1000.
}

\examples{
data(aedes)
f <- freq(aedes)
co <- cutoff(f)
d <- calc.dist(aedes)
fam <- fingers(d,co,make.plot=TRUE)
tf <- true.fams(aedes)
cluster.stat(fam,tf)
}

\keyword{cluster}
