\name{comp.fams}
\alias{comp.fams}

\title{Compare two sets of clusters}

\description{
  Give diagnostic information indicating how well two sets of clusters
  conform. 
}

\usage{
comp.fams(fam1,fam2)
}

\arguments{
 \item{fam1}{A list of clusters; each component in the list is one
   family, containing the indices of the individuals in that family.}
 \item{fam2}{A list, just like \code{fam1}.}
}

\value{
  A list with two components.  The first component is a contingency
  table whose (i,j)th element is the number of individuals in cluster i
  in \code{fam1} and cluster j in \code{fam2}.  The second component is
  a list indicating, for each cluster from \code{fam1}, the cluster
  assignment in \code{fam2}.
}

\author{Karl W Broman \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{cluster.stat}},
  \code{\link{fingers}},
  \code{\link{true.fams}} }

\examples{
data(aedes)
f <- freq(aedes)
co <- cutoff(f)
d <- calc.dist(aedes)
fam <- fingers(d,co,make.plot=TRUE)
tf <- true.fams(aedes)
comp.fams(fam,tf)
}
                    
\keyword{cluster}
