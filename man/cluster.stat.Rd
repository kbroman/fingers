\name{cluster.stat}
\alias{cluster.stat}

\title{Calculate measure of quality of inferred clusters}

\description{
  Calculate a score indicating how well two sets of clusters 
  conform.
}

\usage{
cluster.stat(fam1,fam2,method=c("all","rand","adj","fm","kb"))
}

\arguments{
 \item{fam1}{A list of clusters; each component in the list is one
   family, containing the indices of the individuals in that family.}
 \item{fam2}{A list, just like \code{fam1}.}
 \item{method}{A character string indicating whether to calculate the
   Rand index, the adjusted Rand index, the Fowlkes and Mallows B index,
   or Karl Broman's index. If \code{method=all}, a vector with all four
   indices is returned.}
}

\details{
  In the Rand index (Rand 1971), one considers all pairs of individuals,
  and assigns a 1 to a pair if the individuals are either in the same
  cluster in both \code{fam1} and \code{fam2} or are not in the same
  cluster in both \code{fam1} and \code{fam2}, and assigns a 0 to the pair
  otherwise, and then takes the sum of these, divided by the number of
  pairs of individuals.

  Karl Broman's index (which we don't recommend, but we implement
  here in order to allow comparisons to be made) is just like the Rand
  index, but \code{fam2} is assumed to be the \emph{true} partition, and
  the set of all pairs in the same group (by \code{fam2}) and the set of
  all pairs in different groups (by \code{fam2}), are given equal weight. 

  Let \eqn{n_{ij}}{n[ij]} be the number of individuals in group i by
  partition 1 and group j by partition 2.  Let \eqn{n_{i.} = \sum_{j}
    n_{ij}}{n[i.] = sum j n[ij]} and define \eqn{n_{.j}}{n[.j]}
  similarly.  

  In the adjusted-Rand index (Hubert and Arabie 1985), \dots

  In the Fowlkes and Mallows B index (Fowlkes and Mallows 1983), \dots

}

\value{
  The value of a score for comparing two sets of clusters.
}

\author{Karl W Broman \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{fingers}},
  \code{\link{true.fams}} }

\references{
  WM Rand (1971) Objective criteria for the evaluation of clustering
  methods.  Journal of the American Statistical Association 66:846-850.

  L Hubert and P Arabie (1985) Comparing partitions.  Journal of
  Classification.  2:193-218.

  EB Fowlkes and CL Mallows (1983) A method for comparing two
  hierarchical clusterings.  Journal of the American Statistical
  Association 78:553-584.

  BS Everitt, S Landau and M Leese (2001) Cluster analysis, 4th
  edition.  Arnold, London, pp. 181-3.
}

\examples{
data(aedes)
f <- freq(aedes)
co <- cutoff(f)
d <- calc.dist(aedes)
fam <- fingers(d,co,make.plot=TRUE)
tf <- true.fams(aedes)
cluster.stat(fam,tf)
cluster.stat(fam,tf,method="fm")
}
                    
\keyword{cluster}
