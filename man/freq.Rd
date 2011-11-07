\name{freq}
\alias{freq}

\title{Estimate RAPD allele frequencies}

\description{
  Estimate the frequency of the band allele for a set of RAPD markers.
}

\usage{
freq(dat)
}

\arguments{
 \item{dat}{A matrix of size (n.ind x n.mar) containing RAPD
   phenotypes, with 1 indicating the presence of a band and 0 indicating
   absence.} 
}

\details{
  The RAPDs are assumed to be in Hardy-Weinberg equilibrium, and so the
  frequency of the band allele is estimated as \eqn{\hat{p} =
    1-\sqrt{1-\bar{x}}}{phat = 1 - sqrt(1 - xbar)} where
  \eqn{\bar{x}}{xbar} is the proportion of individuals showing a band.
}  

\value{
  A vector of length n.mar, containing the estimated frequencies of the
  band allele for each RAPD marker.
}

\author{Karl W Broman \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{pull.markers}} }

\references{
  BL Apostol, WC Black IV, BR Miller, P Reiter, BJ Beaty (1993)
  Estimation of the number of full sibling families at an oviposition
  site using RAPD-PCR markers: applications to the mosquito \emph{Aedes
  aegypti}.  Theor Appl Genet 86:991-1000.
}

\examples{
data(aedes)
f <- freq(aedes)
}
                    
\keyword{univar}
