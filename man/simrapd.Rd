\name{simrapd}
\alias{simrapd}

\title{Simulate RAPD data}

\description{
  Simulates RAPD data for a set of sibling families.
}

\usage{
simrapd(n.sib = rep(15,10), p = c(rep(0.125,8),rep(0.175,5),rep(0.225,5),
           rep(0.275,8),rep(0.325,3),rep(0.375,4),
           rep(0.475,4),rep(0.575,3)))
}

\arguments{
 \item{n.sib}{A vector giving the number of siblings per family (length
   is the number of families).}
 \item{p}{A vector of frequencies of the band allele at each marker
   (length is the number of markers).}
}

\details{
  The RAPDs are assumed to be in Hardy-Weinberg equilibrium.
}

\value{
  A matrix of dimension (n.ind x n.mar), giving the RAPD phenotypes for
  each individual at each marker, with 1 indicating a band and 0
  indicating no band.
}

\author{Karl W Broman \email{broman@wisc.edu}}

\seealso{ \code{\link{simulfams}} }

\references{
  BL Apostol, WC Black IV, BR Miller, P Reiter, BJ Beaty (1993)
  Estimation of the number of full sibling families at an oviposition
  site using RAPD-PCR markers: applications to the mosquito \emph{Aedes
  aegypti}.  Theor Appl Genet 86:991-1000.
}

\examples{
data <- simrapd(rep(20,5), p=runif(40, 0.1, 0.6))
}

\keyword{datagen}
