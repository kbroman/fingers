\name{calc.dist}
\alias{calc.dist}
\title{Calculate simple distance matrix}

\description{
  Calculate the simple distance matrix, by the proportion of mismatches,
  for a RAPD data set.
}

\usage{
calc.dist(dat)
}

\arguments{
 \item{dat}{A matrix of size (n.ind x n.mar) containing RAPD
   phenotypes, with 1 indicating the presence of a band and 0 indicating
   absence.}
}

\details{
  For each pair of individuals, we calculate the proportion of RAPD
  markers (among those where both individuals have complete data) at
  which one individual shows a band and the other doesn't.
}

\value{
  A symmetric matrix of dimension (n.ind x n.ind), containing the
  distances between individuals.
}

\author{Karl W Broman \email{broman@wisc.edu}}

\seealso{ \code{\link{llrdist}}, \code{\link{fingers}} }

\references{
  BL Apostol, WC Black IV, BR Miller, P Reiter, BJ Beaty (1993)
  Estimation of the number of full sibling families at an oviposition
  site using RAPD-PCR markers: applications to the mosquito \emph{Aedes
  aegypti}.  Theor Appl Genet 86:991-1000.
}

\examples{
data(aedes)
d <- calc.dist(aedes)
}

\keyword{cluster}
