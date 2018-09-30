\name{aedes}
\alias{aedes}

\title{Data on Aedes aegypti}
\description{
  This is RAPD data for 40 loci typed on a set of 10 full-sibling
  families, with 15 individuals in each family.
}
\usage{data(aedes)}
\format{
  The data is a matrix of 150 rows (the individuals) by 40 columns (the
  RAPD loci).  Each entry is a RAPD phenotype, indicating the presence
  (1) or absence (0) of a band.
}

\references{
  BL Apostol, WC Black IV, BR Miller, P Reiter, BJ Beaty (1993)
  Estimation of the number of full sibling families at an oviposition
  site using RAPD-PCR markers: applications to the mosquito \emph{Aedes
  aegypti}.  Theor Appl Genet 86:991-1000.
}

\source{FINGERS software, WC Black IV, Colorado State University \cr
  \url{http://csu-cvmbs.colostate.edu/academics/mip/aidl/Pages/black-william.aspx}
}

\author{Karl W Broman \email{broman@wisc.edu}}


\examples{
data(aedes)
}

\seealso{\code{\link{shiff1}}, \code{\link{simrapd}} }

\keyword{datasets}
