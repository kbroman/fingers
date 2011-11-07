\name{simulfams}
\alias{simulfams}

\title{Simulate RAPD data}

\description{
  Simulates RAPD data for a set of sibling families.
}

\usage{
simulfams(n.sib=sample(5:20,size=sample(5:20,size=1),repl=TRUE),
          p=runif(sample(5:15,size=1),min=0.1,max=0.6))
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

\author{Laura Plantinga and Karl Broman \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{simrapd}} }

\references{
  BL Apostol, WC Black IV, BR Miller, P Reiter, BJ Beaty (1993)
  Estimation of the number of full sibling families at an oviposition
  site using RAPD-PCR markers: applications to the mosquito \emph{Aedes
  aegypti}.  Theor Appl Genet 86:991-1000.
}

\examples{
data <- simulfams(rep(20,5), p=runif(40, 0.1, 0.6))
}
                    
\keyword{datagen}
