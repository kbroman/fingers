\name{true.fams}

\alias{true.fams}

\title{Identify the true clusters}

\description{
  Use the row names of a RAPD data set to identify the true sets of
  families. 
}

\usage{
true.fams(dat)
}

\arguments{
 \item{dat}{A matrix of size (n.ind x n.mar) containing RAPD
   phenotypes, with 1 indicating the presence of a band and 0 indicating
   absence.  The row names (identifying individuals) are assumed to be
   of the form "family-individual"}
}

\value{
  A list of clusters; each component in the list is one inferred family,
  containing the indices of individuals placed in that family.
}

\author{Karl W Broman \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{aedes}}, \code{\link{simrapd}},
  \code{\link{fingers}},
  \code{\link{cluster.stat}} }

\references{
  BL Apostol, WC Black IV, BR Miller, P Reiter, BJ Beaty (1993)
  Estimation of the number of full sibling families at an oviposition
  site using RAPD-PCR markers: applications to the mosquito \emph{Aedes
  aegypti}.  Theor Appl Genet 86:991-1000.
}

\examples{
data(aedes)
tf <- true.fams(aedes)
}
                    
\keyword{cluster}
