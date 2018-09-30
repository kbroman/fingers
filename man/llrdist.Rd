\name{llrdist}
\alias{llrdist}

\title{Calculate distance matrix based on log likelihood ratio}

\description{
  Calculate a distance matrix, based on the log likelihood ratio
  comparing the hypotheses of full sibling versus unrelated,
  for a RAPD data set.
}

\usage{
llrdist(dat,p=freq(dat))
}

\arguments{
 \item{dat}{A matrix of size (n.ind x n.mar) containing RAPD
   phenotypes, with 1 indicating the presence of a band and 0 indicating
   absence.}
 \item{p}{A vector of band allele frequencies.}
}

\details{
  For each pair of individuals, at each locus, we calculate the log
  likelihood ratio (LLR) comparing the hypotheses \emph{unrelated} with
  \emph{siblings}, with the data being B (both have band), N (neither
  have band) or D (one has band, the other doesn't).  These LLRs are
  \emph{averaged} across individuals.  \bold{Note}: at each
  locus, we re-center the LLRs so that the minimum of the LLRs among
  B/N/D is 0; this makes the resulting distances \eqn{\ge}{>} 0.

  Calculations are performed in a C program.
}

\value{
  A symmetric matrix of dimension (n.ind x n.ind), containing the
  distances between individuals.
}

\author{Karl W Broman \email{broman@wisc.edu}}

\seealso{ \code{\link{calc.dist}}, \code{\link{fingers}} }

\examples{
data(aedes)
f <- freq(aedes)
dis <- llrdist(aedes,f)
}

\keyword{cluster}
