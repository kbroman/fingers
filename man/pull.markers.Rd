\name{pull.markers}
\alias{pull.markers}

\title{Extract markers with allele frequencies in specified range}

\description{
  Extract markers from a RAPD data set that have allele frequencies
  within a specified range.
}

\usage{
pull.markers(dat,lo=0.1,hi=0.6,f=freq(dat))
}

\arguments{
 \item{dat}{A matrix of size (n.ind x n.mar) containing RAPD
   phenotypes, with 1 indicating the presence of a band and 0 indicating
   absence.}
 \item{lo}{Lower bound for band allele frequency.}
 \item{hi}{Upper bound for band allele frequency.}
 \item{f}{Vector of band allele frequencies (included in order to avoid
   recalculating it, if possible).}
}

\value{
  A matrix, like the argument \code{dat}, but containing only those
  markers with band allele frequency between \code{lo} and \code{hi}.
}

\author{Karl W Broman \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{freq}} }

\examples{
data(shiff1)
f <- freq(shiff1)
subset <- pull.markers(shiff1, 0.1, 0.6, f)
}
                    
\keyword{manip}
