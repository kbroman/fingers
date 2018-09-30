\name{cutoff.llr}
\alias{cutoff.llr}

\title{Calculate cutoff for clustering with RAPD markers}

\description{
  Calculate a cutoff (for the LLR distance measure) for hierarchical
  cluster analysis to infer groups of related individuals with RAPD data.
}

\usage{
cutoff.llr(f,method=c("qu","meansib","qs","lr"),value=0.2)
}

\arguments{
 \item{f}{A vector of band allele frequencies for a set of RAPD markers.}
 \item{method}{The method to use to form the cutoff: a quantile of the
   distribution of distances among unrelated (\code{qu}), the mean
   distance between siblings (\code{meansib}), a quantile of the
   distribution of distances among siblings (\code{qs}), or the
   likelihood ratio for unrelateds vs. siblings (\code{lr}).}
 \item{value}{For \code{method="qu"} or \code{method="qs"}, this should
   specify the quantile; for \code{method="lr"}, this should specify the
   likelihood ratio.}
}

\value{
  The cutoff (a single value).
}

\author{Karl W Broman \email{broman@wisc.edu}}

\seealso{ \code{\link{cutoff}},
  \code{\link{llrdist}},
  \code{\link{freq}},
  \code{\link{pull.markers}},
  \code{\link{fingers}} }

\examples{
data(aedes)
f <- freq(aedes)
co1 <- cutoff.llr(f,method="meansib")
co2 <- cutoff.llr(f,method="qu",value=0.2)
co3 <- cutoff.llr(f,method="qs",value=0.9)
co4 <- cutoff.llr(f,method="lr",value=4.0)
}

\keyword{cluster}
