\name{shiff2}
\alias{shiff2}

\title{Schistosome data}
\description{
  This is RAPD data for 10 loci typed on a set of 135 individuals.
  Markers with estimated band allele frequencies outside of the range
  0.1-0.6 have been removed.
}
\usage{data(shiff2)}
\format{
  The data is a matrix of 135 rows (the individuals) by 10 columns (the
  RAPD loci).  Each entry is a RAPD phenotype, indicating the presence
  (1) or absence (0) of a band.
}

\source{Clive Shiff, Molecular Microbiology and Immunology,
  Bloomberg School of Public Health, The Johns Hopkins University
}

\author{Karl W Broman \email{broman@wisc.edu}}

\examples{
data(shiff2)
}

\seealso{\code{\link{shiff1}}, \code{\link{shiff3}},
  \code{\link{aedes}}, \code{\link{simrapd}} }

\keyword{datasets}
