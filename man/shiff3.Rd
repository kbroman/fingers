\name{shiff3}
\alias{shiff3}

\title{Schistosome data}
\description{ 
  This is RAPD data for 10 loci typed on a set of 125 individuals.
  Markers with estimated band allele frequencies outside of the range
  0.1-0.6 have been removed.  Individuals with one or more missing
  values have been removed.
} 
\usage{data(shiff3)}
\format{
  The data is a matrix of 125 rows (the individuals) by 10 columns (the
  RAPD loci).  Each entry is a RAPD phenotype, indicating the presence
  (1) or absence (0) of a band.
}

\source{Clive Shiff, Molecular Microbiology and Immunology,
  Bloomberg School of Public Health, The Johns Hopkins University
}

\author{Karl W Broman \email{kbroman@biostat.wisc.edu}}

\examples{
data(shiff3)
}

\seealso{\code{\link{shiff1}}, \code{\link{shiff2}},
  \code{\link{aedes}}, \code{\link{simrapd}} }

\keyword{datasets}
