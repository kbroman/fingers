######################################################################
#
# freq.R
#
# copyright (c) 2001-2, Karl W Broman
# Last modified Nov, 2002
# First written May, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/fingers package
# Contains: freq, pull.markers
#
######################################################################

# Calculate est'd (observed) frequencies of the band allele
#     dat is a matrix of 0's and 1's where 1 = band
#     (rows = worms; columns = loci)
freq <-
function(dat)
  1 - sqrt(1 - apply(dat,2,mean,na.rm=TRUE))

# Pull out only markers having band allele freq between 0.1 and 0.6
#
#  f = allele frequencies is one of the arguments
#      (if it's not given, it is calculated from the data
#
pull.markers <-
function(dat,lo=0.1,hi=0.6,f=freq(dat)) 
  dat[,f > lo & f < hi]


# end of freq.R
