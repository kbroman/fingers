######################################################################
#
# sims.R
#
# copyright (c) 2001-2, Karl W Broman
# Last modified Nov, 2002
# First written May, 2001
# Licensed under the GNU General Public License version 2 (June, 1991)
# 
# Part of the R/fingers package
# Contains: simrapd, inherit
#           simulfams, inherit2, samppop
#
######################################################################

######################################################################
#
# simrapd: simulate RAPD data for a set of sibling families
#
# input:
#
#   n.sib = vector of length n.fam, containing nos of sibs per family
#
#   p     = vector of frequency of the band allele (length = n.mar)
#
######################################################################

simrapd <-
function(n.sib = rep(15,10),
         p=c(rep(0.125,8),rep(0.175,5),rep(0.225,5),
           rep(0.275,8),rep(0.325,3),rep(0.375,4),
           rep(0.475,4),rep(0.575,3)))
{
  n.fam <- length(n.sib)
  n.mar <- length(p)

  mom.g <- as.data.frame(matrix(rbinom(n.fam*n.mar,2,rep(p,n.fam)),ncol=n.fam))
  dad.g <- as.data.frame(matrix(rbinom(n.fam*n.mar,2,rep(p,n.fam)),ncol=n.fam))

  kid.g <- matrix(unlist(lapply(rbind(mom.g,n.sib),inherit)) +
                  unlist(lapply(rbind(dad.g,n.sib),inherit)),nrow=n.mar)
  kid.g[kid.g == 2] <- 1
  kid.g <- t(kid.g)
  dimnames(kid.g) <- list(paste(rep(1:n.fam,n.sib),unlist(lapply(n.sib,function(x) 1:x)),
                                sep="-"), paste("M",1:n.mar,sep=""))

  kid.g

}

         
######################################################################
#
# inherit: a subroutine for simrapd
#
######################################################################

inherit <-
function(par.g)
{
  n <- par.g[length(par.g)]
  par.g <- par.g[-length(par.g)]
  x <- matrix(rep(par.g,n),ncol=n)
  apply(x,2,function(a) { a[a==1] <- sample(0:1,sum(a==1),repl=TRUE); a[a==2] <- 1; a })
}

         
######################################################################
#
# simulfams: basically the same as sibrapd
#            (written by Laura Plantinga)
#
######################################################################

simulfams <-
function(n.sib=sample(5:20,size=sample(5:20,size=1),repl=TRUE),
         p=runif(sample(5:15,size=1),min=0.1,max=0.6))
{
  # this function simulates a set of full-sibling families 
  # assuming no half-siblings and parents not in sample
  # n.fam is the number of families
  # n.loci is the number of loci evaluable with RAPD analysis
  n.fam <- length(n.sib)
  n.loci <- length(p)

  # gen.prob is a vector of Hardy-Weinberg genotype frequencies
  gen.prob <- cbind((1-p)^2, 2*p*(1-p), p^2)

  # gen will contain genotypes for all individuals in the dataset
  gen <- NULL



  for (i in 1:n.fam) {
			
    # gen.father = genotype of father (0=aa,1=Aa,2=AA)
    # gen.mother = genotype of mother (0=aa,1=Aa,2=AA)
    # lengths of gen.father and gen.mother are n.loci
    # gen.sib is a matrix where rows = siblings and columns = loci
    gen.father <- apply(gen.prob,1,samppop)
    gen.mother <- apply(gen.prob,1,samppop)
    gen.parents <- cbind(gen.father,gen.mother)
    gen.sib <- apply(gen.parents,1,inherit2,n.sib[i])
    if (i==1) gen <- gen.sib else gen<- rbind(gen,gen.sib)
  }
        
  # convert genotypes (0,1,2) to phenotypes (0=absence of band, 1=presence of band)
  gen[gen==2] <- 1

  # name dimensions of gen, columns = # of locus, rows = "family number"."individual number"
  dimnames(gen) <- list(paste(rep(1:n.fam,n.sib),1:sum(n.sib),sep="."),1:n.loci)

  return(gen)
}


# this function simulates Mendelian inheritance of genes
inherit2 <-
function(g,sib)
{

  # g[1] = genotype of father (0=aa,1=Aa,2=AA)
  # g[2] = genotype of mother (0=aa,1=Aa,2=AA)
  # pf = prob. of father passing on A allele
  # pm = prob. of mother passing on A allele
  # g1 = allele from father (0=a, 1=A)
  # g2 = allele from mother (0=a, 1=A)
  # sib = number of siblings in family
  
  pf <- g[1]/2
  pm <- g[2]/2
  g1 <- sample(0:1,prob=c(1-pf,pf),size=sib, repl=TRUE)
  g2 <- sample(0:1,prob=c(1-pm,pm),size=sib, repl=TRUE)
  g1+g2
}

samppop <-
function(p)
{
  # this function assigns a genotype based on values of p for population
  sample(0:2,1,prob=p)		
}

# end of sims.R
