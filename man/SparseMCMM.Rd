\name{SparseMCMM}
\alias{SparseMCMM}
\title{A main function in SparseMCMM framework}
\description{

This function provides estimates of DE, ME, TE, and component-wise MEs. Additionally, it calculates the statistical significances of OME and CME using a permutation procedure grounded on Models (1) and (2).

}


\usage{
SparseMCMM(Treatment, otu.com, outcome, n.split=10,
                    dirichlet.penalty=seq(0,1,0.1),
                    lm.penalty1=seq(0,1,0.1), lm.penalty2=seq(0,2,0.2),
                    low.bound1=NULL, up.bound1=NULL, low.bound2=NULL, up.bound2=NULL,
                    num.per=NULL)
}


\arguments{
  \item{Treatment}{
A numeric vector of the binary treatment (takes the value 1 if it is assigned to the treatment group and takes the value 0 if assigned to the control group) with length = sample size (n).
}
  \item{otu.com}{
  A n*p numeric matrix containing compositional microbiome data. Each row represents a subject, and each column represents a taxon (given the rank, for example, the genus rank) or an OTU. The row sum equals 1.
}

  \item{outcome}{
A numeric vector of the continuous outcome with length = sample size (n).
}

  \item{n.split}{ An integer value. The number of repetitions regarding the split strategy. See \code{Details} for
a more comprehensive discussion on the split strategy. Default=10.}

  \item{dirichlet.penalty}{
A numeric vector that includes potential tuning parameters employed during the estimation process in relation to Model 2. Default=seq(0,1,0.1).
}

  \item{lm.penalty1}{
A numeric vector that includes potential tuning parameters employed during the estimation process in relation to Model 1. Default=seq(0,1,0.1).
}

  \item{lm.penalty2}{
A numeric vector that includes potential tuning parameters employed during the estimation process in relation to Model 1. Default=seq(0,2,0.2).
}

  \item{low.bound1}{
A numeric vector representing the lower bounds of the controls employed during the estimation process in relation to Model 1. Default=NULL, there is no lower bound.
}

  \item{up.bound1}{
A numeric vector representing the upper bounds of the controls employed during the estimation process in relation to Model 1. Default=NULL, there is no upper bound.
}

  \item{low.bound2}{
A numeric vector representing the lower bounds of the controls employed during the estimation process in relation to Model 2. Default=NULL, there is no lower bound.
}

  \item{up.bound2}{
A numeric vector representing the upper bounds of the controls employed during the estimation process in relation to Model 2. Default=NULL, there is no upper bound.
}

  \item{num.per}{
An integer value, the number of permuations. statistical significances of tests TME and CME are calculated based on these permutations. Default=NULL, No calculation for hypothesis test.
}

}

\value{

A list which contains three elements:
  \describe{
  \item{Esitmated Causal Effects}{Estimates of direct effect (DE), mediation effect (ME) and total effect (TE). See \code{Details} for a more comprehensive discussion.}
  \item{Estimated component-wise Mediation Effects}{Estimates of component-wise MEs for all mediators. See \code{Details} for a more comprehensive discussion.}

  \item{Test}{P-values for tests OME and CME if num.per is not NULL, otherwise, no this item.}

  }}

\details{

Within the SparseMCMM framework, regularization techniques are employed to carry out variable selection, aiding in the identification of signature causal microbes. To account for the biases introduced by the regularization techniques employed, we further implement splitting strategy (Rinaldo et al; 2019), which can handle arbitrary penalties and provide asymptotically validated inference. Specifically, we randomly divide the dataset into two equal halves: the first half is utilized for variable selection, while the second half is dedicated to parameter estimation. The estimates of DE, ME, TE, and component-wise MEs were then calculated. We repeated this data splitting procedure multiple times (n.split times) to ensure robustness and accuracy in our estimations and inference.

If n.split > 1, the 'Esitmated Causal Effects' component is a 2 x 3 matrix. The first row contains the average estimates of DE, ME, and TE  based on n.split times of repetitions. Correspondingly, the second row provides the standard errors associated with these averages. If n.split=1,  the 'Estimated Causal Effects' component is a vector of length 3. This vector contains the estimates for DE, ME, and TE, which are calculated based on a single, random data split.

If n.split > 1, the 'Estimated component-wise Mediation Effects' component is a 3 x p matrix.  The first row contains the average estimates of compoment-wise MEs for all p taxa based on n.split times of repetitions. Correspondingly, the second and third rows provide the lower and upper 95\% confidence interval associated with these averages. If n.split=1,  the 'EEstimated component-wise Mediation Effects' component is a vector of length p. This vector contains the estimates for compoment-wise MEs which are calculated based on a single, random data split.
}


\references{
Wang C, Hu J, Blaser MJ, and Li H (2020). Estimating and testing the microbial causal mediation effect with high-dimensional and compositional microbiome data. Bioinformatics. 36(2):347-355.

Wang C, Ahn J, Tarpey T, Stella SY, Hayes RB and Li H (2023). A microbial causal mediation analytic tool for health disparity and applications in body mass index.

Rinaldo A, Wasserman L, G’Sell M (2019). Bootstrapping and sample splitting for high-dimensional, assumption-lean inference.

}
\author{
Chan Wang and Huilin Li.
}

\examples{

# ########### generation data
# ### Sample size and number of mediators
# sample.num=100
# otu.num=10
#
# ###Treatment
# Treatment=rep(c(0,1),each=sample.num/2)
#
# ### parameters
# beta0=c(0.6, -0.3, 0.8, -1.4, -1.2, -1.4, -1.3, -1.0, -0.2, 0.6)
# betaT=rep(0,otu.num)
# betaT[c(1,3)]=c(0.4,0.2)
#
# alpha0=0
# alphaT=1
# alphaZ=alphaC=rep(0,otu.num)
# alphaZ[c(1,3)]=c(0.7,-0.7)
# alphaC[c(1,3)]=c(0.15,-0.15)
#
#
# ############Microbiome data
# library(dirmult)
# X=cbind(rep(1,sample.num),Treatment)
# b=cbind(beta0,betaT)
# gamma.simu=exp(X %*% t(b))
# otu.com=t(apply(gamma.simu,1,rdirichlet,n=1))
#
# ##################Outcome  data
# X=cbind(rep(1,sample.num),Treatment,log(otu.com),log(otu.com)*Treatment)
# b=c(alpha0,alphaT,alphaZ,alphaC)
# outcome=c(b%*%t(X)+rnorm(sample.num,mean = 0, sd =1))
#
#
# ##### SparseMCMM function
# SparseMCMM(Treatment,otu.com,outcome,n.split=1,
#           num.per=10)
}
\keyword{SparseMCMM}







