#require(nloptr)

####### parameter estimation for the log contrast model regression
############ Parameter estimation for the linear log-contrast model
#################################### the sum of squared residuals (SSR) measuring the discrepancy
#####################between the observed outcome and predicted outcome

SSR2=function(alpha,Treatment,Z,outcome,Index_0)
{
  sample.n=nrow(Z)

  if(!is.null(Index_0)) { alphaAll=rep(0, ncol(Z)*2+2); alphaAll[-Index_0]=alpha} else alphaAll=alpha

  X=cbind(rep(1,sample.n),Treatment,Z,Z*Treatment)

  sum((outcome-c(X%*%alphaAll))^2)
}


######The gradient function for SSR with penalties

alphagradient2 = function(alpha,Treatment,Z,outcome,Index_0){

  n=ncol(Z);sample.n=nrow(Z)

  if(!is.null(Index_0)) { alphaAll=rep(0, ncol(Z)*2+2); alphaAll[-Index_0]=alpha} else alphaAll=alpha

  X=cbind(rep(1,sample.n),Treatment,Z,Z*Treatment)

  if(!is.null(Index_0)) return(c(-2*t(X)%*%(outcome-c(X%*%alphaAll)))[-Index_0]) else return( c(-2*t(X)%*%(outcome-c(X%*%alphaAll))))
}



### equality constraints
equal2 <- function(alpha,Treatment,Z,outcome,Index_0) {

  otu.num=ncol(Z);

if(!is.null(Index_0)) { alphaAll=rep(0, ncol(Z)*2+2); alphaAll[-Index_0]=alpha} else alphaAll=alpha

  z1=t(c(rep(0,2),rep(1,otu.num),rep(0,otu.num)))%*%alphaAll
  z2=t(c(rep(0,otu.num+2),rep(1,otu.num)))%*%alphaAll

  return(c(z1,z2))
}


#equality constraint function.  The jacobian is 1 for all variables
eqCon2 = function(alpha,Treatment,Z,outcome,Index_0) {

  otu.num=ncol(Z);
  jacobian=rbind(c(rep(0,2),rep(1,otu.num),rep(0,otu.num)),
                 c(rep(0,2+otu.num),rep(1,otu.num)))

  if(!is.null(Index_0)) jacobian=jacobian[,-Index_0]

  list(constraints=equal2(alpha,Treatment,Z,outcome,Index_0),
       jacobian=jacobian)
}


#function to optimize - a list of objective and gradient
alphatoOpt2 = function(alpha,Treatment,Z,outcome,Index_0){
  list(objective=SSR2(alpha,Treatment,Z,outcome,Index_0),
       gradient=alphagradient2(alpha,Treatment,Z,outcome,Index_0))
}



alpha.estimates2=function(Treatment,otu.com,outcome,Index_0) {
  ### using a small value to replace 0
  pseudo=min(otu.com[otu.com>0])/2
  otu.com=t(apply(otu.com,1,function(x) {if(min(x)==0) return((x+pseudo)/sum(x+pseudo)) else return(x)}))

  sample.num=nrow(otu.com);p=ncol(otu.com)

alpha.initial=c(mean(outcome[Treatment==0]),rep(0,(1+2*p-length(Index_0))))

qq=nloptr(alpha.initial,
          alphatoOpt2,
          eval_g_eq=eqCon2,
           opts =list( "algorithm" = "NLOPT_LD_SLSQP",
                       "xtol_rel" = 1.0e-4,
                       "maxeval"= 5000),
           Treatment=Treatment,Z=log(otu.com),outcome=outcome,
          Index_0=Index_0)


return(qq$solution)

}


#require(Compositional)
#require(nloptr)

#######################Parameter estimation for the Dirichlet regression model
######The log likelihood function for the Dirichlet regression

betalikelihood2.org=function(beta,Treatment,otu.com,Index_0=NULL){

  sample.num=nrow(otu.com); num.otu=ncol(otu.com)

  beta0=beta[1:num.otu];
  betaT=beta[-(1:num.otu)];

  if(!is.null(Index_0)) { betaTAll=rep(0, length(beta0)); betaTAll[-Index_0]=betaT} else betaTAll=betaT

  ##### betaC is  num.otu*q

  X=cbind(rep(1,sample.num),Treatment) #n*(1+q+p)

  b=cbind(beta0,betaTAll) #num.otu*(1+q+p)

  g=exp(X %*% t(b))	# n * num.otu

  ### log likelihood

  -sum(lgamma(rowSums(g)) + rowSums((g-1) * log(otu.com) - lgamma(g)))
}


######The gradient function for the log likelihood

betagradient2 = function(beta,Treatment,otu.com, Index_0=NULL) {

  sample.num=nrow(otu.com); num.otu=ncol(otu.com)
  beta0=beta[1:num.otu];
  betaT=beta[-(1:num.otu)];

  X=cbind(rep(1,sample.num),Treatment) #n*(1+p)

  if(!is.null(Index_0)) { betaTAll=rep(0, length(beta0)); betaTAll[-Index_0]=betaT} else betaTAll=betaT

  b=cbind(beta0,betaTAll) #num.otu*(1+p)

  g=exp(X %*% t(b))	# n * num.otu

  S=t((digamma(rowSums(g))  - digamma(g) + log(otu.com)) * g) %*% X  ## 1+p


  if(!is.null(Index_0)) return(c(-S)[-c(Index_0+num.otu)]) else return(c(-S))

}


beta.estimates2=function(Treatment,otu.com, Index_0=NULL)
{
  ### using a small value to replace 0
  pseudo=min(otu.com[otu.com>0])/2
  otu.com=t(apply(otu.com,1,function(x) {if(min(x)==0) return((x+pseudo)/sum(x+pseudo)) else return(x)}))

  sample.num=nrow(otu.com);p=ncol(otu.com)

  beta0.intial=diri.est(otu.com[Treatment==0,],type = "mle")$param
  beta.initial=c(beta0.intial,rep(0,p-length(Index_0)))

  qq <- optim(beta.initial,
              betalikelihood2.org,
              betagradient2,
              method  = "BFGS",
              control = list(maxit=5000),
             otu.com=otu.com,Treatment=Treatment,
             Index_0=Index_0)

  return(qq$par)

}



