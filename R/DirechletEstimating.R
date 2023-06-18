#######################Parameter estimation for the Dirichlet regression model

######The log likelihood function for the Dirichlet regression

likelihood2.org=function(betaA,Treatment,otu.com){

  num.otu=ncol(otu.com); sample.num=nrow(otu.com);

  beta0=betaA[1:num.otu];
  betaT=betaA[-(1:num.otu)];

  X=cbind(rep(1,sample.num),Treatment) #n*(1+p)
  b=cbind(beta0,betaT) #num.otu*(1+p)
  g=exp(X %*% t(b))	# n * num.otu

  ### log likelihood

  sum(lgamma(rowSums(g)) + rowSums((g-1) * log(otu.com) - lgamma(g)))
}


######The log likelihood function for the Dirichlet regression with lasso penalty for variable selection

likelihood2=function(beta,Treatment,otu.com,lambda.dirichlet){

  num.otu=ncol(otu.com);
  beta.pen=beta[-(1:num.otu)]

  ### -log likelihood+ penalty

  -likelihood2.org(beta,Treatment,otu.com)+lambda.dirichlet*(sum(abs(beta.pen)))
}


######The gradient function for the log likelihood with lasso penalty

gradient2 = function(betaA,Treatment,otu.com,lambda.dirichlet) {

  num.otu=ncol(otu.com); sample.num=nrow(otu.com);

  beta0=betaA[1:num.otu];
  betaT=betaA[-(1:num.otu)];

  X=cbind(rep(1,sample.num),Treatment) #n*2
  b=cbind(beta0,betaT) #num.otu*2
  g=exp(X %*% t(b))	# n * num.otu

  S=t((digamma(rowSums(g))  - digamma(g) + log(otu.com)) * g) %*% X  ## num.otu*2

  der=rep(0,length(betaT))
  if(sum(betaT!=0)>0) der[betaT!=0]=lambda.dirichlet*betaT[betaT!=0]/abs(betaT[betaT!=0])

  c(-S)+c(rep(0,num.otu),der)

}

#### function to optimize - a list of objective and gradient
toOpt2 = function(betaA,Treatment,otu.com,lambda.dirichlet){
  list(objective=likelihood2(betaA,Treatment,otu.com,lambda.dirichlet),
       gradient=gradient2(betaA,Treatment,otu.com,lambda.dirichlet))
}
