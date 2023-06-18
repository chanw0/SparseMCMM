############ Parameter estimation for the linear log-contrast model

#################################### the sum of squared residuals (SSR) measuring the discrepancy
#####################between the observed outcome and predicted outcome

SSR=function(alpha,Treatment,Z,outcome)
{
  sample.n=nrow(Z)

  X=cbind(rep(1,sample.n),Treatment,Z,Z*Treatment)

  sum((outcome-c(X%*%alpha))^2)
}


########### SSR with penalties
SSR.penalty=function(alpha,Treatment,Z,outcome,lambda)
{
  n=ncol(Z);sample.n=nrow(Z)

  alphaZC=alpha[-(1:2)]
  alphaZ=alphaZC[1:n];alphaC=alphaZC[-(1:n)]

  SSR(alpha,Treatment,Z,outcome)+
    lambda[1]*sum(sqrt(alphaZ^2+alphaC^2))+lambda[2]*(sum(abs(alphaC)))
}


######The gradient function for SSR with penalties

gradient = function(alpha,Treatment,Z,outcome,lambda){

  n=ncol(Z);sample.n=nrow(Z)

  alphaZC=alpha[-(1:2)]
  alphaZ=alphaZC[1:n];alphaC=alphaZC[-(1:n)]

  X=cbind(rep(1,sample.n),Treatment,Z,Z*Treatment)

  rr=c(rep(0,2),1/sqrt(alphaZ^2+alphaC^2),1/sqrt(alphaZ^2+alphaC^2))
  rr[which(alphaZ==0)+2]=0
  rr[which(alphaC==0)+2+n]=0

  der=rep(0,n)
  if(sum(alphaC!=0)>0) der[alphaC!=0]=lambda[2]*alphaC[alphaC!=0]/abs(alphaC[alphaC!=0])

  c(-2*t(X)%*%(outcome-c(X%*%alpha))+lambda[1]*alpha*rr)+c(rep(0,2+n),der)
}


#function to optimize - a list of objective and gradient
toOpt = function(alpha,Treatment,Z,outcome,lambda){
  list(objective=SSR.penalty(alpha,Treatment,Z,outcome,lambda),
       gradient=gradient(alpha,Treatment,Z,outcome,lambda))
}


### equality constraints
equal <- function(alpha,Treatment,Z,outcome,lambda) {

  otu.num=ncol(Z);

  z1=t(c(rep(0,2),rep(1,otu.num),rep(0,otu.num)))%*%alpha
  z2=t(c(rep(0,otu.num+2),rep(1,otu.num)))%*%alpha

  return(c(z1,z2))
}

#equality constraint function.  The jacobian is 1 for all variables
eqCon = function(alpha,Treatment,Z,outcome,lambda) {
  otu.num=ncol(Z);

  list(constraints=equal(alpha,Treatment,Z,outcome,lambda),
       jacobian=rbind(c(rep(0,2),rep(1,otu.num),rep(0,otu.num)),
                      c(rep(0,2+otu.num),rep(1,otu.num))))
}

