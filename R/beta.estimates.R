#require(Compositional)
#require(nloptr)

beta.estimates=function(Treatment,otu.com,penalty.lambda=seq(0,1,0.1),
                        low.bound=NULL,up.bound=NULL)
{

  ### using a small value to replace 0
  pseudo=min(otu.com[otu.com>0])/2
  otu.com=t(apply(otu.com,1,function(x) {if(min(x)==0) return((x+pseudo)/sum(x+pseudo)) else return(x)}))

  sample.num=nrow(otu.com);num.otu=ncol(otu.com)


  beta0.intial=diri.est(otu.com[Treatment==0,],type = "mle")$param
  beta.initial=c(beta0.intial,rep(0,num.otu))


  BIC=estimates=NULL

  for(lambda.dirichlet in penalty.lambda)
  {
    qq=nloptr( x0=beta.initial,
               toOpt2,
               lb=low.bound,
               ub=up.bound,
               opts =list("algorithm" = "NLOPT_LD_MMA",
                           "xtol_rel" = 1.0e-4,
                           "maxeval"= 5000),
               otu.com=otu.com,Treatment=Treatment,lambda.dirichlet=lambda.dirichlet)

    if(qq$status %in% c(1,3,4)){

      estimates1=round(qq$solution,3)

      estimates=rbind(estimates,estimates1)

      BIC=c(BIC,log(sample.num)*sum(estimates1!=0)-
              2*likelihood2.org(estimates1,Treatment=Treatment,otu.com=otu.com))
    }}

  if(length(BIC)==0) beta.estimates=beta.initial else beta.estimates=estimates[which.min(BIC),]

  return(beta.estimates)

}
