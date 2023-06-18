
####################### Causal direct effect, mediation effect and total effect

CausalE=function(otu.com,alpha.estimation,beta.estimation)
{
  otu.num=ncol(otu.com)

  cov.1=1; cov.2=c(1,1)
  beta.matrix=matrix(beta.estimation,nrow=otu.num)

  alphaT=alpha.estimation[2]; alphaZC=alpha.estimation[(length(alpha.estimation)-2*otu.num+1):length(alpha.estimation)]

  alphaZ=alphaZC[1:otu.num]; alphaC=alphaZC[-(1:otu.num)]

  gammaT.0=as.vector(exp(cov.1%*%t(beta.matrix[,-ncol(beta.matrix)]))) ; gammaT.1=as.vector(exp(cov.2%*%t(beta.matrix)))

  ET.0=digamma(gammaT.0)-digamma(sum(gammaT.0))
  ET.1=digamma(gammaT.1)-digamma(sum(gammaT.1))


  nie.each=(alphaZ+alphaC)*(ET.1-ET.0)

  DE=alphaT+sum(alphaC*ET.0)
  ME=sum(nie.each)
  TE=DE+ME
  ME.each=nie.each

  AllEffect=c(DE,ME,TE)
  names(AllEffect)=c("DE","ME","TE")

  if(is.null(colnames(otu.com))) names(ME.each)=1:otu.num else names(ME.each)=colnames(otu.com)

  CausalEf=list(AllEffect,ME.each)

  names(CausalEf)=c("Esitmated Causal Effects","Estimated component-wise Mediation Effects")

  return(CausalEf)
}



