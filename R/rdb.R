#' @title Robust Differential Abundance Test
#' 
#' @description \code{rdb} is used to conduct robust differential abundance analysis for compositional data
#' 
#' @details This function returns an indicator vector, where each entry correponds to each component of compoositional data. 
#' 
#' @importFrom stats median
#' @importFrom stats var
#' 
#' @param P A numerical matrix for compositional data, each row represents a sample (the sum should be 1) and each column represents a component. 
#' @param Z A binary vector, 1 means treated group and 0 means control group.
#' @param X A numerical matrix for observed covariates, each row represents a sample and each column represents a covariates.
#' @param alpha A numerical value, indicating the significance level of FWER or FDR. 
#' @param fdr A boolean value, indicating if FDR is the measure of type I error. TRUE means FDR is used and FALSE means FWER is used.
#' 
#' @return  an indicator vector, where each entry correponds to each column of P. TRUE means it is differential component.
#' 
#' @export
#' 
#' @examples
#' m=50
#' d=100 
#' P=matrix(runif(m*d),nrow=m,ncol=d)
#' Z=rep(0,m)
#' Z[1:(m/2)]=1
#' rdb(P,Z)
#' 
#' @author Shulei Wang
rdb <- function(P, Z, X=NULL, alpha=0.1, fdr=FALSE)
{
  if (nrow(P)!=length(Z))
    stop("Please make sure the number of rows in P equal to the length of Z")
  if (!all(Z %in% c(0, 1)))
    stop("Please make sure Z only contain 0 and 1")
  for (i in 1:nrow(P))
  {
    P[i,]=P[i,]/sum(P[i,])
  }
  if (!is.null(X))
  {
    if (is.vector(X)) X<-matrix(X,ncol=1)
    if (nrow(P)!=nrow(X))
      stop("Please make sure the number of rows in P equal to number of rows in X")
  }
  
  filter=apply(P, 2, sum)>0
  P.filter=P[,filter]
  
  treat=Z==1
  mtreat=sum(treat)
  control=Z==0
  mcontrol=sum(control)
  d=ncol(P.filter)
  M=sqrt(2*log(d)/d)
  D=sqrt(2*log(d)-2*log(alpha))
  Dpm=D+0.2*M
  
  meanvar=matrix(0,nrow = 5,ncol = d)
  if (is.null(X)) {
    Ptreat=P.filter[treat,]
    Pcontrol=P.filter[control,]
    meanvar[1,]=apply(Ptreat, 2, mean)
    meanvar[2,]=apply(Pcontrol, 2, mean)
    meanvar[3,]=apply(Ptreat, 2, var)/mtreat
    meanvar[5,]=apply(Pcontrol, 2, var)/mcontrol
  } else {
    for (j in 1:d) {
      tRe=ATE (P.filter[,j], Z, X)
      meanvar[1:2,j]=tRe$est[1:2]
      meanvar[3,j]=tRe$vcov[1,1]
      meanvar[4,j]=tRe$vcov[1,2]
      meanvar[5,j]=tRe$vcov[2,2]
    }
  }
  
  if (fdr == FALSE)
  {
    Vt <- rdb.core(meanvar, M, D, Dpm)
  } else {
    Dlist<-seq(0,D,length.out=100)
    pFDR<-rep(0,length(Dlist))
    for (j in 1:length(Dlist))
    {
      DD = Dlist[j]
      tVt <- rdb.core(meanvar, M, DD, DD+0.2*M)
      R=max(1,sum(!tVt))
      R0=2*d*pnorm(DD,lower.tail =FALSE)
      pFDR[j]=R0/R
    }
    if (pFDR[100]<alpha)
    {
      FThre <- Dlist[min(which(pFDR<=alpha))]
    } else {
      FThre <- Dlist[100]
    }
    Vt <- rdb.core(meanvar, M, FThre, FThre+0.2*M)
  }
  
  ComUt<-rep(FALSE,ncol(P))
  ComUt[filter]=!Vt
  
  return(ComUt)
}

rdb.core <- function(meanvar, M, D, Dpm)
{
  d=ncol(meanvar)
  Vt=rep(TRUE,d)
  while (sum(Vt)>0) {
    Rtreat=sum(meanvar[1,Vt])
    Rcontrol=sum(meanvar[2,Vt])
    if (Rtreat<=0 && Rcontrol>0) {
      tstat=meanvar[2,]/sqrt(meanvar[5,])
    } else if (Rtreat>0 && Rcontrol<=0) {
      tstat=meanvar[1,]/sqrt(meanvar[3,])
    } else if (Rtreat>0 && Rcontrol>0) {
      tstat=(meanvar[1,]/Rtreat-meanvar[2,]/Rcontrol)/sqrt(meanvar[3,]/Rtreat/Rtreat-2*meanvar[4,]/Rtreat/Rcontrol+meanvar[5,]/Rcontrol/Rcontrol)
    } else {
      break
    }
    Mtstat=median(tstat[Vt])
    if (Mtstat>M) {
      Wt=Vt&(tstat<(-D))
    } else if (Mtstat<(-M)) {
      Wt=Vt&(tstat>D)
    } else {
      Wt=Vt&(abs(tstat)>Dpm)
    }
    if (sum(Wt)==0){
      break
    }
    Vt[Wt]=FALSE
  }
  return(Vt)
}
