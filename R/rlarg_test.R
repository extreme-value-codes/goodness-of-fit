# rlarg.test: goodness-of-fit tests for the r largest order statistics model
# (parametric bootstrap). Requires ismev::rlarg.fit.
#
# Returns an r x 8 matrix. Columns 1-4 are the observed statistics
# (probability-plot correlation, Kolmogorov-Smirnov D, Cramer-von Mises W,
# Anderson-Darling A) for each of the r order statistics; columns 5-8 are
# the corresponding bootstrap p-values.
rlarg.test=function(xdat, r = dim(xdat)[2], ydat = NULL, mul = NULL, sigl = NULL,
                    shl = NULL, mulink = identity, siglink = identity, shlink = identity,
                    muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE,
                    method = "Nelder-Mead", maxit = 10000, nsim){
  fit1=rlarg.fit(xdat, r = r, ydat = ydat, mul = mul, sigl = sigl,
                 shl = shl, mulink = mulink, siglink = siglink, shlink = shlink,
                 muinit = muinit, siginit = siginit, shinit = shinit, show = show,
                 method = method, maxit = maxit)
  u1=sort(exp(-(1+fit1$vals[,3]*(xdat[,1]-fit1$vals[,1])/fit1$vals[,2])^(-1/fit1$vals[,3])))
  n=length(u1)
  plotpos=qnorm((1:n-0.5)/n)
  lg=cor(qnorm(u1),plotpos)
  d1<-max(c((1:n)/n-u1,u1-(0:(n-1))/n))
  w1<-sum((u1-((1:n)-0.5)/n)^2)+1/(12*n)
  a1<--sum((2*(1:n)-1)*log(u1)+(2*n+1-2*(1:n))*log(1-u1))/n-n
  ddf=matrix(nrow=r,ncol=4)
  ddf[1,]=c(lg,d1,w1,a1)
  for(k in 1:(r-1)){
    u1=sort(exp(
      -(1+fit1$vals[,3]*(xdat[,k+1]-fit1$vals[,1])/fit1$vals[,2])^(-1/fit1$vals[,3])
      +(1+fit1$vals[,3]*(xdat[,k]-fit1$vals[,1])/fit1$vals[,2])^(-1/fit1$vals[,3])
    ))
    n=length(u1)
    plotpos=qnorm((1:n-0.375)/(n+0.25))
    lg=cor(qnorm(u1),plotpos)
    d1<-max(c((1:n)/n-u1,u1-(0:(n-1))/n))
    w1<-sum((u1-((1:n)-0.5)/n)^2)+1/(12*n)
    a1<--sum((2*(1:n)-1)*log(u1)+(2*n+1-2*(1:n))*log(1-u1))/n-n
    ddf[k+1,]=c(lg,d1,w1,a1)
  }
  # simulation
  n=nrow(xdat)
  xsim=matrix(nrow=n,ncol=r)
  count=matrix(0,ncol=4,nrow=r)
  for(isim in 1:nsim){
    n=nrow(xdat)
    xsim[,1]=fit1$vals[,1]+fit1$vals[,2]*((-log(runif(n)))^(-fit1$vals[,3])-rep(1,n))/fit1$vals[,3]
    for(k in 1:(r-1)){
      n=nrow(xdat)
      xsim[,k+1]=fit1$vals[,1]+fit1$vals[,2]*((((-log(runif(n)))+(1+fit1$vals[,3]*(xsim[,k]-fit1$vals[,1])/fit1$vals[,2])
                                                ^(-1/fit1$vals[,3]))^(-fit1$vals[,3]))-1)/fit1$vals[,3]
    }
    for(i in 1:n){for(j in 1:r){if(is.na(xdat[i,j]))xsim[i,j]=NA}}
    fit2=rlarg.fit(xsim, r = r, ydat = ydat, mul = mul, sigl = sigl,
                   shl = shl, mulink = mulink, siglink = siglink, shlink = shlink,
                   muinit = muinit, siginit = siginit, shinit = shinit, show = show,
                   method = method, maxit = maxit)
    u1=sort(exp(-(1+fit2$vals[,3]*(xsim[,1]-fit2$vals[,1])/fit2$vals[,2])^(-1/fit2$vals[,3])))
    n=length(u1)
    plotpos=qnorm((1:n-0.375)/(n+0.25))
    lg2=cor(qnorm(u1),plotpos)
    d2<-max(c((1:n)/n-u1,u1-(0:(n-1))/n))
    w2<-sum((u1-((1:n)-0.5)/n)^2)+1/(12*n)
    a2<--sum((2*(1:n)-1)*log(u1)+(2*n+1-2*(1:n))*log(1-u1))/n-n
    ddf2=matrix(nrow=r,ncol=4)
    ddf2[1,]=c(lg2,d2,w2,a2)
    for(k in 1:(r-1)){
      u1=sort(exp(
        -(1+fit2$vals[,3]*(xsim[,k+1]-fit2$vals[,1])/fit2$vals[,2])^(-1/fit2$vals[,3])
        +(1+fit2$vals[,3]*(xsim[,k]-fit2$vals[,1])/fit2$vals[,2])^(-1/fit2$vals[,3])
      ))
      n=length(u1)
      plotpos=qnorm((1:n-0.5)/n)
      lg2=cor(qnorm(u1),plotpos)
      d2<-max(c((1:n)/n-u1,u1-(0:(n-1))/n))
      w2<-sum((u1-((1:n)-0.5)/n)^2)+1/(12*n)
      a2<--sum((2*(1:n)-1)*log(u1)+(2*n+1-2*(1:n))*log(1-u1))/n-n
      ddf2[k+1,]=c(lg2,d2,w2,a2)
    }
    count[,1]=count[,1]+(ddf2[,1]<ddf[,1])
    count[,2]=count[,2]+(ddf2[,2]>ddf[,2])
    count[,3]=count[,3]+(ddf2[,3]>ddf[,3])
    count[,4]=count[,4]+(ddf2[,4]>ddf[,4])
  }
  ddf=cbind(ddf,count/nsim)
  return(ddf)
}
