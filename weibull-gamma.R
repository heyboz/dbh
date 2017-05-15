library(MASS)
library(actuar)
library(fitdistrplus)

csfile = 'cs_BA0'
file1 = paste('E:/BaiduYunSych/GitHub/ClusterSize/Outputs/nearest/all/cs_scanner', csfile, '.csv', sep="")
cs = read.csv(file1,header=TRUE)
#sp = c('MACBRE','QUELIT','LITACU',
#       'DISRAC','ARTSTY','DIOMOR',
#       'SYMLAU','ORMPAC','CINPAR')
n = 67
sp = names(cs)[1:n+1]

weibull.scale = c(1:n)
weibull.shape = c(1:n)
weibull.pvalue = c(1:n)
weibull.P = c(1:n)

gamma.rate = c(1:n)
gamma.shape = c(1:n)
gamma.pvalue = c(1:n)
gamma.P = c(1:n)

lnorm.meanlog = c(1:n)
lnorm.sdlog = c(1:n)
lnorm.pvalue = c(1:n)
lnorm.P = c(1:n)

# llogis.rate = c(1:n)
# llogis.shape = c(1:n)
# llogis.P = c(1:n)
# llogis.pvalue = c(1:n)

powerlaw.scale = c(1:n)
powerlaw.shape = c(1:n)
powerlaw.pvalue = c(1:n)
powerlaw.P = c(1:n)

for(i in c(1:n)){
  #text = paste("read.csv('E:/Outputs/Cluster/nearest/", sp[i], "_0.csv',header=TRUE)", sep="")
  #data = eval(parse(text=data))
  cc = eval(parse(text = paste('cs$',sp[i],sep="")))
  #cc = cs$LITELO
  cc = cc[!is.na(cc)]
  #cc = cc[cc>0]
  cc = cc + 1
  fw = fitdistr(cc, 'weibull')
  wtest = ks.test(cc,"pweibull", scale=fw$estimate[2], shape=fw$estimate[1])
  
#   EPS = sqrt(.Machine$double.eps) # "epsilon" for very small numbers
#   llik.weibull <- function(shape, scale, thres, x){
#     sum(dweibull(x - thres, shape, scale, log=T))}
#   thetahat.weibull <- function(x){
#     if(any(x <= 0)) stop("x values must be positive")
#     toptim <- function(theta) -llik.weibull(theta[1], theta[2], theta[3], x)
#     mu = mean(log(x))
#     sigma2 = var(log(x))
#     shape.guess = 1.2 / sqrt(sigma2)
#     scale.guess = exp(mu + (0.572 / shape.guess))
#     thres.guess = 1
#     res = nlminb(c(shape.guess, scale.guess, thres.guess), toptim, lower=EPS)
#     c(shape=res$par[1], scale=res$par[2], thres=res$par[3])}
#   fw3 = thetahat.weibull(cc)
#   w3test = ks.test(cc+fw3[3], 'pweibull', scale=fw3[2], shape=fw3[1])
  
  fg = fitdistr(cc, 'gamma')
  gtest = ks.test(cc,"pgamma", rate=fg$estimate[2], shape=fg$estimate[1])
  
  fl = fitdistr(cc, 'log-normal')
  ltest = ks.test(cc, 'plnorm', meanlog = fl$estimate[1], sdlog = fl$estimate[2])
  
#   fll = fitdistr(cc, dllogis, start=list(shape = 0.1, rate = 1))
#   lltest = ks.test(cc, 'pllogis', shape=fll$estimate[1], rate=fll$estimate[2])
  
  png(file=paste('E:/Outputs/Cluster/nearest/fitting/PDF_',csfile,'/',
                 i,sp[i],'_fitting','.png',sep=""),
      width = 480, height = 480, units = "px")
  h = hist(cc,main=paste('PDF and Hisogram of ',sp[i],sep=""),freq=F,xlab='Cluster Size',
           breaks=40,border='white',col='grey86',cex.lab=1.5,cex.main=1.5)
  lwd = 2
  yhist = h$density
  cx = seq(min(cc),max(cc),length.out=length(h$breaks)-1)
  fpp = nls(yhist~a*cx^(-b),start=list(a=0.4,b=1))
  fp = summary(fpp)
  scalefit = fp$parameters[1]
  shapefit = (-1)*fp$parameters[2]
  xfit = seq(min(cc),max(cc),length.out=64)
  ypowl = scalefit*xfit^shapefit
  curve(scalefit*x^shapefit,col='turquoise3',lwd=lwd,add=T)
  pp = scalefit*cx^shapefit
  ptest = ks.test(yhist, pp)
  
  
  ywei = dweibull(xfit,scale=fw$estimate[2], shape=fw$estimate[1])
  ygam = dgamma(xfit,rate=fg$estimate[2], shape=fg$estimate[1])
  ylnorm = dlnorm(xfit, meanlog=fl$estimate[1], sdlog=fl$estimate[2])
  
  lines(xfit, ylnorm, col=4,lwd=lwd,add=T)
  lines(xfit, ywei, col=2,lwd=lwd,add=T)
  lines(xfit, ygam, col=3,lwd=lwd,add=T)
  #
  #   yfit3 = dllogis(xfit, shape=fll$estimate[1], rate=fl$estimate[2])
  #   lines(xfit, yfit3, col=6,lwd=2)
  legend ('topright', c("weibull","gamma",'log-normal',"power-law"), 
          col=c(2,3,4,'turquoise3'), cex=1.5, lty=1, lwd=2)
  dev.off()
  
#   data = data.frame(cc=cc,ywei=ywei,ygam=ygam,ylnorm=ylnorm,ypowl=ypowl)
#   p1<-ggplot(data, aes(x=cc)) + 
#     geom_histogram(aes(x=cc, y=..density..), alpha=0.3) +
#     geom_line(aes(x=xfit, y=ywei),size=1.2,col='deepskyblue3') +
#     geom_line(aes(x=xfit, y=ygam),size=1.2,col='brown3') +
#     geom_line(aes(x=xfit, y=ylnorm),size=1.2,col='seagreen4') +
#     geom_line(aes(x=xfit, y=ypowl),size=1.2,col='olivedrab4')
#   p1
  
#   ###fail
#   png(file=paste('E:/Outputs/Cluster/nearest/weibull-gamma/ecdf_wgp/',i,sp[i],'_ew','.png',sep=""))
#   x = seq(0,max(cc),0.1)
#   plot(x,pweibull(x,scale=fw$estimate[2], shape=fw$estimate[1]),type="l",col="red", main="ECDF and Weibull, Gamma& Power-law CDF")
#   lines(x,pgamma(x,rate=fg$estimate[2], shape=fg$estimate[1]),type="l",col="blue", add=T)
#   plot(ecdf(scalefit*x^shapefit),type="l",col="green", add=T)
#   plot(ecdf(cc),add=T)
#   dev.off()
#   
#   png(file=paste('E:/Outputs/Cluster/nearest/weibull-gamma/ecdf_gamma/',i,sp[i],'_eg','.png',sep=""))
#   x = seq(0,max(cc),0.1)
#   plot(x,pgamma(x,rate=fg$estimate[2], shape=fg$estimate[1]),type="l",col="blue", main="ECDF and Gamma CDF")
#   plot(ecdf(cc),add=TRUE)
#   dev.off()
#   
#   png(file=paste('E:/Outputs/Cluster/nearest/weibull-gamma/ecdf_exp/',i,sp[i],'_ee','.png',sep=""))
#   x = seq(0,max(cc),0.1)
#   plot(x,pexp(x,rate=fe$estimate[1]),type="l",col="green", main="ECDF and Exponential CDF")
#   plot(ecdf(cc),add=TRUE)
#   dev.off()
  
  weibull.scale[i] = fw$estimate[2]
  weibull.shape[i] = fw$estimate[1]
  weibull.P[i] = wtest$p.value
  weibull.pvalue[i] = 1 * (wtest$p.value>0.05)
  
  lnorm.meanlog[i] = fl$estimate[1]
  lnorm.sdlog[i] = fl$estimate[2]
  lnorm.P[i] = ltest$p.value
  lnorm.pvalue[i] = 1 * (ltest$p.value>0.05)
  
#   llogis.rate[i] = fll$estimate[2]
#   llogis.shape[i] = fll$estimate[1]
#   llogis.P[i] = lltest$p.value
#   llogis.pvalue[i] = 1 * (lltest$p.value>0.05)
  
  gamma.rate[i] = fg$estimate[2]
  gamma.shape[i] = fg$estimate[1]
  gamma.P[i] = gtest$p.value
  gamma.pvalue[i] = 1 * (gtest$p.value>0.05)
  
  powerlaw.scale[i] = scalefit
  powerlaw.shape[i] = shapefit
  powerlaw.pvalue[i] = 1 * (ptest$p.value>0.05)
  powerlaw.P[i] = ptest$p.value
  #weibull.AIC[i] = AIC(fw)
  #gamma.AIC[i] = AIC(fg)
}

df = data.frame(
  sp.code = sp,
  weibull.scale = weibull.scale,
  weibull.shape = weibull.shape,
  gamma.rate = gamma.rate,
  gamma.shape = gamma.shape,
  lnorm.meanlog = lnorm.meanlog,
  lnorm.sdlog = lnorm.sdlog,
#   llogis.rate = llogis.rate,
#   llogis.shape = llogis.shape,
  powerlaw.scale = powerlaw.scale,
  powerlaw.shape = powerlaw.shape,
  weibull.P = weibull.P,
  gamma.P = gamma.P,
  lnorm.P = lnorm.P,
  # llogis.P = llogis.P,
  powerlaw.P = powerlaw.P,
  weibull.sig = weibull.pvalue,
  gamma.sig = gamma.pvalue,
  lnorm.sig = lnorm.pvalue,
  # llogis.sig = llogis.pvalue,
  powerlaw.sig = powerlaw.pvalue
)

write.csv(df, file=paste('E:/BaiduYunSych/GitHub/ClusterSize/Outputs/nearest/monoecious/fitting_parameter/Fitting_',csfile,'.csv',sep=""))

# wg = read.csv('E:/BaiduYunSych/GitHub/ClusterSize/Outputs/nearest/monoecious/fitting_parameter/weibull-gamma_all_ab50.csv',header=TRUE)
# boxplot(wg$weibull.scale, wg$gamma.scale)
# boxplot(wg$weibull.shape, wg$gamma.shape)