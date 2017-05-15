


csfile = 'cs'
file1 = paste('E:/Outputs/Cluster/nearest/Fitting_', csfile, '.csv', sep="")
fitting = read.csv(file1,header=TRUE)



pw = subset(fitting, fitting$weibull.sig==1)
fw = subset(fitting, fitting$weibull.sig==0)
t.test(pw$weibull.scale,fw$weibull.scale)

pg = subset(fitting, fitting$gamma.sig==1)
fg = subset(fitting, fitting$gamma.sig==0)
t.test(pg$gamma.rate,fg$gamma.rate)

pl = subset(fitting, fitting$lnorm.sig==1)
fl = subset(fitting, fitting$lnorm.sig==0)
t.test(pl$lnorm.meanlog, fl$lnorm.meanlog)

t.test(pw$weibull.shape,fw$weibull.shape)
t.test(pg$gamma.shape,fg$gamma.shape)
t.test(pl$lnorm.sdlog, fl$lnorm.sdlog)

par(mfcol=c(2,3))

boxplot(pw$weibull.scale,fw$weibull.scale,names=c('passed','failed'),cex.axis=1.8)
boxplot(pw$weibull.shape,fw$weibull.shape,names=c('passed','failed'),cex.axis=1.8)

boxplot(1/pg$gamma.rate,1/fg$gamma.rate,names=c('passed','failed'),cex.axis=1.8)
boxplot(pg$gamma.shape,fg$gamma.shape,names=c('passed','failed'),cex.axis=1.8)

boxplot(pl$lnorm.meanlog,fl$lnorm.meanlog,names=c('passed','failed'),cex.axis=1.8)
boxplot(pl$lnorm.sdlog,fl$lnorm.sdlog,names=c('passed','failed'),cex.axis=1.8)

