


csfile = 'cs_BA0'
file1 = paste('E:/Outputs/Cluster/nearest/LM_', csfile, '.csv', sep="")
lm = read.csv(file1,header=TRUE)

lm=subset(lm, lm$r_square<2)

attach(lm)
boxplot(r_square, nc_intra,	nc_inter,	dbh,	moisture,	ph,	elev,	convex,	slope,	aspect,
        names=c('r_square', 'nc_intra',	'nc_inter',	'dbh',	'moisture',	'ph',	
        'elev',	'convex',	'slope',	'aspect'),cex.axis=0.8)

boxplot(nc_intra,	nc_inter,	dbh,
        names=c('nc_intra',	'nc_inter','dbh'),cex.axis=1.8)
abline(h=0,lty=3,col='red')

boxplot(nc_intra,	nc_inter,
        names=c('nc_intra',	'nc_inter'),cex.axis=1.8)
abline(h=0,lty=3,col='red')

boxplot(r_square,names='r_square',cex.axis=1.5)
abline(h=0,lty=3,col='red')

boxplot(dbh,names='dbh',cex.axis=1.5)
abline(h=0,lty=3,col='red')

t.test(nc_intra,rnorm(64,mean=0,sd=sd(nc_intra)))
t.test(nc_inter,rnorm(64,mean=0,sd=sd(nc_inter)))
t.test(dbh,rnorm(64,mean=0,sd=sd(dbh)))

boxplot(moisture,	ph,	elev,	convex,	slope,	aspect,
        names=c( 'moisture',	'ph',	
                'elev',	'convex',	'slope',	'aspect'),cex.axis=1.8)
abline(h=0,lty=3,col='red')

abline(v=2.5, lty=1,col='red')

a = nc_intra[nc_intra>0]
length(nc_intra)

b = nc_inter[nc_inter<0]
length(nc_inter)
length(b)
