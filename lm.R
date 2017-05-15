library(MASS)
library(ggplot2)
csfile = 'cs_BA0'
file1 = paste('E:/BaiduYunSych/GitHub/ClusterSize/Outputs/nearest/all/cs_scanner', csfile, '.csv', sep="")
cs = read.csv(file1,header=TRUE)
#sp = c('MACBRE','QUELIT','LITACU',
#       'DISRAC','ARTSTY','DIOMOR',
#       'SYMLAU','ORMPAC','CINPAR')
n = 67
sp = names(cs)[1:n+1]

r_square = c(1:n)
nc_intra = c(1:n)
nc_inter = c(1:n)
dbh = c(1:n)
moisture = c(1:n)
ph = c(1:n)
elev = c(1:n)
convex = c(1:n)
slope = c(1:n)
aspect = c(1:n)


for(i in c(1:n)){
  species = sp[i]
  #species = 'QUELIT'
  file2 = paste('E:/BaiduYunSych/GitHub/ClusterSize/Outputs/nearest/all/ncdata/', species, '.csv', sep="")
  mydata = read.csv(file2, header=TRUE)
  sp.dat = na.omit(mydata[, c(7,4,8,9,13:32)])
  cc = eval(parse(text = paste('cs$',species,sep="")))
  cc = cc[!is.na(cc)]
  cc = cc + 1
  sp.dat[,1] = cc

  #sp.dat = subset(sp.dat, sp.dat[,1]>0)
  
  cor(sp.dat)[cor(sp.dat)>0.9]
  
  #sp.dat$cluster.size = log(sp.dat$cluster.size)
  sp.dat$dbh = log(sp.dat$dbh,2)
  sp.dat$nc_intra = sqrt(sp.dat$nc_intra)+1
  sp.dat$nc_inter = log(sp.dat$nc_inter)
  sp.dat$ph = log(sp.dat$ph)
#   sp.dat$no3 = log(sp.dat$no3)
#   sp.dat$TMn = log(sp.dat$TMn)
#   sp.dat$TCu = log(sp.dat$TCu)
#   sp.dat$TMg = log(sp.dat$TMg)
#   sp.dat$AP = log(sp.dat$AP)
#   sp.dat$AZn = log(sp.dat$AZn)
#   sp.dat$AB = log(sp.dat$AB)
#   sp.dat$ACu = log(sp.dat$ACu)
  #sp.dat$nc_rate = log(mydata$nc_intra/(mydata$nc_total-mydata$nc_intra))
  
  sp.dat1 = scale(sp.dat[,2:24], center=TRUE, scale=TRUE)
  sp.dat1 = cbind(sp.dat[1], sp.dat1)
  
#   par(mfrow=c(4,2))
#   for(i in c(1:8)){
#     x = colnames(sp.dat1)[i]
#     hist(sp.dat1[,i],xlab=x,main=paste('Histogram of',x)) 
#   }
#   
#   par(mfrow=c(4,2))
#   for(i in 17:24){
#     plot(sp.dat1[,i],sp.dat1$cluster.size,cex=.2,xlab=colnames(sp.dat1)[i])
#     lines(lowess(sp.dat1[,i],sp.dat1$cluster.size),col=2) 
#   }
#   
#   
#   f = lm(cluster.size ~ dbh+nc_intra+nc_total+moisure+ph+no3+TP+TB+TMn+TCu+TMg
#           +TAl+TK+AP+AZn+AB+AMn+ACu+AFe+elev+convex+slope+aspect,data=sp.dat1)
#   f.best = step(f)
#   summary(f.best)
  
#   fe = glm(cluster.size ~ dbh+nc_intra+nc_inter+moisure+ph+
#              elev+convex+slope+aspect,data=sp.dat1, family=Gamma(link='log'))
#   fe.best = step(fe)
#   sfe = summary(fe.best)
#   r2 = 1-(sfe$deviance/sfe$null.deviance)
  
  fe = lm(cluster.size ~ dbh+nc_intra+nc_inter+moisure+ph+
             elev+convex+slope+aspect,data=sp.dat1)
  fe.best = step(fe)
  sfe = summary(fe.best)
  r2=sfe$r.squared
  
  if (dim(sfe$coef)[1] > 1) {
  
    eff = data.frame(sfe$coef)
    eff = eff[-1,]
    factors = row.names(eff)
    
    efft = data.frame(t(eff))
    r_square[i] = round(r2,2)
    if ('nc_intra' %in% factors) {nc_intra[i] = efft$nc_intra[1]}
    else {nc_intra[i] = 0}
    if ('nc_inter' %in% factors) {nc_inter[i] = efft$nc_inter[1]}
    else {nc_inter[i] = 0}
    if ('moisure' %in% factors) {
      factors = gsub("moisure","moisture",factors)
      moisture[i] = efft$moisure[1]}
    else {moisture[i] = 0}
    if ('dbh' %in% factors) {dbh[i] = efft$dbh[1]}
    else {dbh[i] = 0}
    if ('ph' %in% factors) {ph[i] = efft$ph[1]}
    else {ph[i] = 0}
    if ('elev' %in% factors) {elev[i] = efft$elev[1]}
    else {elev[i] = 0}
    if ('convex' %in% factors) {convex[i] = efft$convex[1]}
    else {convex[i] = 0}
    if ('slope' %in% factors) {slope[i] = efft$slope[1]}
    else {slope[i] = 0}
    if ('aspect' %in% factors) {aspect[i] = efft$aspect[1]}
    else {aspect[i] = 0}
    
  
    gp1 = ggplot(eff, aes(x = 1:(dim(sfe$coef)[1]-1), y = eff$Estimate, label = factors)) + geom_point(size = 3) 
    gp1 = gp1 + geom_errorbar(aes(ymax=eff$Estimate+1.96*eff$Std..Error, ymin=eff$Estimate-1.96*eff$Std..Error), width = 0.2)
    gp1 = gp1 + geom_text(vjust = -0.5, col = "red", size=10) 
    gp1 = gp1 + geom_hline(yintercept = 0, linetype = 2)
    gp1 = gp1 + xlab("") + ylab(paste("Effect size and 95% confidence interval",' ("Rsquare=',round(r2,2),')',sep=""))#McFadden\'s pseudo-R2
    gp1 = gp1 + theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(),
                                   axis.text.x = element_text(size=20),axis.title.x =element_text(size=23,vjust=-5),
                                   plot.margin = unit(c(1,1,1,1), "cm"))
    gp1 = gp1 + coord_flip()
    ggsave(file=paste('E:/BaiduYunSych/GitHub/ClusterSize/Outputs/nearest/monoecious/effect_size_fig/LM_',csfile,'/',i,sp[i],'effectsize','.png',sep=""), gp1)
  }
}

df = data.frame(
  sp.code = sp,
  r_square = r_square,
  nc_intra = nc_intra,
  nc_inter = nc_inter,
  dbh = dbh,
  moisture = moisture,
  ph = ph,
  elev = elev,
  convex = convex,
  slope = slope,
  aspect = aspect
)

write.csv(df, file=paste('E:/BaiduYunSych/GitHub/ClusterSize/Outputs/nearest/monoecious/effect_size_numeric/LM_',csfile,'.csv',sep=""))

#for(i in c(1:32)){a[i]=sum(is.nan(paste(text='sp.dat[,',i,']',sep="")))}
  