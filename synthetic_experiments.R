library(gbm)
library(tidyverse)
library(parallel)
library(MASS)
library(Rcpp)
Rcpp::sourceCpp('policysearchopt.cpp')

makeweights = function(data, type, pad=0.0){
  return(
    if(type=="ipw"){
      (2*data$t-1)*data$y/data$q
    } else if(type=="ipw-ovlp"){
      (2*data$t-1)*data$y*(1-data$q)
    } else if(type=="unwtd"){
      (2*data$t-1)*data$y
    } else if(type=="dr"){
      data$cate + (2*data$t-1)*data$resids/data$q
    } else if(type=="dr-ovlp"){
      data$cate*data$q*(1-data$q) + (2*data$t-1)*data$resids*(1-data$q)
    } else if(type=="direct"){
      data$cate
    } else if(type=="direct-ovlp"){
      data$cate*data$q*(1-data$q)
    } else if(type=="ovlp-pad"){
      (2*data$t-1)*data$y*(1-data$q) / (1.0 + 2.0*pad*data$q*(1-data$q))
    } else if(type=="ovlp-dr-pad"){
      data$cate/(1/(data$q*(1-data$q)) + 2.0*pad) + (2*data$t-1)*data$resids*(1-data$q) / (1.0 + 2.0*pad*data$q*(1-data$q))
    }
  )
}

addestimates.gbm <- function(data, formula){ 
  data$f0 = tryCatch({
    gbm0 = gbm::gbm(as.formula(paste("y ~ ",formula)),data=filter(data,!t),distribution='gaussian');
    predict(gbm0,data,n.trees=gbm.perf(gbm0, method = "OOB", plot.it = F),type='response')
  }, error=function(e){print(e);mean(data$y[!data$t])})
  data$f1 = tryCatch({
    gbm1 = gbm::gbm(as.formula(paste("y ~ ",formula)),data=filter(data, t),distribution='gaussian');
    predict(gbm1,data,n.trees=gbm.perf(gbm1, method = "OOB", plot.it = F),type='response')
  }, error=function(e){print(e);mean(data$y[ data$t])})
  data$cate = data$f1-data$f0
  data$resids = data$y - (data$f1*(data$t) + data$f0*(!data$t))
  data$q = tryCatch({
    gbm1 = gbm::gbm(as.formula(paste("t ~ ",formula)),data=data,distribution='bernoulli');
    predict(gbm1,data,n.trees=gbm.perf(gbm1, method = "OOB", plot.it = F),type='response')
  }, error=function(e){print(e);mean(data$t)})
  data$q = data$t*data$q + (!data$t)*(1-data$q)
  data
}

addestimates.oracle <- function(data){
  data$cate = data$cate.true
  data$resids = data$y - (data$f1.true*(data$t) + data$f0.true*(!data$t))
  data$q = data$q.true
  data
}

linearpolicy.opt <- function(data, ntheta = 20){
  policysearch(data, ntheta)
}

predict.opt <- function(pol, data) {
  data$x.1 * cos(pol[1]) + data$x.2 * sin(pol[1]) <= pol[2]
}

evalpolicy <- function(policy, data){
  mean( 
    policy*data$f1.true + 
      (!policy)*data$f0.true
  )
}

cate_fn <- function(x){
  rowSums(x)+0.25
}
makedata <- function(n,d=2,pretrans=-1,posttrans=-1,ovlp=3.5){
  x <- matrix(runif(n*d),n,d)*2-1
  if(pretrans>0){
    x <- sign(x)*abs(x)^pretrans
  }
  p <- pnorm(ovlp*x[,1])
  t <- runif(n) <= p
  cate.true <- cate_fn(x)
  f0.true <- x[,1]
  f1.true <- f0.true + cate.true
  y <- (!t)*f0.true + t*f1.true + rnorm(n)
  q.true <- t * p + (!t) * (1-p)
  if(posttrans>0){
    x <- sign(x)*abs(x)^posttrans
  }
  data.frame(x=x,p,t,y,q.true,cate.true,f0.true,f1.true)
}

mytypes = c('ipw','ipw-ovlp','unwtd','dr','dr-ovlp','direct','direct-ovlp')
do.run <- function(n,d,data_test,types,paddings,pretrans_train=F,posttrans_train=F,ovlp=ovlp) {
  best = evalpolicy(data_test$cate.true>0, data_test)
  
  formula = paste(paste("x.", 1:d, sep=""),collapse="+")
  
  df.n     = rep(0, length(types))
  df.ovlp  = rep(0, length(types))
  df.regret = rep(0.0, length(types))
  df.value = rep(0.0, length(types))
  df.type  = rep('', length(types))
  df.pad   = rep(0.0, length(types))
  df.fit   = rep('', length(types))
  i=1
  data = makedata(n,d,pretrans = pretrans_train, posttrans = posttrans_train, ovlp = ovlp)
  data = addestimates.gbm(data,formula)
  for(t in types){
   for(p in paddings){
    data$w=makeweights(data,t,p)
    data$wts0=(data$w<=0)*abs(data$w)
    data$wts1=(data$w>0)*abs(data$w)
    val = evalpolicy(predict.opt(linearpolicy.opt(data), data_test), data_test)
    df.n[i] = n
    df.ovlp[i] = ovlp
    df.value[i] = val
    df.regret[i] = best-val
    df.type[i] = t
    df.fit[i] = 'lr'
    df.pad[i] = p
    i = i+1
   }
  }
  data.frame(n=df.n,ovlp=df.ovlp,value=df.value,regret=df.regret,type=df.type,pad=df.pad,fit=df.fit)
}

do.experiment <- function(ns,d=2,ntest=100000,types=mytypes,paddings=c(0.0),pretrans_train=-1,posttrans_train=-1,pretrans_test=-1,posttrans_test=-1,overlaps=3.5,numCores=1){
  data_test <- makedata(ntest,d,pretrans = pretrans_test, posttrans = posttrans_test)
  do.call(rbind, 
    if(length(overlaps)==1)
      mclapply(ns, function(n){do.run(n,d,data_test,types=types,paddings=paddings,pretrans_train=pretrans_train,posttrans_train=posttrans_train,ovlp=overlaps)},mc.cores=numCores)
    else
      mclapply(overlaps, function(ovlp){do.run(ns,d,data_test,types=types,paddings=paddings,pretrans_train=pretrans_train,posttrans_train=posttrans_train,ovlp=ovlp)},mc.cores=numCores)
  )
}

numCores <- detectCores()

myns = round(rep(exp(seq(log(10), log(10000), (log(10000)-log(10))/20)),5000))
results_correct_stationary = do.experiment(myns,pretrans_train=-1,posttrans_train=-1,pretrans_test=-1,posttrans_test=-1,numCores = numCores) # no misspecification, no covariate shift
results_correct_shift      = do.experiment(myns,pretrans_train=-1,posttrans_train=-1,pretrans_test=2,posttrans_test=-1,numCores = numCores) # no misspecification, yes covariate shift
results_correct_shift2      = do.experiment(myns,pretrans_train=-1,posttrans_train=-1,pretrans_test=0.5,posttrans_test=-1,numCores = numCores) # no misspecification, yes covariate shift
results_wrong_stationary   = do.experiment(myns,pretrans_train=-1,posttrans_train=0.5,pretrans_test=-1,posttrans_test=0.5,numCores = numCores) # yes misspecification, no covariate shift
results_wrong_shift        = do.experiment(myns,pretrans_train=-1,posttrans_train=0.5,pretrans_test=2,posttrans_test=0.5,numCores = numCores) # yes misspecification, yes covariate shift
results_wrong_shift2        = do.experiment(myns,pretrans_train=-1,posttrans_train=0.5,pretrans_test=0.5,posttrans_test=0.5,numCores = numCores) # yes misspecification, yes covariate shift

namemutate <- function(res) {
  res %>% mutate(
    Method = res$type%>%str_replace_all("direct","DM")%>%str_replace_all("ipw","IPW")%>%str_replace_all("dr","DR")%>%str_replace_all("unwtd","Unweighted")%>%str_replace_all("-ovlp",""),
    Retargeted = res$type%>%str_replace_all("(direct|ipw|dr|unwtd)-?","")=="ovlp"
  )
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

plot_correct_stationary = results_correct_stationary%>%namemutate%>%group_by(n,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Well-specified, Stationary") + aes(n,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(2,6))  + ylab("Regret")
plot_correct_shift = results_correct_shift%>%namemutate%>%group_by(n,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Well-specified, Shift inward") + aes(n,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(2,6))  + ylab("Regret")
plot_correct_shift2 = results_correct_shift2%>%namemutate%>%group_by(n,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Well-specified, Shift outward") + aes(n,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(2,6))  + ylab("Regret")
plot_wrong_stationary = results_wrong_stationary%>%namemutate%>%group_by(n,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Misspecified, Stationary") + aes(n,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(2,6))  + ylab("Regret")
plot_wrong_shift = results_wrong_shift%>%namemutate%>%group_by(n,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Misspecified, Shift inward") + aes(n,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(2,6))  + ylab("Regret")
plot_wrong_shift2 = results_wrong_shift2%>%namemutate%>%group_by(n,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Misspecified, Shift outward") + aes(n,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_y_log10() + scale_x_log10() + scale_shape_manual(values=c(2,6))  + ylab("Regret")

ggsave('legend.pdf',plot=g_legend(plot_correct_stationary))
ggsave('plot_correct_stationary.pdf',plot=plot_correct_stationary + theme(legend.position="none",aspect.ratio=1.25))
ggsave('plot_correct_shift.pdf',plot=plot_correct_shift + theme(legend.position="none",aspect.ratio=1.25))
ggsave('plot_correct_shift2.pdf',plot=plot_correct_shift2 + theme(legend.position="none",aspect.ratio=1.25))
ggsave('plot_wrong_stationary.pdf',plot=plot_wrong_stationary + theme(legend.position="none",aspect.ratio=1.25))
ggsave('plot_wrong_shift.pdf',plot=plot_wrong_shift + theme(legend.position="none",aspect.ratio=1.25))
ggsave('plot_wrong_shift2.pdf',plot=plot_wrong_shift2 + theme(legend.position="none",aspect.ratio=1.25))
 
paddings = 2^seq(-12,12,0.75)
results_correct_stationary_pad = do.experiment(rep(2000,50000),types=c('ovlp-dr-pad'),paddings=paddings,pretrans_train=-1,posttrans_train=-1,pretrans_test=-1,posttrans_test=-1,numCores = numCores)
results_correct_shift_pad      = do.experiment(rep(2000,50000),types=c('ovlp-dr-pad'),paddings=paddings,pretrans_train=-1,posttrans_train=-1,pretrans_test=2,posttrans_test=-1,numCores = numCores)
results_correct_shift2_pad      = do.experiment(rep(2000,50000),types=c('ovlp-dr-pad'),paddings=paddings,pretrans_train=-1,posttrans_train=-1,pretrans_test=0.5,posttrans_test=-1,numCores = numCores)
results_wrong_stationary_pad   = do.experiment(rep(2000,50000),types=c('ovlp-dr-pad'),paddings=paddings,pretrans_train=-1,posttrans_train=0.5,pretrans_test=-1,posttrans_test=0.5,numCores = numCores)
results_wrong_shift_pad        = do.experiment(rep(2000,50000),types=c('ovlp-dr-pad'),paddings=paddings,pretrans_train=-1,posttrans_train=0.5,pretrans_test=2,posttrans_test=0.5,numCores = numCores)
results_wrong_shift2_pad        = do.experiment(rep(2000,50000),types=c('ovlp-dr-pad'),paddings=paddings,pretrans_train=-1,posttrans_train=0.5,pretrans_test=0.5,posttrans_test=0.5,numCores = numCores)

plot_correct_stationary_pad = results_correct_stationary_pad%>%filter(type=='ovlp-dr-pad')%>%mutate(c=1/(2*pad))%>%group_by(c,type,n)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Well-specified, Stationary") + aes(c,regret.mean,color=type) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_x_log10() + ylab("Regret")
plot_correct_shift_pad = results_correct_shift_pad%>%filter(type=='ovlp-dr-pad')%>%mutate(c=1/(2*pad))%>%group_by(c,type,n)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Well-specified, Shift inward") + aes(c,regret.mean,color=type) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_x_log10() + ylab("Regret")
plot_correct_shift2_pad = results_correct_shift_pad%>%filter(type=='ovlp-dr-pad')%>%mutate(c=1/(2*pad))%>%group_by(c,type,n)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Well-specified, Shift outward") + aes(c,regret.mean,color=type) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_x_log10() + ylab("Regret")
plot_wrong_stationary_pad = results_wrong_stationary_pad%>%filter(type=='ovlp-dr-pad')%>%mutate(c=1/(2*pad))%>%group_by(c,type,n)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Misspecified, Stationary") + aes(c,regret.mean,color=type) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_x_log10() + ylab("Regret")
plot_wrong_shift_pad = results_wrong_shift_pad%>%filter(type=='ovlp-dr-pad')%>%mutate(c=1/(2*pad))%>%group_by(c,type,n)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Misspecified, Shift inward") + aes(c,regret.mean,color=type) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_x_log10() + ylab("Regret")
plot_wrong_shift2_pad = results_wrong_shift2_pad%>%filter(type=='ovlp-dr-pad')%>%mutate(c=1/(2*pad))%>%group_by(c,type,n)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Misspecified, Shift outward") + aes(c,regret.mean,color=type) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + scale_x_log10() + ylab("Regret")

ggsave('plot_correct_stationary_pad.pdf',plot=plot_correct_stationary_pad + theme(legend.position="none",aspect.ratio=1))
ggsave('plot_correct_shift_pad.pdf',plot=plot_correct_shift_pad + theme(legend.position="none",aspect.ratio=1))
ggsave('plot_correct_shift2_pad.pdf',plot=plot_correct_shift2_pad + theme(legend.position="none",aspect.ratio=1))
ggsave('plot_wrong_stationary_pad.pdf',plot=plot_wrong_stationary_pad + theme(legend.position="none",aspect.ratio=1))
ggsave('plot_wrong_shift_pad.pdf',plot=plot_wrong_shift_pad + theme(legend.position="none",aspect.ratio=1))
ggsave('plot_wrong_shift2_pad.pdf',plot=plot_wrong_shift2_pad + theme(legend.position="none",aspect.ratio=1))

overlaps = 2^seq(-10,5,.5)
results_correct_stationary_ovlp  = do.experiment(10000,overlaps=rep(overlaps,5000),pretrans_train=-1,posttrans_train=-1,pretrans_test=-1,posttrans_test=-1,numCores = numCores)
results_correct_shift_ovlp       = do.experiment(10000,overlaps=rep(overlaps,5000),pretrans_train=-1,posttrans_train=-1,pretrans_test=2,posttrans_test=-1,numCores = numCores)
results_correct_shift2_ovlp      = do.experiment(10000,overlaps=rep(overlaps,5000),pretrans_train=-1,posttrans_train=-1,pretrans_test=0.5,posttrans_test=-1,numCores = numCores)
results_wrong_stationary_ovlp    = do.experiment(10000,overlaps=rep(overlaps,5000),pretrans_train=-1,posttrans_train=0.5,pretrans_test=-1,posttrans_test=0.5,numCores = numCores)
results_wrong_shift_ovlp         = do.experiment(10000,overlaps=rep(overlaps,5000),pretrans_train=-1,posttrans_train=0.5,pretrans_test=2,posttrans_test=0.5,numCores = numCores)
results_wrong_shift2_ovlp        = do.experiment(10000,overlaps=rep(overlaps,5000),pretrans_train=-1,posttrans_train=0.5,pretrans_test=0.5,posttrans_test=0.5,numCores = numCores)

plot_correct_stationary_ovlp = results_correct_stationary_ovlp%>%namemutate%>%group_by(ovlp,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Well-specified, Stationary") + aes(ovlp,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + xlab(expression(beta))  + ylab("Regret") + scale_x_log10() + scale_shape_manual(values=c(2,6))
plot_correct_shift_ovlp = results_correct_shift_ovlp%>%namemutate%>%group_by(ovlp,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Well-specified, Shift inward") + aes(ovlp,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + xlab(expression(beta))  + ylab("Regret") + scale_x_log10() + scale_shape_manual(values=c(2,6))
plot_correct_shift2_ovlp = results_correct_shift2_ovlp%>%namemutate%>%group_by(ovlp,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Well-specified, Shift outward") + aes(ovlp,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + xlab(expression(beta))  + ylab("Regret") + scale_x_log10() + scale_shape_manual(values=c(2,6))
plot_wrong_stationary_ovlp = results_wrong_stationary_ovlp%>%namemutate%>%group_by(ovlp,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Misspecified, Stationary") + aes(ovlp,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + xlab(expression(beta))  + ylab("Regret") + scale_x_log10() + scale_shape_manual(values=c(2,6))
plot_wrong_shift_ovlp = results_wrong_shift_ovlp%>%namemutate%>%group_by(ovlp,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Misspecified, Shift inward") + aes(ovlp,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + xlab(expression(beta))  + ylab("Regret") + scale_x_log10() + scale_shape_manual(values=c(2,6))
plot_wrong_shift2_ovlp = results_wrong_shift2_ovlp%>%namemutate%>%group_by(ovlp,Method,Retargeted)%>%summarise(regret.mean=mean(regret),regret.se=sd(regret)/sqrt(n())) %>% ggplot() + ggtitle("Misspecified, Shift outward") + aes(ovlp,regret.mean,color=Method,shape=Retargeted,linetype=Retargeted) + geom_point() + geom_line() + geom_ribbon(aes(ymin=regret.mean-1.96*regret.se,ymax=regret.mean+1.96*regret.se),alpha=0.15,linetype=0) + xlab(expression(beta))  + ylab("Regret") + scale_x_log10() + scale_shape_manual(values=c(2,6))

ggsave('legend_ovlp.pdf',plot=g_legend(plot_correct_stationary_ovlp))
ggsave('plot_correct_stationary_ovlp.pdf',plot=plot_correct_stationary_ovlp + theme(legend.position="none",aspect.ratio=1.25))
ggsave('plot_correct_shift_ovlp.pdf',plot=plot_correct_shift_ovlp + theme(legend.position="none",aspect.ratio=1.25))
ggsave('plot_correct_shift2_ovlp.pdf',plot=plot_correct_shift2_ovlp + theme(legend.position="none",aspect.ratio=1.25))
ggsave('plot_wrong_stationary_ovlp.pdf',plot=plot_wrong_stationary_ovlp + theme(legend.position="none",aspect.ratio=1.25))
ggsave('plot_wrong_shift_ovlp.pdf',plot=plot_wrong_shift_ovlp + theme(legend.position="none",aspect.ratio=1.25))
ggsave('plot_wrong_shift2_ovlp.pdf',plot=plot_wrong_shift2_ovlp + theme(legend.position="none",aspect.ratio=1.25))
