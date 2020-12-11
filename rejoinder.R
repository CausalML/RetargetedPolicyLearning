library(gbm)
library(tidyverse)
library(parallel)

library(Rcpp)
Rcpp::sourceCpp('policysearchopt_interval.cpp')

addestimates.oracle <- function(data){
  data$cate = data$cate.true
  data$resids = data$y - (data$f1.true*(data$t) + data$f0.true*(!data$t))
  data$q = data$q.true
  data
}

intervalpolicy.opt <- function(data){
  policysearch(data)
}

predict.opt <- function(pol, data) {
  data$x <= pol[1]
}

evalpolicy <- function(policy, data){
  mean( 
    policy*data$f1.true + 
      (!policy)*data$f0.true
  )
}

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
    } else if(type=="dr-ovlp-LLL"){
      (data$cate*data$q*(1-data$q) + (2*data$t-1)*data$resids*(1-data$q)) * abs(data$cate)
    } else if(type=="dr-ovlp-Margin"){
      (data$cate*data$q*(1-data$q) + (2*data$t-1)*data$resids*(1-data$q)) / abs(data$cate)
    } else if(type=="dr-ovlp-LLL2"){
      (data$cate*data$q*(1-data$q) + (2*data$t-1)*data$resids*(1-data$q)) * (data$cate^2)
    } else if(type=="dr-ovlp-Margin2"){
      (data$cate*data$q*(1-data$q) + (2*data$t-1)*data$resids*(1-data$q)) / (data$cate^2)
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

make.cvgroup = function(n, K, right = TRUE) {
  split     = runif(n)
  return(as.numeric(cut(split, quantile(split, probs = seq(0, 1, 1/K)), include.lowest = TRUE, right = right)))
}

make.cvgroup.balanced = function(data, K) {
  cvgroup = numeric(nrow(data))
  cvgroup[data$t==1] = make.cvgroup(sum(data$t==1), K, right = TRUE)
  cvgroup[data$t==0] = make.cvgroup(sum(data$t==0), K, right = FALSE)
  return(cvgroup)
}

addestimates.gbm <- function(data, formula, K = 2){
  cvgroup   = make.cvgroup.balanced(data, K)
  f0 = numeric(nrow(data))
  f1 = numeric(nrow(data))
  q  = numeric(nrow(data))
  for (k in 1:K) {
    f0[cvgroup==k] = tryCatch({
      gbm0 = gbm::gbm(formula=as.formula(paste("y ~ ",formula)),data=filter(data[cvgroup!=k,],!t),distribution='gaussian');
      predict(gbm0,data[cvgroup==k,],n.trees=gbm.perf(gbm0, method = "OOB", plot.it = F),type='response')
    }, error=function(e){print(e);mean(data$y[!data$t])})
    f1[cvgroup==k] = tryCatch({
      gbm1 = gbm::gbm(formula=as.formula(paste("y ~ ",formula)),data=filter(data[cvgroup!=k,], t),distribution='gaussian');
      predict(gbm1,data[cvgroup==k,],n.trees=gbm.perf(gbm1, method = "OOB", plot.it = F),type='response')
    }, error=function(e){print(e);mean(data$y[ data$t])})
    q[cvgroup==k] = tryCatch({
      gbm1 = gbm::gbm(formula=as.formula(paste("t ~ ",formula)),data=data[cvgroup!=k,],distribution='bernoulli');
      predict(gbm1,data[cvgroup==k,],n.trees=gbm.perf(gbm1, method = "OOB", plot.it = F),type='response')
    }, error=function(e){print(e);mean(data$t)})
  }
  data$f0 = f0
  data$f1 = f1
  data$q  = q
  data$cate = data$f1-data$f0
  data$resids = data$y - (data$f1*(data$t) + data$f0*(!data$t))
  data$q = data$t*data$q + (!data$t)*(1-data$q)
  data
}


baseline = function(x){ x }
cate = list()
cate[[1]] = function(x){
  (x-0.5)*(x>0)
}
cate[[3]] = function(x){
  x
}
cate[[4]] = function(x){
  (x-0.3)*(x>-0.4)
}
prop = list()
prop[[1]] = function(x){
  0.5*(x<=0)+pnorm(3.5*x)*(x>0)
}
prop[[3]] = function(x){
  pnorm(3.5*x)
}
prop[[4]] = function(x){
  0.5*(x<=-0.4)+pnorm(2.5*x)*(x>-.04)
}

makedata <- function(n, case=1){
  x <- runif(n)*2-1
  p <- prop[[case]](x)
  t <- runif(n) <= p
  cate.true <- cate[[case]](x)
  f0.true <- baseline(x)
  f1.true <- f0.true + cate.true
  y <- (!t)*f0.true + t*f1.true + rnorm(n)
  q.true <- t * p + (!t) * (1-p)
  x <- -x # flip because my policy is I[x<=theta]
  data.frame(x=x,p,t,y,q.true,cate.true,f0.true,f1.true)
}

mytypes = c('dr','dr-ovlp','dr-ovlp-LLL','dr-ovlp-Margin','dr-ovlp-LLL2','dr-ovlp-Margin2')

do.run <- function(n,data_test,types,case=1) {
 tryCatch({
  best = evalpolicy(data_test$cate.true>0, data_test)
  
  formula = "x"
  
  df.n     = rep(0, length(types))
  df.regret = rep(0.0, length(types))
  df.value = rep(0.0, length(types))
  df.type  = rep('', length(types))
  i=1
  data = makedata(n,case)
  data = addestimates.gbm(data,formula)
  for(t in types){
    data$w=makeweights(data,t)
    data$wts0=(data$w<=0)*abs(data$w)
    data$wts1=(data$w>0)*abs(data$w)
    val = evalpolicy(predict.opt(intervalpolicy.opt(data), data_test), data_test)
    df.n[i] = n
    df.value[i] = val
    df.regret[i] = best-val
    df.type[i] = t
    i = i+1
  }
  data.frame(n=df.n,value=df.value,regret=df.regret,type=df.type)
 }, error = function(e) {
  data.frame()
 })
}

results = list()

ntest=100000
numCores <- detectCores()

for(c in c(1,3,4)) {
  data_test <- makedata(ntest,c)
  res = lapply(1:500, function(j){do.call(rbind, mclapply(1:10, function(i){do.run(2000,data_test,mytypes,c)}, mc.cores=numCores))})
  res = do.call(rbind, res)
  res$case = c
  results[[c]] = res
}

results = do.call(rbind, results)

library(xtable)
print(xtable(
  results %>% group_by(type, case) %>% summarise(reg=mean(regret), se=sd(regret), n=n()) %>% mutate(reg = sprintf('%.3f (%.3f)',reg,se)) %>% select(type,case,reg) %>% pivot_wider(names_from = type, values_from = reg)
  ,include.rownames=F))


