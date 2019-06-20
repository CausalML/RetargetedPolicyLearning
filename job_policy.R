library(tidyverse)
library(gbm)
library(parallel)
library(MASS)
library(Rcpp)
Rcpp::sourceCpp('optlinear.cpp')

optlinearpol.fit <- function(datain, Xform, timelimit=60, verbose=1) {
  X=model.matrix(as.formula(paste('~',Xform)), data=datain)
  colnames(X)=c()
  data = data.frame(wts1=datain$wts1, 
                    wts2=datain$wts2, 
                    wts3=datain$wts3,
                    x=X)
  d = ncol(X)
  lpol = optlinearpol(data, nrow(data), d, timelimit, verbose)
  out = list()
  out$a1 = lpol[1]
  out$b1 = lpol[2:(d+1)]
  out$a2 = lpol[d+2]
  out$b2 = lpol[(d+3):(2*d+2)]
  out$a3 = lpol[2*d+3]
  out$b3 = lpol[(2*d+4):(3*d+3)]
  out$success = lpol[3*d+4]
  out
}

optlinearpol.pred <- function(pol, data, Xform) {
  X=model.matrix(as.formula(paste('~',Xform)), data=data)
  mapply(function(x,y,z){which.max(c(x,y,z))}, pol$a1 + X %*% pol$b1, pol$a2 + X %*% pol$b2, pol$a3 + X %*% pol$b3)
}

job = read.csv('behaghel.csv')
# make ipw weights
job = job %>% mutate(ipw = POIDS_PZ_6MOIS / (CVE*mean(job$CVE)+OPP*mean(job$OPP)+CLA*mean(job$CLA)))
### subtract treatment means
CVEeff = mean(job$ipw * job$EMPLOI_6MOIS * job$CVE) - mean(job$ipw * job$EMPLOI_6MOIS * job$CLA)
OPPeff = mean(job$ipw * job$EMPLOI_6MOIS * job$OPP) - mean(job$ipw * job$EMPLOI_6MOIS * job$CLA)
CLAbase = mean(job$ipw * job$EMPLOI_6MOIS * job$CLA)
job = job %>% mutate(y=EMPLOI_6MOIS - CVE*CVEeff - OPP*OPPeff - CLAbase)

# make personalization X
Xcat = c('temps', 'rsqstat', 'zus',  'College_education', 'One_to_5_years_of_exp_in_the_job',  'Technician', 'Skilled_clerical_worker', 'Skilled_blue_colar', 'Q1', 'Q2', 'Q3', 'Q4')
for(xc in Xcat){
  job[,xc] = factor(job[,xc])
}
Xform = paste(Xcat,collapse=" + ")

job$t = 1*job$CVE + 2*job$OPP + 3*job$CLA

makeweights <- function(data,type,retarget) {
  if(type == 'DM') {
    data$wts1 = data$f1
    data$wts2 = data$f2
    data$wts3 = data$f3
  } else if(type == 'IPW') {
    data$wts1 = (data$t==1)*(data$y)/data$q
    data$wts2 = (data$t==2)*(data$y)/data$q
    data$wts3 = (data$t==3)*(data$y)/data$q
  } else if(type == 'DR') {
    data$wts1 = data$f1 + (data$t==1)*(data$y - data$f1)/data$q
    data$wts2 = data$f2 + (data$t==2)*(data$y - data$f2)/data$q
    data$wts3 = data$f3 + (data$t==3)*(data$y - data$f3)/data$q
  }
  if(retarget==T){
    data$wts1 = data$wts1*data$ovlpw
    data$wts2 = data$wts2*data$ovlpw
    data$wts3 = data$wts3*data$ovlpw
  }
  data
}

do.run <- function(o) {
  trnidx = sample(nrow(job),floor(nrow(job)/5))
  jobtst = job[-trnidx,]
  jobtrn = job[trnidx,]
  
  kp = (2/3)/o
  
  jobtrn$z.keep = ((jobtrn$One_to_5_years_of_exp_in_the_job==1)-mean(jobtrn$One_to_5_years_of_exp_in_the_job==1))/sd(jobtrn$One_to_5_years_of_exp_in_the_job==1) + ((jobtrn$Skilled_blue_colar==1)-mean(jobtrn$Skilled_blue_colar==1))/sd(jobtrn$Skilled_blue_colar==1) +  ((jobtrn$rsqstat=='RS2')-mean(jobtrn$rsqstat=='RS2'))/sd(jobtrn$rsqstat=='RS2')
  jobtrn$p.keep =(
    (jobtrn$z.keep<quantile(jobtrn$z.keep,1/3))*(kp/2+jobtrn$CLA*(1-3*kp/2))+
    (jobtrn$z.keep>=quantile(jobtrn$z.keep,2/3))*(kp/2+jobtrn$CVE*(1-3*kp/2))+
    ((jobtrn$z.keep>=quantile(jobtrn$z.keep,1/3)) & (jobtrn$z.keep<quantile(jobtrn$z.keep,2/3)))*(kp/2+jobtrn$OPP*(1-3*kp/2))
  )
  jobtrn$keep = runif(nrow(jobtrn)) <= jobtrn$p.keep
  jobtrn = jobtrn[jobtrn$keep,]
  
  gbmPS = gbm(as.formula(paste("t ~ One_to_5_years_of_exp_in_the_job + Skilled_blue_colar + rsqstat + ipw")),data=jobtrn,distribution='multinomial') #propensity
  ps = predict(gbmPS,jobtrn,n.trees=gbm.perf(gbmPS, method = "OOB", plot.it = F),type='response')
  jobtrn$q = rowSums(ps*jobtrn[,c('CVE','OPP','CLA')])
  jobtrn$suminvps = rowSums(1/ps)
  
  gbmCVE = gbm(as.formula(paste("y ~ ",Xform)),data=filter(jobtrn,CVE==1),distribution='gaussian') #public
  gbmOPP = gbm(as.formula(paste("y ~ ",Xform)),data=filter(jobtrn,OPP==1),distribution='gaussian') #private
  gbmCLA = gbm(as.formula(paste("y ~ ",Xform)),data=filter(jobtrn,CLA==1),distribution='gaussian') #standard
  
  jobtrn$f1 = predict(gbmCVE,jobtrn,n.trees=gbm.perf(gbmCVE, method = "OOB", plot.it = F),type='response')
  jobtrn$f2 = predict(gbmOPP,jobtrn,n.trees=gbm.perf(gbmOPP, method = "OOB", plot.it = F),type='response')
  jobtrn$f3 = predict(gbmCLA,jobtrn,n.trees=gbm.perf(gbmCLA, method = "OOB", plot.it = F),type='response')
  
  jobtrn$ovlpw = 1/(jobtrn$suminvps+1/2)
  jobtrn$ovlpw = jobtrn$ovlpw/mean(jobtrn$ovlpw)
  
  m = 10
  df.value = rep(0.0, m)
  df.type  = rep('', m)
  df.retarget  = rep('', m)
  df.success  = rep('', m)
  df.ovlp  = rep('', m)
  i=1
  for(t in c('DM','DR','IPW')) {
    for(retarget in c(F,T)) {
      jobtrn = makeweights(jobtrn,t,retarget)
      
      jobtrnfold <- jobtrn%>%group_by_at(Xcat)%>%summarise(wts1=sum(wts1),wts2=sum(wts2),wts3=sum(wts3))
      pol = optlinearpol.fit(jobtrnfold, Xform, 300, 0)
      
      df.value[i] = mean(jobtst$ipw * jobtst$y * (jobtst$t==optlinearpol.pred(pol, jobtst, Xform)))
      df.type[i] = t
      df.retarget[i] = retarget
      df.ovlp[i] = o
      df.success[i] = pol$success
      i = i+1
    }
  }
  
  df.value[i] = mean(jobtst$ipw * jobtst$y * (jobtst$t==1))
  df.type[i]  = 'all1'
  df.retarget[i] = F
  df.ovlp[i] = o
  df.success[i] = 1
  i = i+1
  
  df.value[i] = mean(jobtst$ipw * jobtst$y * (jobtst$t==2))
  df.type[i]  = 'all2'
  df.retarget[i] = F
  df.ovlp[i] = o
  df.success[i] = 1
  i = i+1
  
  df.value[i] = mean(jobtst$ipw * jobtst$y * (jobtst$t==3))
  df.type[i]  = 'all3'
  df.retarget[i] = F
  df.ovlp[i] = o
  df.success[i] = 1
  i = i+1
  
  predOPP = predict(gbmOPP,jobtst,n.trees=gbm.perf(gbmOPP, method = "OOB", plot.it = F),type='response')
  predCVE = predict(gbmCVE,jobtst,n.trees=gbm.perf(gbmCVE, method = "OOB", plot.it = F),type='response')
  predCLA = predict(gbmCLA,jobtst,n.trees=gbm.perf(gbmCLA, method = "OOB", plot.it = F),type='response')
  rec = mapply(function(x,y,z){which.max(c(x,y,z))}, predCVE, predOPP, predCLA)
  
  df.value[i] = mean(jobtst$ipw * jobtst$y * (jobtst$t==rec))
  df.type[i]  = 'DC'
  df.retarget[i] = F
  df.ovlp[i] = o
  df.success[i] = 1
    
  data.frame(value=df.value,type=df.type,retarget=df.retarget,ovlp=df.ovlp,success=df.success)
}

numCores <- detectCores()
res = mclapply(rep(4,1440),do.run,mc.cores=numCores)
