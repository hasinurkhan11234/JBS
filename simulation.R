######################################
#######simulation#####################
######################################

  beta <- c(0.50,1.15)
  Z <- matrix(rnorm(n*2), ncol = 2) #covariates
  u <- runif(n) #for inverse transform simulation
  #time simulation based on inverse CDF transform method - refer to paper
  alpha <- 0.75
  tau <- c(2.50,6)
  lambda <-c(0.05,0.35,0.55)
  #two change points
  time <- ifelse(-log(u) <= (lambda[1]*exp(Z%*%beta)*tau[1]^alpha), (-log(u)/(lambda[1]*exp(Z%*%beta)))^(1/alpha),
                 ifelse(-log(u) > ((lambda[1]*exp(Z%*%beta)*tau[1]^alpha) + (lambda[2]*exp(Z%*%beta)*(tau[2]^alpha - tau[1]^alpha))),
                        ((-log(u) - (lambda[1]*exp(Z%*%beta)*tau[1]^alpha) -
                            (lambda[2]*exp(Z%*%beta)*(tau[2]^alpha - tau[1]^alpha)) + (lambda[3]*exp(Z%*%beta)*tau[2]^alpha))/(lambda[3]*exp(Z%*%beta)))^(1/alpha),
                        ((-log(u) - (lambda[1]*exp(Z%*%beta)*tau[1]^alpha) + (lambda[2]*exp(Z%*%beta)*tau[1]^alpha))/(lambda[2]*exp(Z%*%beta)))^(1/alpha))
  )
  
time<-((-log(u) - (lambda[1]*exp(Z%*%beta)*tau[1]^alpha) -
                            (lambda[2]*exp(Z%*%beta)*(tau[2]^alpha - tau[1]^alpha)) + (lambda[3]*exp(Z%*%beta)*tau[2]^alpha))/(lambda[3]*exp(Z%*%beta)))^(1/alpha)  
  n=500
  cen_time<-NULL
  if(censor == 20){
    cen_time<- runif(n,0,quantile(time,0.947))
  }
  t1<-pmin(time,cen_time)
  status <- ifelse(time <= cen_time, 1, 0)
  prop.table(table(status))
  delta <- ifelse(time <= cen_time, 1, 0)
  
  delta <- rep(1, times = n) #censoring 'indicator' for no censoring option
  
  perc.cens <-length(which(delta == 0))/length(delta) #percent censored
  
  # time <- ifelse(time0 < censor, time0, censor) #time is min(censor, time0)
  dat <- data.frame(t1, delta)
  dat1<-data.frame(ceiling(t1),delta)   # to make the time round figure
  
  # applying survfit command
  s=survfit(Surv(ceiling(t1),delta)~1,data=dat1)
  s1=summary(s)
  ts.plot(s1$n.event)
Dt=s1$n.event

At <- 2*sqrt( Dt + 3/8)
ts.plot(At)
nsp_poly(At,deg=1) -> At.n

rr<-(At[1:7])
mn.rr<-mean(rr)
rr1<-as.list(rr)
rr1[[1]]
S <- list()
S[[1]] <- 0
for(i in 1:7){
  S[[i+1]]<-S[[i]] + (rr1[[i]]-mn.rr)  
}
rr[[34]]
s1=S %>% unlist()
max.S<-max(s1)
S.diff<-max.S-min.S
c=which(s1==max.S)
# so we got our first cp at 3

nsp_poly(At[4:28],deg=1) -> At1.n
rr2<-(At[4:8])
mn.rr2<-mean(rr2)
rr2<-as.list(rr2)
rr1[[1]]
S2 <- list()
S2[[1]] <- 0
for(i in 1:5){
  S2[[i+1]]<-S2[[i]] + (rr2[[i]]-mn.rr2)  
}
rr[[34]]
s_2=S2 %>% unlist()
max.S2<-max(s_2)
S.diff<-max.S-min.S
c=which(s_2==max.S2)
# c=3; so,4+2=6
# so we get our second change point at time=6 


set.seed(100)

### for generating 10 datasets.
nsim=1000
n=500
 beta <- c(0.50,1.15)
  Z <- matrix(rnorm(n*2), ncol = 2) #covariates
 #for inverse transform simulation
  #time simulation based on inverse CDF transform method - refer to paper
  alpha <- 0.75
  tau <- c(2.50,6)
  lambda <-c(0.05,0.35,0.55)

#### simulated data generation
n=500
dat<- list()
censor=20
for(i in 1:nsim){
  dat[[i]] <- tibble(
    u = runif(n),
    t = as.numeric(((-log(u) - (lambda[1]*exp(Z%*%beta)*tau[1]^alpha) -
                            (lambda[2]*exp(Z%*%beta)*(tau[2]^alpha - tau[1]^alpha)) +
             (lambda[3]*exp(Z%*%beta)*tau[2]^alpha))/(lambda[3]*exp(Z%*%beta)))^(1/alpha ) ) ,
  
    cen_time= runif(n,0,quantile(t,0.95)),
    t1=pmin(t,cen_time),
  status = ifelse(t <= cen_time, 1, 0),
  delta = ifelse(t <= cen_time, 1, 0),
  t2=ceiling(t1)
  
 )
}

#### finding the censoring percentages for 100 simulated datasets.
prop.table(table(dat[[100]]$status))
pr.tb<-list()
for(i in 1:nsim){
  pr.tb[[i]]=prop.table(table(dat[[i]]$status))
}

pr.tb[[1]]
#########################################
### performiong NSP in two simulated data
#########################################
library(survival)
s=list()
s2=list()
s3=list()
Dt<-list()
At<-list()
nsp_cal<-list()
interval<-list()
int_u<-list()
int_l<-list()
i_u=NULL
i_l=NULL
i_u1=NULL
i_l1=NULL
conf=NULL
for(i in 1:nsim){
  s[[i]]=survfit(Surv(t2,delta)~1,data=dat[[i]])
  s2[[i]]=summary(s[[i]])
  s3[[i]]=s2[[i]]$n.event
  Dt[[i]]<-s3[[i]]
  At[[i]]<-2*sqrt( Dt[[i]] + 3/8)
  nsp_cal[[i]]<-nsp_poly(At[[i]],deg = 1)
  interval[[i]]<-nsp_cal[[i]]$intervals
  int_u[[i]]<-interval[[i]]$ends
  int_l[[i]]<-interval[[i]]$starts
  i_u<-int_u %>% unlist()
  i_l<-int_l%>%unlist()
  i_u<-i_u[i_u<2]
  i_l<-i_l[i_l>6]
  conf <- (100*(1-(length(i_u)/nsim+length(i_l)/nsim)))

}
conf
int_l
# So, when we are working with 100 datasets,92% time the simulated confidence interval 
#fall on true confidence interval.

At[[10]]

############################################################
############ Estimating the change points ##################
############################################################


#######################################
#######################################
#######################################
rr1<-list()
mn.rr1<-NULL
S1<-list()
max.S1<-list()
c1<-list()
rr2<-list()
mn.rr2<-NULL
S2<-list()
max.S2<-list()
c2<-list()
 for(i in 1:nsim){
  rr1[[i]]=At[[i]][1:5]
  mn.rr1[i]=mean(rr1[[i]])
  S1[i][1]=0
  for(j in 1:5){
    S1[[i]][j+1] = S1[[i]][j] + (rr1[[i]][j] - mn.rr1[i])
  }
  max.S1[[i]]<-max(S1[[i]])
  c1[[i]]=which(S1[[i]]==max.S1[[i]])
  cp1<-mean(unlist(c1))
  ##### 2nd cp #####
  rr2[[i]]=At[[i]][cp1:9]
  mn.rr2[i]<-mean(rr2[[i]])
  S2[i][1]=0
  for(k in 1:length(rr2[[i]])){
    S2[[i]][k+1] = S2[[i]][k] + (rr2[[i]][k] - mn.rr2[k])
  }
  max.S2[[i]]<-max(S2[[i]])
  c2[[i]]=which(S2[[i]]==max.S2[[i]])
  cp2<-mean(unlist(c2)) 
  
  # cp2 = 3.7 ~ 4
  # THAT MEANS IF 1ST CP IS 2.5 THEN 
  # 2ND CP IS cP:9 = 6
  cp22<-6
  }


library(SurvRegCensCov)
dat1=list()
dat_cp1<-list()
WR1<-list()
lambda1<-list()
alpha1<-list()
SR<-list()
beta1<-list()
beta2<-list()

############
dat_cp2<-list()
WR2<-list()
lambda2<-list()
alpha2<-list()
SR2<-list()
beta.1<-list()
beta.2<-list()

############
dat_cp3<-list()
WR3<-list()
lambda3<-list()
alpha3<-list()
SR3<-list()
beta_1<-list()
beta_2<-list()
max_time<-list()
sim.tab<-tibble()
for(i in 1:nsim){
  dat1[[i]]<-dat[[i]]%>%
    arrange(t1)
  max_time[[i]]<-max(dat1[[i]]$t1)
  ###########
  # for 1st cp
  ###########
  
  dat_cp1[[i]]<-dat1[[i]]%>%
    filter(t1<=2.5)
  WR1[[i]]<-WeibullReg(Surv(dat_cp1[[i]]$t1, dat_cp1[[i]]$delta) ~ Z[1:nrow(dat_cp1[[i]]),], data=dat_cp1[[i]])
  lambda1[[i]]<- WR1[[i]]$coef[1]
  lm1<-mean(unlist(lambda1))
  alpha1[[i]]<-WR1[[i]]$coef[2]
  al1<-mean(unlist(alpha1))
  SR[[i]]<-WR1[[i]]$summary
  beta1[[i]]<-SR[[i]]$coefficients[2]
  beta2[[i]]<-SR[[i]]$coefficients[3]
  b1<-mean(unlist(beta1))
  b2<-mean(unlist(beta2))
  
  #############
  # for 2nd cp
  ###########
  
  dat_cp2[[i]]<-dat1[[i]]%>%
    filter(t1>2.5)%>%
    filter(t1<=6)

  WR2[[i]]<-WeibullReg(Surv(dat_cp2[[i]]$t1, dat_cp2[[i]]$delta) ~ Z[1:nrow(dat_cp2[[i]]),], data=dat_cp2[[i]])
    lambda2[[i]]<- WR2[[i]]$coef[1]

  lm2<-mean(unlist(lambda2))

  alpha2[[i]]<-WR2[[i]]$coef[2]
  al2<-mean(unlist(alpha2))

  SR2[[i]]<-WR2[[i]]$summary

  beta.1[[i]]<-SR2[[i]]$coefficients[2]
  beta.2[[i]]<-SR2[[i]]$coefficients[3]
  b.1<-mean(unlist(beta.1))
  b.2<-mean(unlist(beta.2))
  
  
  ############################
  #### 3rd segment ##########
  ##########################
  
  dat_cp3[[i]]<-dat1[[i]]%>%
    filter(t1>6)%>%
    filter(t1<=max_time[[i]])

  WR3[[i]]<-WeibullReg(Surv(dat_cp3[[i]]$t1, dat_cp3[[i]]$delta) ~ Z[1:nrow(dat_cp3[[i]]),], data=dat_cp3[[i]])

  lambda3[[i]]<- WR3[[i]]$coef[1]
  
  lm3<-mean(unlist(lambda3))

  alpha3[[i]]<-WR3[[i]]$coef[2]
  al3<-mean(unlist(alpha3))

  SR3[[i]]<-WR3[[i]]$summary

  beta_1[[i]]<-SR3[[i]]$coefficients[2]
  beta_2[[i]]<-SR3[[i]]$coefficients[3]
  b_1<-mean(unlist(beta_1))
  b_2<-mean(unlist(beta_2))
  
  sim.tab<-tibble(lambda = c(lm1,lm2,lm3),
                         alpha = c(al1,al2,al3),
                         beta1=c(b1,b.1,b_1),
                         beta2 = c(b2,b.2,b_2))
  
  
}
