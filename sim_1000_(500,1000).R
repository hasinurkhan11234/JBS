
set.seed(100)

### for generating 1000 datasets and sample size =500.
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

pr.tb<-list()
for(i in 1:nsim){
  pr.tb[[i]]=prop.table(table(dat[[i]]$status))
}


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




########################################################
#######################################################
###################################################

#### simulating 1000 times with sample size 1000
### for generating 1000 datasets and sample size =500.



nsim=1000
n=1000
beta <- c(0.50,1.15)
Z <- matrix(rnorm(n*2), ncol = 2) #covariates
#for inverse transform simulation
#time simulation based on inverse CDF transform method - refer to paper
alpha <- 0.75
tau <- c(2.50,6)
lambda <-c(0.05,0.35,0.55)

#### simulated data generation
n=1000
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

pr.tb<-list()
for(i in 1:nsim){
  pr.tb[[i]]=prop.table(table(dat[[i]]$status))
}


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
























