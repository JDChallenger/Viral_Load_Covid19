library(ggplot2)
library(dde)
source('ode_model.R')

#Read in samples for the posterior distribution
#This dataframe contains subject-level, study-level, and global parameters
smpls <- readRDS('posterior_samples.rds')
W <- dim(smpls)[1]
ics <- readRDS('initial_conditions.rds')
#hist(log10(ics))
M <- 155 # no. of subjects
S <- 16 #no. of studies

#For now, just look at samples for the population-level immune responses (early & late respectively)
smpls2 <- smpls[,c((2*M + 2*S + 1), (2*M+2*S+4))]
#samples for the early immune response parameter
hist(smpls2[,1])
#samples for the late immune response parameter
hist(smpls2[,2])
mean(smpls2[,1])
mean(smpls2[,2])
sd(smpls2[,1])
sd(smpls2[,2])

#put these samples through the ODE model. 
#Then, calculate median trajectory, and central 90% interval
l <- dim(smpls2)[1]
tt <- 20
simm <- matrix(0, nrow = l, ncol = tt)
simmB <- matrix(0, nrow = l, ncol = tt)
yy <- rep(0,tt)
for(i in 1:l){
  
  F_kaQ <- smpls2[i,1]
  F_ImaxQ <- smpls2[i,2]
  y0 <- c(0, median(ics), 0,0,0) #initial conditions for ODE model at t=-5
  ttt <- seq(0, tt, length.out = tt+1)
  p <- c(F_kaQ, F_ImaxQ)
  y <- dde::dopri(y0, ttt, viral4dS, p, return_time = F)
  for(j in 1:tt){
    #incorporate sigma
    yy[j] <- rlnorm(1,meanlog = log(y[1+j,2]), sdlog = smpls[i,349])
  }
  simm[i,] <-y[2:(tt+1),2]
  simmB[i,] <-yy
  
}

#Calculate median trajectory (and the interval)
medd <- rep(0,tt)
lwer <- rep(0,tt)
uppr <- rep(0,tt)
lwerB <- rep(0,tt)
upprB <- rep(0,tt)
for(k in 1:tt){
  medd[k] <- quantile(simm[,k],0.5)
  lwer[k] <- quantile(simm[,k],0.05)
  uppr[k] <- quantile(simm[,k],0.95)
  lwerB[k] <- quantile(simmB[,k],0.05)
  upprB[k] <- quantile(simmB[,k],0.95)
}
dfa <- data.frame('median' = medd, 'lower'=lwer, 'upper'=uppr, 'lowerB'=lwerB, 'upperB'=upprB,
                  'day'=-5+seq(1,tt,1))
#plot
ggplot(dfa) + geom_line(aes(x=day,y=median)) + scale_y_log10() + theme_bw() + 
  geom_ribbon(aes(x=day,ymin=lower,ymax=upper),alpha=.45,fill = 'cyan') + 
  geom_ribbon(aes(x=day,ymin=lowerB,ymax=upperB),alpha=.2,fill = 'cyan') + 
  xlab('Days after symptom onset') + ylab('log10 (Viral copies / ml)') + 
  #Add data. Just showing data from patients  with >2 samples (as used to fit the model)
  geom_point(data=df4, aes(x=Day, y=value,color = factor(StudyNum))) + labs(color = 'Study')
