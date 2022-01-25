#Load df4 by running 'load_dataset.R'
#Please note that running the Stan models can be slow!

library(rstan)
library(dplyr)
#rstan_options(auto_write=TRUE)
options(mc.cores = 4) #Remove this if you don't have a multicore machine
library(rethinking)
library(ggplot2)
library(gridExtra)
#df4 loaded elsewhere

#First do all data: (basic & severe only)
df4$logvalue <- log10(df4$value)
str(df4)

#For stan, StudyNum needs to be continuous integers
lsst <- unique(df4$StudyNum) #NOTE: We've lost Study P
df4$StudyNum2 <- NA
for(i in 1:length(lsst)){
  df4$StudyNum2[df4$StudyNum == lsst[i]] <- i
}

#Ditto for PatID
lsst <- unique(df4$PatientID)
df4$PatID <- NA
for(i in 1:length(lsst)){
  df4$PatID[df4$PatientID == lsst[i]] <- i
}

df9 <- select(df4,'Day','Estimated','logvalue','StudyNum2','PatID','SevMax3','LOD')
#make mild default
#df9$ModYN <- ifelse(df9$SevMax3 == 'Moderate',1,0)
#df9$SevYN <- ifelse(df9$SevMax3 == 'Severe',1,0)
#Not mild
df9$NotMild <- ifelse(df9$SevMax3 != 'Mild', 1, 0)
df9$Severity2 <- NA
df9[df9$NotMild==0,]$Severity2 <- 'Mild'
df9[df9$NotMild==1,]$Severity2 <- 'Moderate or Severe'

get_waic <- function(stanfit){
  # compute WAIC from the returned object from Stan
  # the log likelihood must be named as 'log_lik'
  waic <- function(log_lik) {
    Tn <- - mean(log(colMeans(exp(log_lik))))
    fV <- mean(colMeans(log_lik^2) - colMeans(log_lik)^2)
    waic <- Tn + fV
    waic
  }
  
  stanfit %>% rstan::extract() %>% .$log_lik %>% waic()
}

#data for stan. 
dfst <- list(M =length(unique(df9$PatID)),
             S =length(unique(df9$StudyNum2)),
             L = dim(df9)[1],
             LOD = df9$LOD,
             #SevYN = df9$SevYN,
             #ModYN = df9$ModYN,
             NotMild = df9$NotMild, 
             logvalue = df9$logvalue,
             Day = df9$Day,
             StudyNum = df9$StudyNum2,
             PatID = df9$PatID
)
fit_mA1 <- stan(file = 'report_regression_A1.stan', data = dfst,
                chains = 4, cores = 4, iter = 6000, warmup = 1000, control = list(adapt_delta = 0.95)) 
fit_mA2 <- stan(file = 'report_regression_A2nm.stan', data = dfst,
                chains = 4, cores = 4, iter = 6000, warmup = 1000, control = list(adapt_delta = 0.95)) 

#Output for Suppl. Table 3
smpl <- extract(fit_mA1)
mean(smpl$b0) #av. slope
mean(smpl$a0) #av. intercept
mean(smpl$sigma) #av. sigma
mean(smpl$sigma_study[,1]) #av. sigma
mean(smpl$sigma_study[,2]) #av. sigma
mean(smpl$sigma_PatID[,1]) #av. sigma
mean(smpl$sigma_PatID[,2]) #av. sigma
mean(smpl$sigma) #av. sigma

#Use fit_mA1 to generate the multi-panel regression figure.

S <- length(unique(df9$StudyNum))
dfe <- data.frame('StudyNum' = integer(), 'Mean'=numeric(), 'lwr' = numeric(),
                  'uppr'=numeric(), 'Day' = integer(), 'LOD' = numeric())
W <- 1000 #how many samples to use
simm <- matrix(0, nrow=W, ncol= 16)
simmA <- matrix(0, nrow=W, ncol= 16)
for(i in 1:S){
  for(j in 1:W){
    ls <- seq(0,15,1)
    ls2 <- (smpl$a0[j] + smpl$a1[j,i]) + (smpl$b0[j] + smpl$b1[j,i])*ls
    ls3 <- (smpl$a0[j]) + (smpl$b0[j])*ls
    simm[j,] <- ls2
    simmA[j,] <- ls3
  }
  muX <- apply(simm, 2, mean)
  mu.PIX <- apply(simm, 2, PI, prob=0.95)
  aux <- data.frame('StudyNum' = rep(i,16), 'Mean'=muX, 'lwr' = mu.PIX[1,],
                    'uppr' = mu.PIX[2,], 'Day' = seq(0,15,1), 
                    'LOD' = rep(df9[df9$StudyNum == i,]$LOD[1],16))
  dfe <- rbind(dfe,aux)
}
muXA <- apply(simmA, 2, mean)
mu.PIXA <- apply(simmA, 2, PI, prob=0.95)
dfaa <- data.frame('Day' = seq(0,15,1), 'Mean'=muXA, 'lwr' = mu.PIXA[1,],
                   'uppr' = mu.PIXA[2,])

#UPDATE STUDY NUMBERS TO MATCH TABLE 1. In dfe & df9. 1-12 the same. Then skip 13.
df9$StudyNum3 <- df9$StudyNum2
df9[df9$StudyNum2==16,]$StudyNum3 <- 17
df9[df9$StudyNum2==15,]$StudyNum3 <- 16
df9[df9$StudyNum2==14,]$StudyNum3 <- 15
df9[df9$StudyNum2==13,]$StudyNum3 <- 14

dfe$StudyNum3 <- dfe$StudyNum
dfe[dfe$StudyNum==16,]$StudyNum3 <- 17
dfe[dfe$StudyNum==15,]$StudyNum3 <- 16
dfe[dfe$StudyNum==14,]$StudyNum3 <- 15
dfe[dfe$StudyNum==13,]$StudyNum3 <- 14

#Figure2
ggplot() + geom_line(data = dfe, aes(x=Day, y=Mean), color = 'black') + facet_wrap(~StudyNum3) +  
  ylab("log10 (Viral copies / ml)") + xlab("Days after symptom onset") +
  geom_point(data = df9, aes(x=Day, y=logvalue, color = factor(StudyNum3), shape = factor(Estimated)),
             alpha = .4) + themeJDC +
  geom_ribbon(data = dfaa, aes(x=Day, ymin = lwr, ymax = uppr), alpha=.1, fill = 'black') + 
  geom_ribbon(data = dfe, aes(x=Day, ymin = lwr, ymax = uppr, fill = factor(StudyNum3)), alpha=.2) + 
  theme(strip.background = element_rect(fill = 'white'), legend.position = 'none') +
  geom_line(data=dfaa, aes(x=Day,y=Mean),linetype='dashed',alpha=.5) + 
  geom_line(data=dfe, aes(x=Day, y=log10(LOD)),color = 'grey',alpha=.99) +
  coord_cartesian(ylim=c(-0.2, NA)) + scale_shape_manual(values = c(20,18))


########### Make the severity plot (fit_mA2) #################
smpl <- extract(fit_mA2)
dfe <- data.frame('StudyNum' = integer(), 'Mean'=numeric(), 'lwr' = numeric(),
                  'uppr'=numeric(), 'Day' = integer(), 'LOD' = numeric())
W <- 1000 #how many samples to use
simm <- matrix(0, nrow=W, ncol= 16) #mild
simmA <- matrix(0, nrow=W, ncol= 16) #not mild
#simmB <- matrix(0, nrow=W, ncol= 16) #severe
for(j in 1:W){
  ls <- seq(0,15,1)
  ls1 <- (smpl$a0[j]) + (smpl$b0[j])*ls
  ls2 <- (smpl$a0[j] + smpl$aSnM[j]) + (smpl$b0[j] + smpl$bSnM[j])*ls
  #ls3 <- (smpl$a0[j] + smpl$aSS[j]) + (smpl$b0[j] + smpl$bSS[j])*ls
  simm[j,] <- ls1
  simmA[j,] <- ls2
  #simmB[j,] <- ls3
}
mu <- apply(simm, 2, mean)
mu.PI <- apply(simm, 2, PI, prob=0.95)
aux <- data.frame('Mean'=mu, 'lwr' = mu.PI[1,],
                  'uppr' = mu.PI[2,], 'Day' = seq(0,15,1), 'Severity' = rep('Mild'))
dfe <- rbind(dfe,aux)
muA <- apply(simmA, 2, mean)
mu.PA <- apply(simmA, 2, PI, prob=0.95)
aux <- data.frame('Mean'=muA, 'lwr' = mu.PA[1,],
                  'uppr' = mu.PA[2,], 'Day' = seq(0,15,1), 'Severity' = rep('Moderate or Severe'))
dfe <- rbind(dfe,aux)
#muB <- apply(simmB, 2, mean)
#mu.PB <- apply(simmB, 2, PI, prob=0.95)
#aux <- data.frame('Mean'=muB, 'lwr' = mu.PB[1,],
#                  'uppr' = mu.PB[2,], 'Day' = seq(0,15,1), 'Severity' = rep('Severe'))
#dfe <- rbind(dfe,aux)
pw1 <- ggplot() + theme(legend.position = c(0.85,0.85)) +  
  geom_point(data=df9, aes(x=Day,y=logvalue,color = Severity2), alpha=.75) +
  geom_ribbon(data=dfe, aes(x=Day, ymin = lwr, ymax = uppr, fill = Severity, group=Severity),alpha=.4) +
  geom_line(data=dfe, aes(x=Day, y=Mean, color = Severity, group=Severity)) +
  #scale_fill_manual(values = c('red','purple','blue'), guide = 'none') + 
  scale_fill_manual(values = c("#00A6A6","#F08700"), labels = c('Mild','Moderate\nor Severe')) + 
  scale_color_manual(values = c('#00A6A6','#F08700'), guide = 'none') + themeJDC + 
  #scale_color_manual(values = c('red','purple','blue'), guide = 'none') + themeJDC + 
  #guides(fill = guide_legend(NULL), color = guide_legend(NULL))   + 
  ylab("log10 (Viral copies / ml)") + xlab("Days after symptom onset")
pw1


########### Studies with all variables present ################
table(df4[is.na(df4$Sex)==T,]$Study)
table(df4[is.na(df4$Age)==T,]$Study)
table(df4[is.na(df4$SevMax)==T,]$Study)

df4 <- df4[is.na(df4$Sex)==F,]

#For ulam, r.e. should be continuous integers?
lsst <- unique(df4$StudyNum)
df4$StudyNum2 <- NA
for(i in 1:length(lsst)){
  df4$StudyNum2[df4$StudyNum == lsst[i]] <- i
}

#Ditto for PatID
lsst <- unique(df4$PatientID)
df4$PatID <- NA
for(i in 1:length(lsst)){
  df4$PatID[df4$PatientID == lsst[i]] <- i
}

df10 <- select(df4,'Day','Estimated','logvalue','StudyNum2','PatID','SevMax3','LOD','Age','Sex')
#make mild default
#df10$ModYN <- ifelse(df10$SevMax3 == 'Moderate',1,0)
#df10$SevYN <- ifelse(df10$SevMax3 == 'Severe',1,0)
df10$NotMild <- ifelse(df10$SevMax3 != 'Mild', 1, 0)
df10$Severity2 <- NA
df10[df10$NotMild==0,]$Severity2 <- 'Mild'
df10[df10$NotMild==1,]$Severity2 <- 'NotMild'

#make <40 default
df10$AGmid <- ifelse(df10$Age > 40 & df10$Age < 60,1,0)
df10$AGold <- ifelse(df10$Age > 59,1,0)

df10$Male <- ifelse(df10$Sex=='M',1,0)

# new data for stan. 
dfst2 <- list(M =length(unique(df10$PatID)),
              S =length(unique(df10$StudyNum)),
              L = dim(df10)[1],
              LOD = df10$LOD,
              #SevYN = df10$SevYN,
              #ModYN = df10$ModYN,
              NotMild = df10$NotMild, 
              AGmid = df10$AGmid,
              AGold = df10$AGold,
              Male = df10$Male,
              logvalue = df10$logvalue,
              Day = df10$Day,
              StudyNum = df10$StudyNum2,
              PatID = df10$PatID
)
# new data for stan. (same but for the 'A' models) 
dfst2a <- list(M =length(unique(df10$PatID)),
               S =length(unique(df10$StudyNum2)),
               L = dim(df10)[1],
               LOD = df10$LOD,
               NotMild = df10$NotMild, 
               #SevYN = df10$SevYN,
               #ModYN = df10$ModYN,
               #AGmid = df10$AGmid,
               #AGold = df10$AGold,
               #Male = df10$Male,
               logvalue = df10$logvalue,
               Day = df10$Day,
               StudyNum = df10$StudyNum2,
               PatID = df10$PatID
)
# new data for stan, now age is continuous
# dfst2y <- list(M =length(unique(df10$PatID)),
#                S =length(unique(df10$StudyNum2)),
#                L = dim(df10)[1],
#                LOD = df10$LOD,
#                SevYN = df10$SevYN,
#                ModYN = df10$ModYN,
#                AgeLG = (df10$Age - mean(df10$Age))/sd(df10$Age),#log(df10$Age), #log and substract mean??
#                Male = df10$Male,
#                logvalue = df10$logvalue,
#                Day = df10$Day,
#                StudyNum = df10$StudyNum2,
#                PatID = df10$PatID
# )

fit_mB1 <- stan(file = 'report_regression_A1.stan', data = dfst2a,
                chains = 3, cores = 3, iter = 6000, warmup = 1000) 
print(fit_mB1, pars = c("a0", "b0", "sigma"))
fit_mB2 <- stan(file = 'report_regression_A2nm.stan', data = dfst2a,
                chains = 4, cores = 4, iter = 8000, warmup = 1000, control = list(adapt_delta = 0.95)) 
print(fit_mB2, pars = c("a0", "b0", "sigma"))
# fit_mB3 <- stan(file = 'report_regression_B2.stan', data = dfst2,
#                 chains = 3, cores = 3, iter = 6000, warmup = 1000) 
# print(fit_mB3, pars = c("a0", "b0", "sigma"))
fit_mB3x <- stan(file = 'report_regression_B2x.stan', data = dfst2,
                 chains = 4, cores = 4, iter = 6000, warmup = 1000, control = list(adapt_delta = 0.96)) 
print(fit_mB3x, pars = c("a0", "b0","aA2","bA2", "sigma"))
#fit_mB3y <- stan(file = 'report_regression_B2y.stan', data = dfst2y,
#                 chains = 4, cores = 4, iter = 6000, warmup = 1000, control = list(adapt_delta = 0.96)) 
#print(fit_mB3y, pars = c("a0", "b0","aa2","ba2", "sigma"))
#compare(fit_mB3x,fit_mB3y)
fit_mB4 <- stan(file = 'report_regression_B3.stan', data = dfst2,
                chains = 4, cores = 4, iter = 6000, warmup = 1000, 
                control = list(adapt_delta = 0.95, max_treedepth = 15)) 
print(fit_mB4, pars = c("a0", "b0", "sigma"))
# fit_mB5 <- stan(file = 'report_regression_B4.stan', data = dfst2,
#                 chains = 3, cores = 3, iter = 6000, warmup = 1000) 
# print(fit_mB5, pars = c("a0", "b0", "sigma"))
fit_mB5x <- stan(file = 'report_regression_B4xnm.stan', data = dfst2,
                 chains = 4, cores = 4, iter = 6000, warmup = 1000, 
                 control = list(adapt_delta = 0.95, max_treedepth = 15)) 
print(fit_mB5x, pars = c("a0", "b0", "sigma"))
#compare(fit_mB5,fit_mB5x)
fit_mB6 <- stan(file = 'report_regression_B5.stan', data = dfst2,
                chains = 4, cores = 4, iter = 6000, warmup = 1000, control = list(adapt_delta = 0.95)) 
print(fit_mB6, pars = c("a0", "b0", "sigma"))
# fit_mB7 <- stan(file = 'report_regression_B6.stan', data = dfst2,
#                 chains = 3, cores = 3, iter = 6000, warmup = 1000) 
# print(fit_mB7, pars = c("a0", "b0","aGM","bGM",'aA1','aA2','bA1','bA2', "sigma"))
fit_mB7x <- stan(file = 'report_regression_B6x.stan', data = dfst2,
                 chains = 4, cores = 4, iter = 6000, warmup = 1000, control = list(adapt_delta = 0.95)) 
print(fit_mB7x, pars = c("a0", "b0","aGM","bGM",'aA2','bA2', "sigma"))
#compare(fit_mB7,fit_mB7x)
# fit_mB8 <- stan(file = 'report_regression_B7.stan', data = dfst2,
#                 chains = 3, cores = 3, iter = 6000, warmup = 1000) 
# print(fit_mB8, pars = c("a0", "b0", "sigma"))
fit_mB8x <- stan(file = 'report_regression_B7x.stan', data = dfst2,
                 chains = 4, cores = 4, iter = 14000, warmup = 1000, 
                 control = list(adapt_delta = 0.98,max_treedepth = 17)) 
print(fit_mB8x, pars = c("a0", "b0","aGM","bGM",'aA2','bA2', 'aSnM', 'bSnM' , "sigma"))
#pairs(fit_mB8x, pars = c("a0", "b0","aGM","bGM",'aA2','bA2', 'aSM', 'aSS', 'bSM' ,'bSS', "sigma"))
#compare(fit_mB8,fit_mB8x)

compare(fit_mB1, fit_mB2, fit_mB3x, fit_mB4, fit_mB5x, fit_mB6,
        fit_mB7x, fit_mB8x)

########### Make the (ALTERNATIVE) age plot (fit_mB3x) #################
df10$AG <- NA
df10$AG[df10$Age < 60] <- 0
#df10$AG[df10$Age > 40 & df10$Age <60] <- 1
df10$AG[df10$Age > 59] <- 1

smpl <- extract(fit_mB3x)
dfe <- data.frame('Age' = integer(), 'Mean'=numeric(), 'lwr' = numeric(),
                  'uppr'=numeric(), 'Day' = integer())
W <- 1000 #how many samples to use
simm <- matrix(0, nrow=W, ncol= 16) #<60
simmA <- matrix(0, nrow=W, ncol= 16) #60+
for(j in 1:W){
  ls <- seq(0,15,1)
  ls1 <- (smpl$a0[j]) + (smpl$b0[j])*ls #u60s
  ls2 <- (smpl$a0[j] + smpl$aA2[j]) + (smpl$b0[j] + smpl$bA2[j])*ls #60+
  simm[j,] <- ls1
  simmA[j,] <- ls2
  
}
mu <- apply(simm, 2, mean)
mu.PI <- apply(simm, 2, PI, prob=0.95)
aux <- data.frame('Mean'=mu, 'lwr' = mu.PI[1,],
                  'uppr' = mu.PI[2,], 'Day' = seq(0,15,1), 'AG' = rep("0",16))
dfe <- rbind(dfe,aux)
muA <- apply(simmA, 2, mean)
mu.PA <- apply(simmA, 2, PI, prob=0.95)
aux <- data.frame('Mean'=muA, 'lwr' = mu.PA[1,],
                  'uppr' = mu.PA[2,], 'Day' = seq(0,15,1), 'AG' = rep("1",16))
dfe <- rbind(dfe,aux)
pw2 <- ggplot() + theme(legend.position = c(0.85,0.85)) +
  geom_point(data=df10, aes(x=Day,y=logvalue,color = factor(AG)), alpha=.4) +
  geom_ribbon(data=dfe, aes(x=Day, ymin = lwr, ymax = uppr, fill = factor(AG),
                            group=factor(AG)),alpha=.4) +
  geom_line(data=dfe, aes(x=Day, y=Mean, color = factor(AG), group=factor(AG))) +
  #scale_fill_manual(values = c('red','purple','blue'), guide = 'none') + 
  scale_fill_manual(values = c('orange','darkblue'), labels = c('Under 60','Over 60')) + 
  scale_color_manual(values = c('orange','darkblue'), guide = 'none') + themeJDC + 
  #scale_color_manual(values = c('red','purple','blue'), guide = 'none') + themeJDC + 
  guides(fill = guide_legend(title='Age (years)'))   + 
  ylab("log10 (Viral copies / ml)") + xlab("Days after symptom onset")
pw2

########### Make the Sex plot (fit_mB4) #################
smpl <- extract(fit_mB4)
dfe <- data.frame('Sex' = integer(), 'Mean'=numeric(), 'lwr' = numeric(),
                  'uppr'=numeric(), 'Day' = integer())
W <- 1000 #how many samples to use
simm <- matrix(0, nrow=W, ncol= 16) #F
simmA <- matrix(0, nrow=W, ncol= 16) #M
for(j in 1:W){
  ls <- seq(0,15,1)
  ls1 <- (smpl$a0[j]) + (smpl$b0[j])*ls #female
  ls2 <- (smpl$a0[j] + smpl$aGM[j]) + (smpl$b0[j] + smpl$bGM[j])*ls
  simm[j,] <- ls1
  simmA[j,] <- ls2
}
mu <- apply(simm, 2, mean)
mu.PI <- apply(simm, 2, PI, prob=0.95)
aux <- data.frame('Mean'=mu, 'lwr' = mu.PI[1,],
                  'uppr' = mu.PI[2,], 'Day' = seq(0,15,1), 'Male' = rep(0,16))
dfe <- rbind(dfe,aux)
muA <- apply(simmA, 2, mean)
mu.PA <- apply(simmA, 2, PI, prob=0.95)
aux <- data.frame('Mean'=muA, 'lwr' = mu.PA[1,],
                  'uppr' = mu.PA[2,], 'Day' = seq(0,15,1), 'Male' = rep(1,16))
dfe <- rbind(dfe,aux)
pw3 <- ggplot() + theme(legend.position = c(0.85,0.85)) + 
  geom_point(data=df10, aes(x=Day,y=logvalue,color = factor(Male)), alpha=.5) +
  geom_ribbon(data=dfe, aes(x=Day, ymin = lwr, ymax = uppr, fill = factor(Male),
                            group=Male),alpha=.4) +
  geom_line(data=dfe, aes(x=Day, y=Mean, color = factor(Male), group=factor(Male))) +
  #scale_fill_manual(values = c('red','purple','blue'), guide = 'none') + 
  scale_fill_manual(values = c("#00B4D8","#03045E"), labels = c('Female','Male')) + 
  scale_color_manual(values = c("#00B4D8","#03045E"), guide = 'none') + themeJDC + 
  #scale_color_manual(values = c('red','purple','blue'), guide = 'none') + themeJDC + 
  guides(fill = guide_legend(title=NULL))   + 
  ylab("log10 (Viral copies / ml)") + xlab("Days after symptom onset")
pw3

#Figure 3
grid.arrange(pw1,pw2,pw3,nrow=1)

#g2 <- gridExtra::arrangeGrob(pw1,pw2,pw3,nrow = 1)
#ggsave(file="Report/Report_figures/Three_panel_regression_atleast3.pdf", g2,height = 6.0, width = 14.0)