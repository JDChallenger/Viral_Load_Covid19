library(readxl)
library(ggplot2)
dsi <- read_excel('CombinedDataset.xlsx',skip = 4, sheet = 'Imm_Data')

dsi$Patient <- as.factor(dsi$ID)
levels(dsi$Patient) <- seq(1,12,1)
dsi2 <- as.data.frame(dsi)

#dsi2 <- dsi2[,c(1,2,5,12,13,14,15,16,17,18,19,20)]
#Change % inhibition to a decimal
dsi2$inhD <- (1/100)*dsi2$`% inhibition`

dsi3 <- dsi2[dsi2$`Days from symptom onset` < 40,]
dsi3 <- dsi3[,c(1,2,5,6,7,8)]
colnames(dsi3) <- c('ID','day','PC_inh','Tcell_tot','Patient','inhD')
dsi3$Patient <- as.integer(dsi3$Patient)

#logistic growth fn
xx3 <- function(t, a, b){
  rt <- 1/(1+exp(-(a+b*t)))
  return(rt)
}

#####################################
# Model neutralising antibody response using a logistic curve

dsi3$inhF <- dsi3$inhD
dsi3$inhF <- ifelse(dsi3$inhF > 0.005,dsi3$inhF, 0.005)
dsi3$inhG <- log(dsi3$inhF/(1-dsi3$inhF))

dfr <- data.frame('Patient'= integer(), 'a'=numeric(),'b'=numeric())
dfrG <- data.frame('Patient'= integer(), 'day'=integer(), 'inhH' = numeric())
#Loop over the 12 patients
for(j in 1:12){
  md <- lm(inhG ~ day, data = dsi3[dsi3$Patient==j & dsi3$day < 17,])
  #summary(md)['(Intercept)']
  dfr2 <- data.frame('Patient'=j,'a'=md$coefficients['(Intercept)'],'b'=md$coefficients['day'])
  dfr <- rbind(dfr,dfr2)
  
  v5 <- sapply(seq(0,30,1), xx3, b = dfr2$b, a = dfr2$a)
  df3x <- data.frame('Patient' = rep(j,31),'day'=seq(0,30,1), 'inhH' = v5)
  dfrG <- rbind(dfrG, df3x)
}

pAUC1 <- ggplot() + geom_line(data = dfrG, aes(x=day,y=inhH)) + facet_wrap(~Patient) + themeJDC +
  theme(strip.background = element_rect(fill = 'white', color = 'black')) +
  geom_point(data = dsi3, aes(x=day, y=inhF),color = 'purple') + xlim(c(0,30)) + 
  annotate('segment', x = 16.5, xend = 16.5, y = 0, yend = 1, linetype = 'dashed', alpha = .3) +
  ylab('Neutralising Antibody Response') + xlab('Days after symptom onset') +
  geom_area(data = dfrG[dfrG$day < 17,], aes(x=day,y=inhH), fill = 'purple', alpha = .2)
pAUC1

#Can now calculate AUCs. Need to choose time period to use
aucz <- rep(0,12)
for(j in 1:12){
  aux <- dfrG[dfrG$Patient==j,]
  sm <- 0
  for(i in 1:16){ #?
    sm <- sm + 0.5*(aux$inhH[i] + aux$inhH[i+1])*1 #dt = 1
  }
  aucz[j] <- sm
}
aucz

#####################################
# Model total T cell response using a logistic curve. (Need to provide a ceiling this time)

dsi3$Tcelltot2 <- (dsi3$Tcell_tot + 0.3) / (max(dsi3$Tcell_tot) + 1)
dsi3$logTcelltot3 <- log((dsi3$Tcelltot2)/(1-(dsi3$Tcelltot2)))

dfr <- data.frame('Patient'= integer(), 'a'=numeric(),'b'=numeric())
dfrT2 <- data.frame('Patient'= integer(), 'day'=integer(), 'lgTc' = numeric(),
                    'lgTc2' = numeric())
for(j in 1:12){
  md <- lm(logTcelltot3 ~ day, data = dsi3[dsi3$Patient==j & dsi3$day < 17,])
  if(j==6){
    print(summary(md))#['(Intercept)']
  }
  dfr2 <- data.frame('Patient'=j,'a'=md$coefficients['(Intercept)'],'b'=md$coefficients['day'])
  dfr <- rbind(dfr,dfr2)
  
  v5 <- sapply(seq(0,30,1), xx3, b = dfr2$b, a = dfr2$a)
  v52 <- dfr2$a + seq(0,30,1)*dfr2$b
  df3x <- data.frame('Patient' = rep(j,31),'day'=seq(0,30,1), 'lgTc' = v5, 'lgTc2' = v52 )
  dfrT2 <- rbind(dfrT2, df3x)
}

pAUC2 <- ggplot() + geom_line(data = dfrT2, aes(x=day,y=lgTc)) + facet_wrap(~Patient) + themeJDC + 
  geom_point(data = dsi3, aes(x=day, y=Tcelltot2), color = 'blue',shape = 7) + 
  xlim(c(0,30)) + theme(strip.background = element_rect(fill = 'white', color = 'black')) +
  annotate('segment', x = 1.7, xend = 1.7, y = 0, yend = 1, linetype = 'dashed', alpha = .3) +
  annotate('segment', x = 16.5, xend = 16.5, y = 0, yend = 1, linetype = 'dashed', alpha = .3) + 
  geom_area(data = dfrT2[dfrT2$day < 17 & dfrT2$day > 1 ,],
            aes(x = day, y = lgTc), fill = 'blue', alpha = .2) +
  ylab('Total T cell resp. (rescaled)') + 
  xlab('Days after symptom onset') 
pAUC2

#Can now calculate AUCs. Need to choose time period to use
aucT <- rep(0,12)
for(j in 1:12){
  aux <- dfrT2[dfrT2$Patient==j,]
  sm <- 0
  for(i in 3:16){ #? where to start??
    sm <- sm + 0.5*(aux$lgTc[i] + aux$lgTc[i+1])*1 #dt = 1
  }
  aucT[j] <- sm
}
aucT
