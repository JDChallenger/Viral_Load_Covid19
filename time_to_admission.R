#Extract first VL measure, to see if it varies with severity..
#Practice on study G, 11 patients
# str(dfG)
# table(dfG$PatientID)
# lst <- unique(dfG$PatientID)
# lst[1]
# l <- length(lst)

admission_day <- function(dfG,l){
  lst <- unique(dfG$PatientID)
  store_sev <- rep(0,l)
  store_time <- rep(0,l)
  store_est <- rep(0,l)
for(i in 1:l){
  aux <- dfG[dfG$PatientID==lst[i],]
  store_time[i] <- min(aux$Day)
  store_sev[i] <- unique(aux$SevMax3)[1]
  #store_est[i] <- unique(aux$Estimated)[1]
  if(length(unique(aux$SevMax3)) != 1 ){
    print(paste0('oopsy ',lst[i]))
  }
  if(is.na(store_time[i])==T ){
    print(paste0('oopsy ',lst[i]))
  }
}
aux2 <- data.frame('Severity' = store_sev, 'Day' = store_time)
return(aux2)
}
#admission_day(dfG,l)

first_VL <- function(dfG,l){
  lst <- unique(dfG$PatientID)
  store_sev <- rep(0,l)
  store_est <- rep(0,l)
  store_VL <- rep(0,l)
  store_AG <- rep(0,l)
  for(i in 1:l){
    aux <- dfG[dfG$PatientID==lst[i],]
    store_VL[i] <- aux$value[aux$Day == min(aux$Day)]
    store_sev[i] <- unique(aux$SevMax3)[1]
    store_AG[i] <- unique(aux$AG)[1]
    store_est[i] <- unique(aux$Estimated)[1]
    if(length(unique(aux$Severity)) > 1){
      print(paste0('oopsy ',lst[i]))
    }
  }
  aux2 <- data.frame('Severity' = store_sev, 'VL' = store_VL, 
                     'AgeGroup' = store_AG, 'Est' = store_est)
  return(aux2)
}

max_VL <- function(dfG,l){
  lst <- unique(dfG$PatientID)
  store_sev <- rep(0,l)
  store_VL <- rep(0,l)
  store_AG <- rep(0,l)
  store_est <- rep(0,l)
  for(i in 1:l){
    aux <- dfG[dfG$PatientID==lst[i],]
    store_VL[i] <- max(aux$value)
    store_sev[i] <- unique(aux$SevMax3)[1]
    store_AG[i] <- unique(aux$AG)[1]
    store_est[i] <- unique(aux$Estimated)[1]
    if(length(unique(aux$Severity)) > 1){
      print('oopsy')
    }
  }
  aux2 <- data.frame('Severity' = store_sev, 'VL' = store_VL,
                     'AgeGroup' = store_AG, 'Est' = store_est)
  return(aux2)
}


######
#Is the first VL the max one?
first_maxQ <- function(dfG,l){
  lst <- unique(dfG$PatientID)
  storeY <- rep(0,l)
  store_sev <- rep(0,l)
  store_AG <- rep(0,l)
  for(i in 1:l){
    aux <- dfG[dfG$PatientID==lst[i],]
    if(aux$value[which.min(aux$Day)]==max(aux$value)){
      storeY[i] <- 1
    }
    store_sev[i] <- unique(aux$Severity)[1]
    store_AG[i] <- unique(aux$AG)[1]
    if(length(unique(aux$Severity)) > 1){
      print('oopsy')
    }
    
  }
  aux2 <- data.frame('Severity' = store_sev, 'YN' = storeY, 'AgeGroup' = store_AG)
  return(aux2)
}
#tst <- first_maxQ(df,l)
#mean(tst$YN)


#########################
#Get medians and quantiles

med_and_quantiles <- function(dfG){
  md <- rep(0,30) #ALL
  q1 <- rep(0,30)
  q2 <- rep(0,30)
  mdM <- rep(0,30) #Mild
  q1M <- rep(0,30)
  q2M <- rep(0,30)
  mdMo <- rep(0,30) #Moderate
  q1Mo <- rep(0,30)
  q2Mo <- rep(0,30)
  mdS <- rep(0,30) #Severe
  q1S <- rep(0,30)
  q2S <- rep(0,30)
  for(i in 1:30){
    aux <- dfG[dfG$Day==(i-1),]
    q1[i] <- quantile(aux$value, 0.25,na.rm=T)
    q2[i] <- quantile(aux$value, 0.75,na.rm=T)
    md[i] <- median(aux$value,na.rm=T)
    aux <- dfG[dfG$Day==(i-1)& dfG$SevMax3=='Mild',]#dfG$Severity=='Mild',]
    q1M[i] <- quantile(aux$value, 0.25,na.rm=T)
    q2M[i] <- quantile(aux$value, 0.75,na.rm=T)
    mdM[i] <- median(aux$value,na.rm=T)
    aux <- dfG[dfG$Day==(i-1)& dfG$SevMax3=='Moderate',]#dfG$Severity=='Mild',]
    q1Mo[i] <- quantile(aux$value, 0.25,na.rm=T)
    q2Mo[i] <- quantile(aux$value, 0.75,na.rm=T)
    mdMo[i] <- median(aux$value,na.rm=T)
    aux <- dfG[dfG$Day==(i-1)& dfG$SevMax3=='Severe',]#dfG$Severity=='Severe',]
    q1S[i] <- quantile(aux$value, 0.25,na.rm=T)
    q2S[i] <- quantile(aux$value, 0.75,na.rm=T)
    mdS[i] <- median(aux$value,na.rm=T)
  }
  aux2 <- data.frame('Day' = seq(0,29,1),'md'=md,'q1'=q1,'q2'=q2,'mdM'=mdM,'q1M'=q1M,'q2M'=q2M,
                     'mdMo'=mdMo,'q1Mo'=q1Mo,'q2Mo'=q2Mo, 'mdS'=mdS,'q1S'=q1S,'q2S'=q2S)
  return(aux2)
}
#med_and_quantiles(df)

#re-do but log10 first???
med_and_quantiles_lg10 <- function(dfG){
  md <- rep(0,30) #ALL
  q1 <- rep(0,30)
  q2 <- rep(0,30)
  mdM <- rep(0,30) #Mild
  q1M <- rep(0,30)
  q2M <- rep(0,30)
  mdMo <- rep(0,30) #Moderate
  q1Mo <- rep(0,30)
  q2Mo <- rep(0,30)
  mdS <- rep(0,30) #Severe
  q1S <- rep(0,30)
  q2S <- rep(0,30)
  for(i in 1:30){
    aux <- dfG[dfG$Day==(i-1),]
    q1[i] <- quantile(log10(aux$value), 0.25,na.rm=T)
    q2[i] <- quantile(log10(aux$value), 0.75,na.rm=T)
    md[i] <- median(log10(aux$value),na.rm=T)
    aux <- dfG[dfG$Day==(i-1)& dfG$SevMax3=='Mild',]#dfG$Severity=='Mild',]
    q1M[i] <- quantile(log10(aux$value), 0.25,na.rm=T)
    q2M[i] <- quantile(log10(aux$value), 0.75,na.rm=T)
    mdM[i] <- median(log10(aux$value),na.rm=T)
    aux <- dfG[dfG$Day==(i-1)& dfG$SevMax3=='Moderate',]#dfG$Severity=='Mild',]
    q1Mo[i] <- quantile(log10(aux$value), 0.25,na.rm=T)
    q2Mo[i] <- quantile(log10(aux$value), 0.75,na.rm=T)
    mdMo[i] <- median(log10(aux$value),na.rm=T)
    aux <- dfG[dfG$Day==(i-1)& dfG$SevMax3=='Severe',]#dfG$Severity=='Severe',]
    q1S[i] <- quantile(log10(aux$value), 0.25,na.rm=T)
    q2S[i] <- quantile(log10(aux$value), 0.75,na.rm=T)
    mdS[i] <- median(log10(aux$value),na.rm=T)
  }
  aux2 <- data.frame('Day' = seq(0,29,1),'md'=md,'q1'=q1,'q2'=q2,'mdM'=mdM,'q1M'=q1M,'q2M'=q2M,
                     'mdMo'=mdMo,'q1Mo'=q1Mo,'q2Mo'=q2Mo, 'mdS'=mdS,'q1S'=q1S,'q2S'=q2S)
  return(aux2)
}



med_and_quantiles_age <- function(dfG, age = 40){
  md <- rep(0,20) #ALL
  q1 <- rep(0,20)
  q2 <- rep(0,20)
  mdU <- rep(0,20) #U40
  q1U <- rep(0,20)
  q2U <- rep(0,20)
  mdO <- rep(0,20) #Over 40
  q1O <- rep(0,20)
  q2O <- rep(0,20)
  for(i in 1:20){
    aux <- dfG[dfG$Day==(i-1),]
    q1[i] <- quantile(aux$value, 0.25,na.rm=T)
    q2[i] <- quantile(aux$value, 0.75,na.rm=T)
    md[i] <- median(aux$value,na.rm=T)
    aux <- dfG[dfG$Day==(i-1)&dfG$Age <= age,]
    q1U[i] <- quantile(aux$value, 0.25,na.rm=T)
    q2U[i] <- quantile(aux$value, 0.75,na.rm=T)
    mdU[i] <- median(aux$value,na.rm=T)
    aux <- dfG[dfG$Day==(i-1)&dfG$Age > age,]
    q1O[i] <- quantile(aux$value, 0.25,na.rm=T)
    q2O[i] <- quantile(aux$value, 0.75,na.rm=T)
    mdO[i] <- median(aux$value,na.rm=T)
  }
  aux2 <- data.frame('Day' = seq(0,19,1),'md'=md,'q1'=q1,'q2'=q2,'mdU'=mdU,'q1U'=q1U,'q2U'=q2U,
                     'mdO'=mdO,'q1O'=q1O,'q2O'=q2O)
  return(aux2)
}
#med_and_quantiles_age(df)

