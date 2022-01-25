#df loaded
source("time_to_admission.R")
mycolors <- c("#FFC20A","#0C7BDC")

#Suppl. Figure 2
ggplot(df, aes(x=Day, y=log10(value),group = PatientID,color = factor(Estimated))) + 
  geom_line(alpha=.6) + geom_point(alpha=.25, size = .9) + xlab('Days after symptom onset') +
  facet_wrap(~StudyNum) + 
  ylab("log10 (Viral copies per ml)") + theme_bw() +
  theme(strip.background = element_rect(colour="black", fill="white")) +
  scale_color_manual(name='Original data',labels=c("Viral load","Cycle\nthreshold"),values = mycolors)

dfT <- med_and_quantiles(df)
dfTlg <- med_and_quantiles_lg10(df)
dfT0 <- med_and_quantiles(df[df$Estimated==0,])
dfT1 <- med_and_quantiles(df[df$Estimated==1,])
dfT0lg <- med_and_quantiles_lg10(df[df$Estimated==0,])
dfT1lg <- med_and_quantiles_lg10(df[df$Estimated==1,])

pl_lg <- ggplot() + theme_bw() + xlim(-2.3,23) + #scale_y_log10() +
  geom_ribbon(data = dfTlg, aes(x=Day, ymin = q1, ymax = q2), alpha=.4,fill = 'purple') + 
  xlab('Days after symptom onset') + 
  theme(legend.position = c(0.83,0.83), legend.background = element_rect(colour = 1)) + 
  ylab("log10 (Viral copies / mL)")  + scale_shape_manual(values = c(19,1), labels = c("Viral load","Cycle\nthreshold")) +
  geom_point(data=df, aes(x=Day, y=log10(value), shape = factor(Estimated)), color = 'black', alpha=.6) +
  geom_line(data = dfTlg, aes(x=Day, y = md),color = 'purple', size = 1) + 
  guides(shape=guide_legend("Original data"), color = FALSE)

pl10lg <- ggplot() + theme_bw() + xlim(-2.3,23) +
  geom_ribbon(data = dfT1lg, aes(x=Day, ymin = q1, ymax = q2), alpha=.4,fill = mycolors[2]) + 
  geom_ribbon(data = dfT0lg, aes(x=Day, ymin = q1, ymax = q2), alpha=.48,fill = mycolors[1]) +
  geom_point(data=df, aes(x=Day, y=log10(value), color = factor(Estimated)), shape=16, alpha=.63) +
  geom_line(data = dfT1lg, aes(x=Day, y = md),color = mycolors[2], size=1) +
  geom_line(data = dfT0lg, aes(x=Day, y = md),color = mycolors[1], size=1) + 
  xlab('Days after symptom onset') + ylab("log10 (Viral copies / mL") + 
  theme(legend.position = c(0.83,0.81),legend.background = element_rect(colour = 1)) + 
  scale_color_manual(name='Original data',labels=c("Viral load","Cycle\nthreshold"),values = mycolors)
pl10lg

ff <- seq(-2,23,1)

dff <- data.frame('dy'=ff,'cou'=table(df$Day)[1:26])
plb <-ggplot(dff, aes(x=dy,y=cou.Freq)) + geom_bar(stat = 'identity',fill='purple',alpha=.6) + 
  theme_bw() + xlab('Days after symptom onset') + ylab('No. of data points')
plb
dff2 <- data.frame('couVL'=table(df[df$Estimated==0,]$Day)[1:24])
dff2$couVL.Var1 <- as.numeric(as.character(dff2$couVL.Var1))
dff3 <- data.frame('couCt'=table(df[df$Estimated==1,]$Day)[1:26])
dff3$couCt.Var1 <- as.numeric(as.character(dff3$couCt.Var1))
dff2 <- rbind(dff2, data.frame('couVL.Var1'=c(-2,-1), 'couVL.Freq'=c(0,0)))
dff4 <- dff2[order(dff2$couVL.Var1),]
dff5 <- cbind(dff4,dff3)
pl10b <- ggplot(dff5) + theme_bw() + 
  geom_bar(aes(x=couCt.Var1,y=couVL.Freq+couCt.Freq),stat='identity',fill = mycolors[2],alpha=.6) +
  geom_bar(aes(x=couVL.Var1,y=couVL.Freq),stat='identity',fill = mycolors[1],alpha=.8) + 
  theme_bw() + xlab('Days after symptom onset') + ylab('No. of data points')
pl10b

#Figure 1
cowplot::plot_grid(pl_lg,pl10lg,plb,pl10b, nrow=2, labels = c('a','b','c','d'),
                   rel_heights = c(1,0.5,1,0.5))
