library(readxl)
library(ggplot2)

#The two datasets we'll use for the analysis are:
# (i) df (everything); 
# (ii) df4 (samples from no more than 15 days after symptom onset, and from patients with >2 samples)

df <- as.data.frame(read_excel('CombinedDataset.xlsx', sheet = 'Viral_Load'))

df2 <- df[df$Day<16 & df$Day > -3,] #because -9 is tough to deal with!

dim(df2)

#Remove patients with no +ve samples
lid <- unique(df2$PatientID)
l <- length(lid)
lid2 <- c()
for(i in 1:l){
  aux <- df2[df2$PatientID==lid[i],]
  if(max(aux$value) < 1.1){
    print(paste0('ooh ',lid[i], ' ', max(aux$value)))
    lid2 <- c(lid2, lid[i])
  }
}

df3 <- df2
for(j in 1:(length(lid2))){
  print(df3[df3$PatientID == lid2[j],]$value)
  df3 <- df3[df3$PatientID != lid2[j],]
}
dim(df3)

# Re-do: only use subjects with >2 samples?
lid <- unique(df2$PatientID)
l <- length(lid)
lid3 <- c()
innf <- rep(0,l)
for(i in 1:l){
  aux <- df2[df2$PatientID==lid[i],]
  innf[i] <- dim(aux)[1]
  if(dim(aux)[1] < 3){
    #print(paste0('ooh ',lid[i], ' ', max(aux$value)))
    lid3 <- c(lid3, lid[i])
  }
}

df4 <- df2
for(j in 1:length(lid3)){
  df4 <- df4[df4$PatientID != lid3[j],]
}

dim(df3)
dim(df4)
length(unique(df4$PatientID))

#Also, load theme for plots
themeJDC <- theme(panel.background = element_blank(),panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(colour = "black" ,fill = NA), 
                  legend.background = element_blank(), 
                  strip.background = element_rect(fill = 'white', color = 'black'),
                  legend.key = element_blank(), axis.text = element_text(size = 10.0),
                  axis.title = element_text(size=10.6), legend.text = element_text(size = 10.2),
                  strip.text.x = element_text(size = 10.2))
