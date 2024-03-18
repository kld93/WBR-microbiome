# Manure Antibiotic absorbance UV-vis analysis
# July 19 2018
# Christine B. Georgakakos

# The purpose of this code is to clean the data, add columns, and prepare for analysis. 
# At the bottom of this code is breif plotting for data visualization

# Load Data
setwd('~/Dropbox/Soil and Water/antibiotics')

# This data is already zeroed according on day of measurment
#Data <- read.csv("Data/Trial1_0.45um_JerryDellFarm.csv") # Contains data for T1,T2, T3, CM, CA, CM2 filtered with 0.45um
Data <- read.csv("Data/Trial1_0.2um_JerryDellFarm.csv") # Contains data for T4, CM2, CA2 filtered with 0.2um

#Load SPAN DATA - recorded by computer
OneHour <- read.csv("Data/1hr_0.2um_T4,CM2,CA2_20180723.csv") #Survey Speed
ThreeHour <-read.csv("Data/3hr_0.2um_T4,CM2,CA2_20180723.csv") #Survey Speed
EightHour <-read.csv("Data/8hr_0.2um_T4,CM2,CA2_20180723.csv") #Medium Speed
TwentyFourHour <- read.csv("Data/24hr_0.2um_T4,CM2,CA2_20180724.csv") #Medium Speed
FortyEightHour <- read.csv("Data/48hr_0.2um_t4,CM2,CA2_20180725.csv") #Medium Speed

#Zero SPAN readings according to the MQ readings obtained in 3hr and 48hr
for (j in 1:nrow(ThreeHour)){
  ThreeHour$MQ_avg_Survey[j] <- mean(ThreeHour$MQ_3hr[j],ThreeHour$MQ_3hr1[j])
  for (k in 2:10){
    ThreeHour[,k][j] <- ThreeHour[,k][j] - ThreeHour$MQ_avg_Survey[j]
  }
  for (l in 2:5){
    OneHour[,l][j] <- OneHour[,l][j] - ThreeHour$MQ_avg_Survey[j]
  }
    
}
j<-1
k<-2
for (j in 1:nrow(EightHour)){
  for (k in 2:10){
    EightHour[,k][j] <- EightHour[,k][j] - FortyEightHour$MQ[j]
    TwentyFourHour[,k][j] <- TwentyFourHour[,k][j] - FortyEightHour$MQ[j]
    FortyEightHour[,k][j] <- FortyEightHour[,k][j] - FortyEightHour$MQ[j]
  }
  
}


# step 1: Add descriptive data to the sheet
# Write a for loop to assign manure and erythromycsin concentrations

for (i in 1:nrow(Data)){
  if (Data$Treatment[i] == "CM") {
    Data$Manure_mgL[i] <- 1000 # mg/L
    Data$Eryth_mgL[i] <- 0
  }else if (Data$Treatment[i] == "CA") {
    Data$Manure_mgL[i] <- 0
    Data$Eryth_mgL[i] <- 250 #mg/L
  }else if (Data$Treatment[i] == "T1"){
    Data$Manure_mgL[i] <- 1000
    Data$Eryth_mgL[i] <- 250 #mg/L
  }else if (Data$Treatment[i] == "T2") {
    Data$Manure_mgL[i] <- 1000
    Data$Eryth_mgL[i] <- 500 #mg/L
  }else if (Data$Treatment[i] =="CM2"){
    Data$Manure_mgL[i] <- 250 #mg/L
    Data$Eryth_mgL[i] <- 0
  }else if (Data$Treatment[i] == "T3"){
    Data$Manure_mgL[i] <- 250 #mg/L
    Data$Eryth_mgL[i] <- 250 #mg/L
  }else if (Data$Treatment[i] == "MQ"){
    Data$Manure_mgL[i] <- 0
    Data$Eryth_mgL[i] <- 0
  }else if (Data$Treatment[i] == "CM2"){
    Data$Manure_mgL[i] <- 250 #mg/L
    Data$Eryth_mgL[i] <- 0 
  }else if (Data$Treatment[i] == "T4"){
    Data$Manure_mgL[i] <- 250 #mg/L
    Data$Eryth_mgL[i] <- 625 #mg/L
  }else if (Data$Treatment[i] == "CA2"){
    Data$Manure_mgL[i] <- 0
    Data$Eryth_mgL[i] <- 625 #mg/L
  }
}

#***********End Data Calculations*************

#Plots for 0.2um filter pore
plot(Data$Hrs_Approx[Data$Treatment == "CM2"], Data$AbsorbAvg[Data$Treatment == "CM2"], pch = 3, xlab = "Hours", ylab = "Absorbance", ylim = c(0.04,0.18), main = "280 nm Absorbance") #plus sign
points(Data$Hrs_Approx[Data$Treatment == "CA2"], Data$AbsorbAvg[Data$Treatment == "CA2"], pch = 18)
points(Data$Hrs_Approx[Data$Treatment == "T4"], Data$AbsorbAvg[Data$Treatment == "T4"], pch = 10, col = "red") #cirlc with cross
legend("topleft", inset = 0.05,legend = c( "Manure Control", "Antibiotic Control", "Antibiotic with Manure"), col = c("black", "black", "red"), pch = c(3,18,10), cex = 0.8, box.lty = 1)


# Plots for 0.45um filter pore
plot(Data$Hrs_Approx[Data$Treatment == "CM"], Data$AbsorbAvg[Data$Treatment == "CM"], pch = 3, xlab = "Hours", ylab = "Absorbance", ylim =c(.01,.35), main = "0.45um Filter") #plus sign
points(Data$Hrs_Approx[Data$Treatment == "T1"], Data$AbsorbAvg[Data$Treatment == "T1"], pch = 1, col = "red") #open circle
points(Data$Hrs_Approx[Data$Treatment == "T2"], Data$AbsorbAvg[Data$Treatment == "T2"], pch = 10, col = "blue") #cirlc with cross
points(Data$Hrs_Approx[Data$Treatment == "CA"], Data$AbsorbAvg[Data$Treatment == "CA"], pch = 18)

plot(Data$Hrs_Approx[Data$Treatment == "CA"], Data$AbsorbAvg[Data$Treatment == "CA"], pch = 18, ylim = c(0,.15), xlab = "Hours", ylab = "Absorbance")
points(Data$Hrs_Approx[Data$Treatment == "T3"], Data$AbsorbAvg[Data$Treatment == "T3"], pch = 1, col = "red")
points(Data$Hrs_Approx[Data$Treatment == "CM2"], Data$AbsorbAvg[Data$Treatment == "CM2"], pch = 3)


#Plot the change in shape of the absorbtion spectra over the 24 hours

par(mfrow = c(2,2))
par(mar = c(2.5,2,1,0.5))
#1hr
plot(OneHour$Wavelength_nm, OneHour$CM2.5_1hr, type = "l", ylim = c(.04,.4), xlim = c(200,500), main = "1 Hour")
lines(OneHour$Wavelength_nm, OneHour$CM2.2_1hr)
lines(OneHour$Wavelength_nm, OneHour$CA2.53_1hr, col = "grey")
lines(OneHour$Wavelength_nm, OneHour$T4.1_1hr, col = "red")

#3hr
plot(ThreeHour$Wavelength_nm, ThreeHour$T4.1_3hr, col = "red", type = "l", ylim = c(.04,.4), xlim = c(200,500), main = "3 Hours")
lines(ThreeHour$Wavelength_nm, ThreeHour$CM2.5_3hr)
lines(ThreeHour$Wavelength_nm, ThreeHour$CM2.3_3hr)
lines(ThreeHour$Wavelength_nm, ThreeHour$CM2.1_3hr)
lines(ThreeHour$Wavelength_nm, ThreeHour$CA2.1_3hr, col = "grey")
lines(ThreeHour$Wavelength_nm, ThreeHour$CA2.3_3hr, col = "grey")
lines(ThreeHour$Wavelength_nm, ThreeHour$CA2.5_3hr, col = "grey")
lines(ThreeHour$Wavelength_nm, ThreeHour$T4.3_3hr, col = "red")
lines(ThreeHour$Wavelength_nm, ThreeHour$T4.5_3hr, col = "red")

#8hr
plot(EightHour$Wavelength_nm, EightHour$T4.1_8hr, col = "red", type = "l", ylim = c(.04,.4), xlim = c(200,500), main = "8 Hours")
lines(EightHour$Wavelength_nm, EightHour$CM2.5_8hr)
lines(EightHour$Wavelength_nm, EightHour$CM2.3_8hr)
lines(EightHour$Wavelength_nm, EightHour$CM2.1_8hr)
lines(EightHour$Wavelength_nm, EightHour$CA2.1_8hr, col = "grey")
lines(EightHour$Wavelength_nm, EightHour$CA2.3_8hr, col = "grey")
lines(EightHour$Wavelength_nm, EightHour$CA2.5_8hr, col = "grey")
lines(EightHour$Wavelength_nm, EightHour$T4.3_8hr, col = "red")
lines(EightHour$Wavelength_nm, EightHour$T4.5_8hr, col = "red")

#24hr
plot(TwentyFourHour$Wavelength_nm, TwentyFourHour$T4.1_24hr, col = "red", type = "l", ylim = c(.04,.4), xlim = c(200,500), main = "24 Hours")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CM2.5_24hr)
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CM2.3_24hr)
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CM2.1_24hr)
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CA2.1_24hr, col = "grey")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CA2.3_24hr, col = "grey")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CA2.5_24hr, col = "grey")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$T4.3_24hr, col = "red")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$T4.5_24hr, col = "red")

#48hr
plot(FortyEightHour$Wavelength_nm, FortyEightHour$T4.1_48hr, col = "red", type = "l", ylim = c(.04,.4), xlim = c(200,500), main = "48 Hours")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CM2.5_48hr)
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CM2.3_48hr)
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CM2.1_48hr)
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CA2.1_48hr, col = "grey")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CA2.3_48hr, col = "grey")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CA2.2_48hr, col = "grey") # CA2.5 did not have enought sample left to run SPAN, therefore took CA2.2 instead
lines(FortyEightHour$Wavelength_nm, FortyEightHour$T4.3_48hr, col = "red")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$T4.5_48hr, col = "red")

#Time lapse for each treatment

par(mfrow = c(2,2))
par(mar = c(2.5,2,1,0.5))
#CM
plot(OneHour$Wavelength_nm, OneHour$CM2.5_1hr, type = "l", ylim = c(0,.3), xlim = c(200,500), main = "Manure Control")
lines(OneHour$Wavelength_nm, OneHour$CM2.2_1hr)
lines(ThreeHour$Wavelength_nm, ThreeHour$CM2.5_3hr, col = "blue")
lines(ThreeHour$Wavelength_nm, ThreeHour$CM2.3_3hr, col = "blue")
lines(ThreeHour$Wavelength_nm, ThreeHour$CM2.1_3hr, col = "blue")
lines(EightHour$Wavelength_nm, EightHour$CM2.5_8hr, col = "green")
lines(EightHour$Wavelength_nm, EightHour$CM2.3_8hr, col = "green")
lines(EightHour$Wavelength_nm, EightHour$CM2.1_8hr, col = "green")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CM2.5_24hr, col = "orange")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CM2.3_24hr, col = "orange")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CM2.1_24hr, col = "orange")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CM2.5_48hr, col = "red")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CM2.3_48hr, col = "red")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CM2.1_48hr, col = "red")

#T4
plot(ThreeHour$Wavelength_nm, ThreeHour$T4.1_3hr, col = "blue", type = "l", ylim = c(0,.3), xlim = c(200,500), main = "Antibiotic + Manure")
lines(OneHour$Wavelength_nm, OneHour$T4.1_1hr)
lines(ThreeHour$Wavelength_nm, ThreeHour$T4.3_3hr, col = "blue")
lines(ThreeHour$Wavelength_nm, ThreeHour$T4.5_3hr, col = "blue")
lines(EightHour$Wavelength_nm, EightHour$T4.1_8hr, col = "green")
lines(EightHour$Wavelength_nm, EightHour$T4.3_8hr, col = "green")
lines(EightHour$Wavelength_nm, EightHour$T4.5_8hr, col = "green")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$T4.3_24hr, col = "orange")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$T4.5_24hr, col = "orange")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$T4.1_24hr, col = "orange")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$T4.3_48hr, col = "red")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$T4.5_48hr, col = "red")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$T4.1_48hr, col = "red")

#CA2
plot(OneHour$Wavelength_nm, OneHour$CA2.53_1hr, type = "l", ylim = c(0,.3), xlim = c(200,500), main = "Antibiotic Control")
lines(ThreeHour$Wavelength_nm, ThreeHour$CA2.1_3hr, col = "blue")
lines(ThreeHour$Wavelength_nm, ThreeHour$CA2.3_3hr, col = "blue")
lines(ThreeHour$Wavelength_nm, ThreeHour$CA2.5_3hr, col = "blue")
lines(EightHour$Wavelength_nm, EightHour$CA2.1_8hr, col = "green")
lines(EightHour$Wavelength_nm, EightHour$CA2.3_8hr, col = "green")
lines(EightHour$Wavelength_nm, EightHour$CA2.5_8hr, col = "green")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CA2.1_24hr, col = "orange")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CA2.3_24hr, col = "orange")
lines(TwentyFourHour$Wavelength_nm, TwentyFourHour$CA2.5_24hr, col = "orange")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CA2.1_48hr, col = "red")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CA2.3_48hr, col = "red")
lines(FortyEightHour$Wavelength_nm, FortyEightHour$CA2.2_48hr, col = "red")



