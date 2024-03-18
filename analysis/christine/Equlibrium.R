# Antibiotic Adsorption to Manure and Soil Project
# Christine B. Georgakakos
# February 19, 2021

# Purpose of Code: 
# Equlibrium Experiment

# Load data
setwd('~/Dropbox/Soil and Water/antibiotics')
Data <- read.csv('Experiments/Data/RawData/Erythromycin_20210402.csv') # Soil Equilibrium and Manure filtrate with Soil

# Calculate time elapsed for each sample from start to filter in hours

Data$StartDateTime <- as.Date(ISOdatetime(Data$StartYear, Data$StartMonth, Data$StartDay, Data$StartHour, Data$StartMin,0))
Data$StopDateTime <- as.Date(ISOdatetime(Data$StopYear, Data$StopMonth, Data$StopDay, Data$StopHour, Data$StopMin,0))
Data$StartHourTime<- Data$StartHour + Data$StartMin/60 # hours
Data$StopHourTime<- Data$StopHour + Data$StopMin/60 # hours

# Calculated hours elapsed
Data$DayElapsed <- as.numeric(Data$StopDateTime - Data$StartDateTime) # days elapsed

Data$HoursElapsed <- Data$DayElapsed*24 + Data$StopHourTime - Data$StartHourTime

# Round hours elapsed up, to use as x-axis on eventual plot
Data$RoundedHours <-ceiling(Data$HoursElapsed) # ceiling() rounds up

#Average the multiple microplate reads, and add the dual wavelengths together
for (i in 1:nrow(Data)){
  Data$Average450nm[i] <- mean(Data$Read1_450nm[i],Data$Read2_450nm[i],Data$Read3_450nm[i],Data$Read4_450nm[i], Data$Read5_450nm[i], na.rm = TRUE)
  Data$Average630nm[i] <- mean(Data$Read1_630nm[i],Data$Read2_630nm[i],Data$Read3_630nm[i], Data$Read4_630nm[i],Data$Read5_630nm[i], na.rm = TRUE)
}

# Add total absorbance column for both 450 and 630 nm if machine does not autoclaculate
Data$TotalAbsorbance <- Data$Average450nm + Data$Average630nm

# Calculate percent absorbance by normalizing by 'S0' the 0ppb standard. 
Data$PercentAbsorbance <- Data$TotalAbsorbance/Data$TotalAbsorbance[1]*100


## STANDARD CURVE ANALYSIS
# Plot standard curve
Standard<- Data[which(Data$SampleType == 'Standard'),]

plot( Standard$TotalAbsorbance, Standard$StdConcen_ppb, log = 'xy')

# Take only the linear range
StandardLinear <- Standard[which(Standard$StdConcen_ppb > 0.1 & Standard$StdConcen_ppb < 50),]

# log log: log(cocentration) = m* log(percentAbsorbance)+ b
StandardCurve<- lm(log(StandardLinear$StdConcen_ppb)~log(StandardLinear$TotalAbsorbance))
summary(StandardCurve)

modelIntercept <- StandardCurve$coefficients[1]
modelSlope<-StandardCurve$coefficients[2]

#calculate concentration from new curve
Data$MeasuredConcentration <- exp(modelSlope*log(Data$TotalAbsorbance) + modelIntercept)

Data$ActualConcentration_ppb <- Data$MeasuredConcentration *Data$Dillution # I don't think this needs to be in a loop

# Calculate Average Values between 3 replicates

# Average replicates and calculate max & min of x & y direction so can create double error bars. 
# x-variable = aqueous concentration, y-variable = q
Data$ActualConcentration_mgL <- Data$ActualConcentration_ppb/1000
Data$q_mgmg <- (667/1000 - Data$ActualConcentration_mgL)*0.150/(Data$Soil_g*1000) # total volume of batches = 150 mL 

library(dplyr)

AverageEqulibrium<- Data %>% group_by(RoundedHours, Treatment) %>% 
  summarize(ActualConcentration_ppb.mean = mean(ActualConcentration_ppb),
            ymin = min(ActualConcentration_ppb), 
            ymax = max(ActualConcentration_ppb), 
            HoursElapsed.mean = mean(HoursElapsed), 
            xmin = min(HoursElapsed), 
            xmax = max(HoursElapsed))

# Same as above, but for q instead of c on y axis. 
AverageEqulibrium<- Data %>% group_by(RoundedHours, Treatment) %>% 
  summarize(q_mgmg.mean = mean(q_mgmg),
            ymin = min(q_mgmg), 
            ymax = max(q_mgmg), 
            HoursElapsed.mean = mean(HoursElapsed), 
            xmin = min(HoursElapsed), 
            xmax = max(HoursElapsed))

# Run this for q not for C
AverageEqulibrium <- AverageEqulibrium[which(AverageEqulibrium$Treatment == 'S4'),] # q not calc for controls

# Remove entries that are not equilibrium experiment.

AverageEqulibrium <- AverageEqulibrium[which(AverageEqulibrium$Treatment == 'S4' | AverageEqulibrium$Treatment == 'AC'),]


# plot Concentration left  with error bars
ggplot(AverageEqulibrium, aes(x = HoursElapsed.mean)) + 
  geom_point(size = 8, aes(x = HoursElapsed.mean, y = ActualConcentration_ppb.mean, color = factor(Treatment, labels = c("Control", "Soil-only")))) + 
  labs(x="Elapsed Time [hr]", y = "Aqueous Concentration [ppb]") + 
  theme_bw() + # remove background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  # remove grid lines
  geom_errorbar(aes(ymin = ymin,ymax = ymax))  +  # this is currently the data range
  geom_errorbarh(aes(xmin = xmin,xmax = xmax, y = ActualConcentration_ppb.mean))+    # this is currently the data range
  scale_color_manual("Treatment", values=c( "#E69F00", "#56B4E9"))+
  scale_fill_discrete( labels = c("Control", "Soil-only"))+
  theme(text = element_text(size = 20))+
  geom_abline(slope = 0, intercept = 667, linetype = "longdash")

# plot q data with error bars
ggplot(AverageEqulibrium, aes(x = HoursElapsed.mean)) + 
  geom_point(aes(x = HoursElapsed.mean, y = q_mgmg.mean), size = 2) + 
  labs(x="Elapsed Time [hr]", y = "Adsorption density [mg/mg]") + 
  theme_bw() + # remove background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  # remove grid lines
  geom_errorbar(aes(ymin = ymin,ymax = ymax))  +  # this is currently the data range
  geom_errorbarh(aes(xmin = xmin,xmax = xmax, y = q_mgmg.mean)) +   # this is currently the data range 
  theme(text = element_text(size = 20))


FirstOrderKinetics = nls(q_mgmg.mean~q_e*(1-exp(-k*HoursElapsed.mean)), start = list(q_e = 1.0e-5, k=.1), data = AverageEqulibrium, control = nls.control(maxiter =100))

summary(FirstOrderKinetics)

# Identify equilibrium time statistically as the first sample that is not statistically different from both (t-1) and (t+1). 

EqulibriumSamples <- Data[which(Data$Treatment =='S4'),] 

EqulibTimes<-unique(EqulibriumSamples$RoundedHours)

EqulibTimes<- EqulibTimes[order(EqulibTimes)] # reorder in increasing order so that lines() works


for (i in 2:length(EqulibTimes)){
  TestTime <- EqulibriumSamples[which(EqulibriumSamples$RoundedHours ==EqulibTimes[i]),]
  PreviousTime <- EqulibriumSamples[which(EqulibriumSamples$RoundedHours ==EqulibTimes[i-1]),]
  NextTime <- EqulibriumSamples[which(EqulibriumSamples$RoundedHours ==EqulibTimes[i+1]),]
  print(c(EqulibTimes[i], 'hours'))
  PreviousTest<- wilcox.test(TestTime$ActualConcentration_ppb,PreviousTime$ActualConcentration_ppb, exact = TRUE)
  print(PreviousTest$p.value)
  NextTest<- wilcox.test(TestTime$ActualConcentration_ppb,NextTime$ActualConcentration_ppb, exact = TRUE)
  print(NextTest$p.value)
  PreviousNextTest<-wilcox.test(PreviousTime$ActualConcentration_ppb,NextTime$ActualConcentration_ppb, exact = TRUE)
  print(PreviousNextTest$p.value)
} 

# Are controls different from one another throughout experiment? Test first to last (12 hr, 96 hr). 

# Conclusion: not statistically different. 

ControlFirst<- Data[which(Data$Treatment == "AC"& Data$RoundedHours == 13),]
ControlLast<- Data[which(Data$Treatment == "AC"& Data$RoundedHours == 96),]

ControlTest<- wilcox.test(ControlFirst$ActualConcentration_ppb,ControlLast$ActualConcentration_ppb, exact = TRUE)
summary(ControlTest)
ControlTest

