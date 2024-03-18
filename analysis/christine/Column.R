# Antibiotic Adsorption to Manure and Soil Project
# Christine B. Georgakakos
# April 19, 2021

# Purpose of Code: 
# Column Experiment Analysis

# Load data
setwd('~/Dropbox/Soil and Water/antibiotics')

Data <- read.csv('Experiments/Data/RawData/Erythromycin_20210512.csv') # Column Test
Data <- read.csv('Experiments/Data/RawData/Erythromycin_20210514.csv') # Column Test 2
Data <- read.csv('Experiments/Data/RawData/Erythromycin_20210519.csv') # Soil-only Columns, 3 replicates
Data <- read.csv('Experiments/Data/RawData/Erythromycin_20210521.csv') # Soil-Manure Columns, 3 replicates

VolumeTimeData <- read.csv('Experiments/Data/RawData/ColumnTimeVolumeData.csv') # All column volume/time data: 6 columns, 3 replicates for soil-only and soil+manure

#Average the multiple reads, and add the dual wavelengths together
for (i in 1:nrow(Data)){
  Data$Average450nm[i] <- mean(Data$Read1_450nm[i],Data$Read2_450nm[i],Data$Read3_450nm[i],Data$Read4_450nm[i], Data$Read5_450nm[i], na.rm = TRUE)
  Data$Average630nm[i] <- mean(Data$Read1_630nm[i],Data$Read2_630nm[i],Data$Read3_630nm[i], Data$Read4_630nm[i],Data$Read5_630nm[i], na.rm = TRUE)
}


# Add total absorbance column for both 450 and 630 nm if machine does not autoclaculate
Data$TotalAbsorbance <- Data$Average450nm + Data$Average630nm

# Calculate percent absorbance by normalizing by 'S0' the 0ppb standard. 
Data$PercentAbsorbance <- Data$TotalAbsorbance/Data$TotalAbsorbance[1]*100

# Plot standard curve
Standard<- Data[which(Data$SampleType == 'Standard'),]

plot( Standard$TotalAbsorbance, Standard$StdConcen_ppb, log = 'xy')

# Take only the linear range
StandardLinear <- Standard[which(Standard$StdConcen_ppb > 0.1 & Standard$StdConcen_ppb < 50 & Standard$StdConcen_ppb != 10),]

# Find equation of standard curve
# log log: log(cocentration) = m* log(percentAbsorbance)+ b
StandardCurve<- lm(log(StandardLinear$StdConcen_ppb)~log(StandardLinear$TotalAbsorbance))
summary(StandardCurve)

modelIntercept <- StandardCurve$coefficients[1]
modelSlope<-StandardCurve$coefficients[2]

#calculate Measured concentration from new curve
Data$MeasuredConcentration_ppb <- exp(modelSlope*log(Data$TotalAbsorbance) + modelIntercept)

# Calculate actual concentration using dillution ratio
Data$ActualConcentration_ppb <- Data$MeasuredConcentration_ppb *Data$Dillution # I don't think this needs to be in a loop

ColumnTest <- Data[which(Data$SampleType == 'Column'),]
Column1 <- Data[which(Data$ColumnNumber == 1),]
Column2 <- Data[which(Data$ColumnNumber == 2),]
Column3 <- Data[which(Data$ColumnNumber == 3),]
Column4 <- Data[which(Data$ColumnNumber == 4),] # columns 4, 5, 6 are manure + soil
Column5 <- Data[which(Data$ColumnNumber == 5),]
Column6 <- Data[which(Data$ColumnNumber == 6),]

# Bind all antibiotic data without standards etc for gg plotting. 
AllColumns <- rbind(Column1, Column2, Column3, Column4, Column5, Column6)


# Pull the cumulative volume passed through the column from the VolumeTimeData into the antibiotic data dataframe
for (i in 1:nrow(AllColumns)){
  TargetRow <- VolumeTimeData[which(VolumeTimeData$SampleName == AllColumns$SampleName[i]),]
  AllColumns$CumulativeVolume_mL[i] <-  TargetRow$CumulativeVolume_mL[1]
}



# write all column data to single CV
write.csv(AllColumns, 'Experiments/Data/CalculatedData/ColumnExperiment.csv')


AllColumns<- read.csv('Experiments/Data/CalculatedData/ColumnExperiment.csv')
AllColumns$ColumnNumber <- as.character(AllColumns$ColumnNumber)

ggplot(AllColumns, aes(x = CumulativeVolume_mL)) + 
  geom_point(aes(x = CumulativeVolume_mL, y = ActualConcentration_ppb, color = SolidType, shape = ColumnNumber, size = 2)) + 
  labs(x="Cumulative Volume [mL]", y = "Aqueous Concentration [ppb]") + 
  theme_bw() + # remove background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+# remove grid lines
  geom_smooth(aes(x=CumulativeVolume_mL, y = ActualConcentration_ppb, color = SolidType), )+
  theme(text = element_text(size = 20)) # increase font size




### For plotting individual columns. 
par(mfrow=c(3,1))
par(mar=c(4,4,1,1))

AllColumns$ColumnNumber <- as.numeric(AllColumns$ColumnNumber)
for (j in 1:6){
  OneColumn <- AllColumns[which(AllColumns$ColumnNumber == j),]
  plot(OneColumn$CumulativeVolume_mL,OneColumn$ActualConcentration_ppb, main = j, xlab = "Cumulative Volume [mL]", ylab = "Aqueous Concentration [ppb]", ylim = c(0,500))
}

#plot(ColumnTest$SampleNumber ,ColumnTest$ActualConcentration_ppb, ylab = 'Aqueous Concentration', xlab = 'Sample Number', pch =20, cex = 2, main = 'Test', ylim=c(0,400))
plot(Column1$SampleNumber ,Column1$ActualConcentration_ppb, ylab = 'Aqueous Concentration', xlab = 'Sample Number', pch =20, cex = 2, main = '1',  ylim=c(0,400))
plot(Column2$SampleNumber ,Column2$ActualConcentration_ppb, ylab = 'Aqueous Concentration', xlab = 'Sample Number', pch =20, cex = 2, main = '2',  ylim=c(0,400))
plot(Column3$SampleNumber ,Column3$ActualConcentration_ppb, ylab = 'Aqueous Concentration', xlab = 'Sample Number', pch =20, cex = 2, , main = '3',  ylim=c(0,400))

plot(Column4$SampleNumber ,Column4$ActualConcentration_ppb, ylab = 'Aqueous Concentration', xlab = 'Sample Number', pch =20, cex = 2, main = '4',  ylim=c(0,400))
plot(Column5$SampleNumber ,Column5$ActualConcentration_ppb, ylab = 'Aqueous Concentration', xlab = 'Sample Number', pch =20, cex = 2, main = '5',  ylim=c(0,400))
plot(Column6$SampleNumber ,Column6$ActualConcentration_ppb, ylab = 'Aqueous Concentration', xlab = 'Sample Number', pch =20, cex = 2, , main = '6',  ylim=c(0,400))


VolumeTimeData$CollectionTime_min<-as.numeric(VolumeTimeData$CollectionTime_min)
VolumeTimeData$SampleVolume_mL<-as.numeric(VolumeTimeData$SampleVolume_mL)

### Calculate flow rates ###
for (k in 1:6){
  ColumnSubset<- VolumeTimeData[which(VolumeTimeData$ColumnNumber == k),]
  TotalCollectionTime <-sum(ColumnSubset$CollectionTime_min, na.rm = TRUE)
  TotalCalcFlowRate <- sum(ColumnSubset$SampleVolume_mL, na.rm=TRUE) / TotalCollectionTime
  print(c('Column', k, 'TotalCalcFlowRate=', TotalCalcFlowRate))
  
}




