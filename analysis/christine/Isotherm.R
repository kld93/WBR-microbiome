# Antibiotic Adsorption to Manure and Soil Project
# Christine B. Georgakakos
# April 19, 2021

# Purpose of Code: 
# Isotherm analysis

setwd('~/Dropbox/Soil and Water/antibiotics')

# Load in Isotherm data
SoilIsotherm <- read.csv('Experiments/Data/CalculatedData/SoilIsotherm.csv')
SoilManureIsotherm <- read.csv('Experiments/Data/CalculatedData/SoilManureIsotherm.csv')


# Convert aqueous concentration into mg/L
SoilIsotherm$ActualConcentration_mgL <- SoilIsotherm$ActualConcentration_ppb/1000
SoilManureIsotherm$ActualConcentration_mgL <- SoilManureIsotherm$ActualConcentration_ppb/1000

# Calculate q (mg/mg) adsorbed 

SoilIsotherm$q_mgmg <- (SoilIsotherm$Initial_ppb/1000 - SoilIsotherm$ActualConcentration_mgL)*0.150/(SoilIsotherm$Soil_g*1000) # total volume of batches = 150 mL 

SoilManureIsotherm$q_mgmg <- (SoilManureIsotherm$Initial_ppb/1000 - SoilManureIsotherm$ActualConcentration_mgL)*0.150/((SoilManureIsotherm$Soil_g+SoilManureIsotherm$Manure_g)*1000) # total volume of batches = 150 mL 

#bind the manure and soil data into one data frame
CombinedIsotherm <- rbind(SoilIsotherm,SoilManureIsotherm)

#Remove negative adsorption values
for (i in 1:nrow(CombinedIsotherm)){
  if (is.na(CombinedIsotherm$q_mgmg[i]) == TRUE){
    CombinedIsotherm$q_mgmg[i] <- NA
  }else if (CombinedIsotherm$q_mgmg[i] < 0){
    CombinedIsotherm$q_mgmg[i] <- NA
  }
}

# Average replicates and calculate max & min of x & y direction so can create double error bars. 
# x-variable = aqueous concentration, y-variable = q

library(dplyr)

AverageIsotherm<- CombinedIsotherm %>% group_by(Treatment, SolidType, Initial_ppb) %>% 
  summarize(ActualConcentration_mgL.mean = mean(ActualConcentration_mgL),
            xmin = min(ActualConcentration_mgL), 
            xmax = max(ActualConcentration_mgL), 
            q_mgmg.mean = mean(q_mgmg), 
            ymin = min(q_mgmg), 
            ymax = max(q_mgmg))


# Define groups within the data 
# Run just one isotherm to evaulate models for Soil or Soil+Manure

#Soil
Isotherm <- AverageIsotherm[which(AverageIsotherm$SolidType == 'Soil' & is.na(AverageIsotherm$Initial_ppb) == FALSE & AverageIsotherm$q_mgmg.mean > 0),]

#Soil+Manure -- should RENAME this so it doens't overwrite things
Isotherm <- AverageIsotherm[which(AverageIsotherm$SolidType == 'SoilManure' & AverageIsotherm$Initial_ppb > 0 & AverageIsotherm$q_mgmg.mean > 0),]


#Controls
AntibioticControl <- AverageIsotherm[which(is.na(AverageIsotherm$SolidType) == TRUE),] 

SoilControl <- AverageIsotherm[which(AverageIsotherm$Treatment == 'SC'),]

# initial visualization
# Isotherm
plot(Isotherm$ActualConcentration_mgL.mean, Isotherm$q_mgmg.mean, ylab = 'q [mg/mg]', xlab = 'Aqueous Concentration [mg/L]')

plot(Isotherm$ActualConcentration_mgL.mean, Isotherm$q_mgmg.mean, log = 'xy', ylab = 'q [mg/mg]', xlab = 'Aqueous Concentration [mg/L]')





### Running models ####

### 1. Freundlich Isotherm ###
# q = k_f*c^n

Isotherm <- as.data.frame(Isotherm)

# Calculate parameters using non-linear regression: 
# plot and overlay with start guesses. 
#Soil:
Freundlich = nls(q_mgmg.mean~k_f*(ActualConcentration_mgL.mean^n), start = list(k_f = 1.0e-4, n=.1), data = Isotherm, control = nls.control(maxiter =100))

#Manure:
Freundlich = nls(q_mgmg.mean~k_f*(ActualConcentration_mgL.mean^n), start = list(k_f = 1.0e-5, n=5), data = Isotherm, control = nls.control(maxiter =100))

summary(Freundlich)

# Soil
k_f = 2.169e-4  # Insert new parameters from non linear model after model generation
n = 4.522e-1
c_mod2 <-seq(1,60,1)
q_mod2 <-k_f*(c_mod2)^n

# Manure
k_f =  1.873e-15 # Insert new parameters from non linear model after model generation; k_f not significant
n = 5.975 # (significant predictor, p-value = 0.0341)


### 2. Freundlich Liniearization  ###
# log10 q = a + b* log10 c
# q = 10^(a) * c^b

# Log transform the data (Note: log = ln in R; log10 = log base 10)
Isotherm$Logq <- log10(Isotherm$q_mgmg.mean)
Isotherm$Logc <- log10(Isotherm$ActualConcentration_mgL.mean)

plot(Isotherm$Logc, Isotherm$Logq, ylab = "log(q) [mg/g]", xlab = "log(c) [mg/L]")
FreundlichLM <- lm(Isotherm$Logq~Isotherm$Logc)

logc_mod3<- seq(-7,5,0.5)
logq_mod3<- FreundlichLM$coefficients[1] + logc_mod3*FreundlichLM$coefficients[2]
lines(logc_mod3,logq_mod3)


k_f2 <- 10^(FreundlichLM$coefficients[1])
n2<- FreundlichLM$coefficients[2]


### 3. Langmuir Isotherm ###

# 3A. Calculate using non-linear regression
# q = qmax*K*c/(1+K*c)
#Soil
LangmuirNLR = nls(Isotherm$q_mgmg.mean~qmax*K*(Isotherm$ActualConcentration_mgL.mean)/(1+K*Isotherm$ActualConcentration_mgL.mean), start = list(qmax = .001, K=0.04), data = Isotherm)
summary(LangmuirNLR)

#Manure
LangmuirNLR = nls(Isotherm$q_mgmg.mean~qmax*K*(Isotherm$ActualConcentration_mgL.mean)/(1+K*Isotherm$ActualConcentration_mgL.mean), start = list(qmax = 1e-3, K=10), data = Isotherm)
summary(LangmuirNLR)

# Save coefficients as variables for later plotting
# Soil
qmax <- 0.0015300 
K <- 0.0801433

# Manure
qmax <- 0.0002144
K <- -7.5902545

# Janurary Testing for odd shape. 

LangNLMTestC <- seq(.00000001,100,.5)
LangNLMTestC <- Isotherm$ActualConcentration_mgL.mean[1:8]
TESTQ <- Isotherm$q_mgmg.mean[1:8]


LangmuirNLRTEST = nls(TESTQ~qmax*K*(LangNLMTestC)/(1+K*LangNLMTestC), start = list(qmax = 1e-3, K=10))
summary(LangmuirNLRTEST)

TESTK <- 0.0006231
TESTQMAX <- 0.425193


LangNLMTestQ <- TESTQMAX*TESTK*LangNLMTestC / (1+ TESTK*LangNLMTestC)
plot((Isotherm$ActualConcentration_mgL.mean), (Isotherm$q_mgmg.mean))
lines((LangNLMTestC),(LangNLMTestQ))
  
  
  
# 3B. Linearize langmuir 1/c vs 1/q
# 1/q = a + b*1/c
# q = c/(a*c + b)

Isotherm$Cinverse <- 1/Isotherm$ActualConcentration_mgL.mean
Isotherm$Qinverse <- 1/Isotherm$q_mgmg.mean

plot(Isotherm$Cinverse, Isotherm$Qinverse, ylab = "1/q", xlab = "1/c")
LangmuirInvCInvQ<- lm(Isotherm$Qinverse~Isotherm$Cinverse)

cinv<- seq(0,400, 1)
qinv <- LangmuirInvCInvQ$coefficients[2]*cinv + LangmuirInvCInvQ$coefficients[1] 

lines(cinv, qinv)

qmax_inv <- 1/LangmuirInvCInvQ$coefficients[1]
K_inv<- 1/(qmax_inv*LangmuirInvCInvQ$coefficients[2])
  

# 3C. Linearize langmuir as c vs c/q
# c/q = a + b*c
# q = c/ (a+b*c)

Isotherm$CoverQ <- Isotherm$ActualConcentration_mgL.mean/Isotherm$q_mgmg.mean
plot(Isotherm$ActualConcentration_mgL.mean, Isotherm$CoverQ, ylab = "c/q", xlab = "c")
LangumirCCoverQ <- lm(Isotherm$CoverQ~Isotherm$ActualConcentration_mgL.mean)

c_mod3C<-seq(1,150,1)
coverq <- LangumirCCoverQ$coefficients[1]+ c_mod3C*LangumirCCoverQ$coefficients[2]
lines(c_mod3C,coverq)  

qmax_coverq <- 1/LangumirCCoverQ$coefficients[2]
K_coverq  <- 1/(qmax_coverq*LangumirCCoverQ$coefficients[1])


# Plot all models on the same graph to compare SSE

Isotherm$qFreundlichNLR <-k_f*(Isotherm$ActualConcentration_mgL.mean)^n
Isotherm$qFreundlichNLR[1:5] <- NA 
Isotherm$qFreundlichNLR[6:9] <-k_f*(Isotherm$ActualConcentration_mgL.mean[6:9])^n
Isotherm$qFreundlichLM <- k_f2*(Isotherm$ActualConcentration_mgL.mean)^n2
Isotherm$qLangmuirNLR <- qmax*K*(Isotherm$ActualConcentration_mgL.mean)/(1+K*Isotherm$ActualConcentration_mgL.mean)
Isotherm$qLangmuirLMinv <- 1/(LangmuirInvCInvQ$coefficients[1]+ LangmuirInvCInvQ$coefficients[2]/Isotherm$ActualConcentration_mgL.mean )
Isotherm$qLangmuirLMcoverq <- Isotherm$ActualConcentration_mgL.mean/(LangumirCCoverQ$coefficients[1] +LangumirCCoverQ$coefficients[2]*Isotherm$ActualConcentration_mgL.mean)


Isotherm<- Isotherm[order(Isotherm$Initial_ppb),] # reorder in increasing order so that lines() works

plot(Isotherm$ActualConcentration_mgL.mean,Isotherm$q_mgmg.mean, ylab = "q [mg/mg]", xlab = "c [mg/L]")
# plot(Isotherm$ActualConcentration_mgL,Isotherm$q_mgmg, ylab = "q [mg/mg]", xlab = "c [mg/L]", log = 'xy')
lines(Isotherm$ActualConcentration_mgL.mean, Isotherm$qFreundlichNLR, lty = 2)
lines(Isotherm$ActualConcentration_mgL.mean, Isotherm$qFreundlichLM, lty = 9, col = "red")
lines(Isotherm$ActualConcentration_mgL.mean, Isotherm$qLangmuirLMinv, lty = 4)
lines(Isotherm$ActualConcentration_mgL.mean, Isotherm$qLangmuirLMcoverq, lty = 1)
lines(Isotherm$ActualConcentration_mgL.mean, Isotherm$qLangmuirNLR, lty = 3)
legend("bottomright", inset = 0.05,legend = c( "Data", "Freundlich NLR", "Freundlich LM", "Langmuir LM [1/c vs 1/q]", "Langmuir LM [c vs c/q]", "Langmuir NLR"), col = c("black", "black", "red", "black", "black", "black"), lty = c(NA, 2,9,4,1,3), cex = 0.8, box.lty = 1)

# Calculate SSE = Σ (x-xhat)^2
FreundNLR_SSE = sum((Isotherm$q_mgmg.mean - Isotherm$qFreundlichNLR)^2)
FreundLM_SSE = sum((Isotherm$q_mgmg.mean - Isotherm$qFreundlichLM)^2)
LangmuirNLR_SSE = sum((Isotherm$q_mgmg.mean - Isotherm$qLangmuirNLR)^2)
LangmuirLMinv_SSE = sum((Isotherm$q_mgmg.mean - Isotherm$qLangmuirLMinv)^2)
LangmuirLMcoverq_SSE = sum((Isotherm$q_mgmg.mean - Isotherm$qLangmuirLMcoverq)^2)

# Save SSE values into a table for comparison. 
SSEValues <- data.frame(matrix( nrow = 5, ncol = 2))
SSEValues$Isotherm <- c("Freundlich NLR","Freundlich LM", "Langmuir NLR", "Langmuir LM (inv)", "Langmuir LM (c over q)")
SSEValues$SSE<- c(FreundNLR_SSE,FreundLM_SSE, LangmuirNLR_SSE, LangmuirLMinv_SSE, LangmuirLMcoverq_SSE)
SSEValues[,1:2]<- NULL

# model parameters into a table 
ModelParameters <- data.frame(matrix( nrow = 5, ncol = 1))
ModelParameters$Isotherm <- c("Freundlich NLR","Freundlich LM", "Langmuir NLR", "Langmuir LM (inv)", "Langmuir LM (c over q)")
ModelParameters$K <- c(k_f, k_f2, K, K_inv, K_coverq)
ModelParameters$NorQmax <- c(n, n2, qmax, qmax_inv, qmax_coverq)
ModelParameters[,1] <- NULL



### plotting using ggplot###

#install.packages('tidyverse')
#install.packages('ggplot2')
library("ggplot2")



ggplot(Isotherm, aes(x = ActualConcentration_mgL.mean)) + 
  geom_point(aes(x = ActualConcentration_mgL.mean, y = q_mgmg.mean)) + 
  labs(x="Aqueous Concentration [mg/L]", y = "q [mg/mg]") + 
  geom_line(aes(x=ActualConcentration_mgL.mean, y = qFreundlichNLR), linetype = "dashed") + 
  geom_line(aes(x=ActualConcentration_mgL.mean, y = qFreundlichLM), linetype = "dotted") + 
  geom_line(aes(x=ActualConcentration_mgL.mean, y = qLangmuirLMinv), linetype = "longdash") +
  geom_line(aes(x=ActualConcentration_mgL.mean, y = qLangmuirLMcoverq), linetype = "solid") +
  geom_line(aes(x=ActualConcentration_mgL.mean, y = qLangmuirNLR), linetype = "dotdash", color = "red") + 
  geom_errorbar(aes(ymin = ymin,ymax = ymax))  +  # this is currently the data range
  geom_errorbarh(aes(xmin = xmin,xmax = xmax, y = q_mgmg.mean)) +   # this is currently the data range
  theme_bw() + # remove background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  # remove grid lines
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') 


ggplot(Isotherm, aes(x = ActualConcentration_mgL.mean)) + 
    geom_point(aes(x = ActualConcentration_mgL.mean, y = q_mgmg.mean, size = 2)) + 
    labs(x="Aqueous Concentration [mg/L]", y = "Adsorption Density (q) [mg/mg]") + 
    geom_line(aes(x=ActualConcentration_mgL.mean, y = qFreundlichNLR, linetype = "Eq.S1", color = "Eq.S1")) + 
    geom_line(aes(x=ActualConcentration_mgL.mean, y = qFreundlichLM, linetype = "Eq.S2", color = "Eq.S2")) + 
  geom_line(aes(x=ActualConcentration_mgL.mean, y = qLangmuirNLR, linetype = "Eq.2", color = "Eq.2")) +
   geom_line(aes(x=ActualConcentration_mgL.mean, y = qLangmuirLMinv, linetype = "Eq.S3", color = "Eq.S3")) +
    geom_line(aes(x=ActualConcentration_mgL.mean, y = qLangmuirLMcoverq, linetype = "Eq.3", color = "Eq.3")) +
    geom_errorbar(aes(ymin = ymin,ymax = ymax))  +  # this is currently the data range
    geom_errorbarh(aes(xmin = xmin,xmax = xmax, y = q_mgmg.mean)) +   # this is currently the data range
    theme_bw() + # remove background color
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +  # remove grid lines
    #scale_x_continuous(trans='log10') +
    #scale_y_continuous(trans='log10') + 
    scale_color_manual("Model", values = c("Eq.S1" = "black", "Eq.S2" = "black","Eq.2" = "black", "Eq.S3" = "black", "Eq.3" = "red" ))+
    scale_linetype_manual("Model", values = c("Eq.S1" = 1, "Eq.S2" = 2, "Eq.2" = 5, "Eq.S3" = 3, "Eq.3" = 4))+
    # Soil - Langmuir NLR is best (eq.2), color = red; soil+manure - Langmuir c/q is best (eq.3), color = red
    theme(legend.direction = "vertical", legend.position = "right", legend.key = element_blank(), 
          legend.background = element_rect(fill = "white", colour = "gray30")) + 
    guides(fill = guide_legend(keywidth = 1, keyheight = 1), linetype=guide_legend(keywidth = 3, keyheight = 1),colour=guide_legend(keywidth = 3, keyheight = 1)) 
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') 

  
plot(AntibioticControl$Initial_ppb, AntibioticControl$ActualConcentration_mgL.mean)  

# both isotherms together, colored by solid type 
ggplot(AverageIsotherm, aes(x = ActualConcentration_mgL.mean)) + 
  geom_point(aes(x = ActualConcentration_mgL.mean, y = q_mgmg.mean, color = SolidType, size = 1.5)) + 
  labs(x="Equlibrium Aqueous Concentration (c) [mg/L]", y = "Adsorption (q) [mg/mg]") + 
  geom_errorbar(aes(ymin = ymin,ymax = ymax))  +  # this is currently the data range
  geom_errorbarh(aes(xmin = xmin,xmax = xmax, y = q_mgmg.mean)) +   # this is currently the data range
  theme_bw() + # remove background color
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  # remove grid lines
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme(text = element_text(size = 20))
  
# Estimate of required volume of solution for given soil mass
c = 500/1000 # ppb -> mg/L
m_soil <- 150*1000 # g -> mg

q = qmax*K*(c)/(1+K*c) # Langmuir NLR
q = k_f2*c^n2 # Freundlich LM 

m = q*m_soil
V = m/c




  
  