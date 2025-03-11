library(tidyverse) # data wrangling and Exploratory data analysis (EDA)
library(dplyr)
library(MASS) # EDA
library(corrr) # EDA
library (corrplot) # VIF plot
library(naniar) # EDA # checking for missing values 
library(MuMIn)  # model selection
library(AER)
library(ggplot2) # for producing graphs
library(visreg)  # model plotting 
library(GGally)  # for producing graphs
library(ggtext) # for producing graphs
library(ggpubr) # for producing graphs
library(confintr) # Cramer's V
library(epiR) # odds ratio calculation
library(pwr)  # power, effect size and sample size calculation
library(MuMIn) # model selection
library(parameters) # model visualization
library(modelsummary) # smooth model summary output
library(DHARMa) # model check
##############################################################################

##############################################################################
## STEP ONE ##
##############################################################################
# data = Nest_Data %>%
#   filter(`SpeciesB` == "Little Cormorant") %>%
#   dplyr::select(-c(`SpeciesB`, `Plant`, `District`, `PlotNo`,
#                    `GPSN`, `GPSE`)) %>%
#   filter(`PlantN` > 0) %>%
#   filter(NestN <= 100)
# save(data, file = "Cormorant.RData")

# load R data file
load("/Users/lynx025/Desktop/Cormorant.RData")
# check for NA values
any(is.na(data))
# remove rows with NA values 
data <- na.omit(data)

##############################################################################
## STEP TWO: DATA VISUALISATION ##
##############################################################################
# boxplot (data$NestN, main = 'Nest No.')
# boxplot(data$PlantN, main = "Total Plants No.")
# boxplot(data$PSpecies, main = "Plants with Species")
# boxplot(data$HeightA, main = "Average Tree Height")
# boxplot(data$HeightL, main = "Lowest Height")
# boxplot(data$HeightH, main = "Highest Height")
# 
# plot(data$NestN, data$HeightA)
# plot(data$NestN, data$HeightL)
# plot(data$NestN, data$HeightH)
# plot(data$NestN, data$PlantN)
# plot(data$NestN, data$PSpecies)
# 
# ggplot(data, aes(x = Region, y = NestN)) +
#   geom_boxplot()
# 
# ggplot(data, aes(x = PType, y = NestN)) +
#   geom_boxplot()


##############################################################################
## STEP THREE: CORRELATION and VIF CHECKING ##
##############################################################################
# data = data %>% dplyr::select(PlantN,
#                               PSpecies,
#                               HeightA, 
#                               HeightL, 
#                               HeightH)
# # two correlated variables can not be in the same model
# corD = cor(data) # now, run it
# colnames(corD) = colnames(data)
# rownames (corD) = colnames(data) # save it as dataframe
# write.csv(corD, 'matrix1.csv')
# 
# # VIF values
# corrD = ginv(cor(data)) # ginv is a function of MASS
# colnames(corrD) = colnames(data)
# rownames (corrD) = colnames(data)
# corrplot(corr = corrD, method = "number", is.corr= F)


##############################################################################
## STEP FOUR: MODEL RUN ##
##############################################################################
# Poisson Regression
model1  =  glm(NestN ~ PlantN + PType + Region + HeightH, 
               family = poisson(link="log"), data = data)

step(model1, direction = 'backward')
parameters::parameters(model1) 
DHARMa::testDispersion(model1) 
DHARMa::testZeroInflation(model1) 
DHARMa::dispersiontest(model1) 

# NOTE:
# overdisperion exists
# lets try with negative binomal

model_nb = glm.nb(NestN ~ PlantN + PType + Region + HeightA, data = data)
summary(model_nb)$theta
parameters(model_nb) # package parameters
testDispersion(model_nb) # package DHARMa
testZeroInflation(model_nb) # package DHARMa # no zero inflation 
# dispersiontest(model_nb) # package AER # not necessary as overdispersion is accounted for


##############################################################################
## STEP FIVE: MODEL SELECTION (RELATIVE) ##
##############################################################################
step(model_nb, direction = 'backward')
nb_residuals <- simulateResiduals(model_nb)
plot(nb_residuals)

data$PType <- as.factor(data$PType)
data$PType <- relevel(data$PType, 
                      ref = "Native") # for organising model summary 
                                      # table
 
# likelihood-ratio tests for models below 2 delta AIC
model2  =  glm(NestN ~
                 PlantN + PType + Region + HeightA,
               family = negative.binomial(theta = 1.73, link="log"), data = data)
step(model2, direction = 'backward')

model3 = glm(NestN ~ PlantN + PType + HeightA,
             family = 
               negative.binomial(theta = 1.73, 
                                 link="log"), 
             data = data)

model4 = glm(NestN ~ PlantN +  HeightA,
             family = 
               negative.binomial(theta = 1.73, 
                                 link="log"), 
             data = data)
model5 = glm(NestN ~ 1,
             family = 
               negative.binomial(theta = 1.73, 
                                 link="log"), 
             data = data)

# perform likelihood ratio tests
anova(model3, model2, test = "Chisq")  # testing model2 vs model3
anova(model4, model3, test = "Chisq")  # testing model3 vs model4
anova(model4, model2, test = "Chisq")  # testing model2 vs model4 (full vs. most reduced)

# AIC and BIC values
AIC(model2, model3, model4, model5)
BIC(model2, model3, model4, model5)

# check modelsummary 
modelsummary::modelsummary(list(model2, model3, model4, model5), 
             stars = TRUE, 
             statistic = c("estimate", "std.error", "p.value"))
parameters::parameters(model2)
parameters::parameters(model3)
parameters::parameters(model4)

##############################################################################
## STEP SIX: MODEL SELECTION (ABSOLUTE) ##
##############################################################################
residuals <- simulateResiduals(model2)
plot(residuals)
plotResiduals(residuals)
# no quantile deviation detected

residuals <- simulateResiduals(model3)
plot(residuals)
plotResiduals(residuals)
# quantile deviation detected

residuals <- simulateResiduals(model4)
plot(residuals)
plotResiduals(residuals)
# more quantile deviation detected

# NOTE
# model 2 (full model) should be used
# for inference

# more checking on model 2
# checking autocorrelation
acf(residuals(model2)) # independence is OK
# checking overdispersion and zero-inflation
parameters::parameters(model2) 
DHARMa::testDispersion(model2) # no dispersion
DHARMa::testZeroInflation(model2) # no zero inflation
DHARMa::testOutliers(model2) # no outlier
# dispersiontest(model_nb) # package AER 
#                          # not necessary as overdispersion 
                           # is accounted for
# plot(model2)

##############################################################################
## STEP SEVEN: PLOTTING THE EFFECT ##
##############################################################################
visreg(model_nb, "PlantN", scale = "response", rug= 0)

# Set the font to Times New Roman, font size to 14, and color to gray50
par(family = "Times New Roman", cex.axis = 1.2, col.axis = "gray50", cex.lab = 1.2, col.lab = "gray50")

?visreg
citation("DHARMa")


# Create the visreg plot
visreg(model_nb, "HeightA", scale = "response", rug = 0,
       ylab = "no. of nests per colony site", 
       xlab = "average tree height (m) per colony site",
       line=list(col="steelblue"),      # Change the curve line color
       fill=list(col="gray90")) # Change the confidence interval shade)

# Set the font to Times New Roman, font size to 14, and color to gray50
par(family = "Times New Roman", cex.axis = 1.6, 
    col.axis = "gray60", cex.lab = 1.4, 
    col.lab = "gray40",
    mgp = c(3, 0.9, 0))