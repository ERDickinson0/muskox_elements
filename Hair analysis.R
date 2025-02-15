###############################################################################
#################### Muskox elements and populations analyses #################
####################### March 2024 -  Eleanor Dickinson #######################
###############################################################################

################ Read in the data and load the packages #######################
rm(list=ls()) # clear the environment
setwd("C:/Users/elean/OneDrive - University of Calgary/Kutz Postdoc/Hair elements/Data analysis")

# load libraries
library(dplyr) 
library(devtools) 
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(ggthemes)
library(lme4)
library(ggfortify)
library(MuMIn)
library(reshape)
library(meantables)
library(car)
library(caret)
library(MASS)
library(klaR)
library(rrcov)
library(scales)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# read data
data <- read.csv("Pop_minerals_fulldata.csv")
head(data)
str(data)

##
data$Location <- as.factor(data$Location)
data$Sex <- as.factor(data$Sex)
data$Age <- as.factor(data$Age)
data$Trend <- as.factor(data$Trend)
data$Cd <- as.numeric(data$Cd)

#### Editing columns for consistency ####
# Age 
levels(data$Age)

for(i in 1:nrow(data)){
  if(data$Age[i] == "unknown"){
    data$Age[i] <- "Unknown"
  }else if(data$Age[i] == "0"){
    data$Age[i] <- "Unknown"
  }else if(data$Age[i] == "calf/yearling"){
    data$Age[i] <- "yearling"
  }else if(data$Age[i] == "18 mos."){
    data$Age[i] <- "yearling"
  }else if(data$Age[i] == "calf"){
    data$Age[i] <- "Calf"
  }else if(data$Age[i] == "6 mos."){
    data$Age[i] <- "Calf"
  }else {
    data$Age[i] <- "Adult"
  }
}
data$Age <- as.factor(data$Age)
levels(data$Age)
data$Age <- factor(data$Age, levels = c("Calf", "yearling", "Adult", "Unknown"))

## Sex ##
levels(data$Sex)

for(i in 1:nrow(data)){
  if(data$Sex[i] == "M"){
    data$Sex[i] <- "Male"
  }else if(data$Sex[i] == "male"){
    data$Sex[i] <- "Male"
  }else if(data$Sex[i] == "F"){
    data$Sex[i] <- "Female"
  }else if(data$Sex[i] == "female"){
    data$Sex[i] <- "Female"
  }else {
    data$Sex[i] <- "Unknown"
  }
}

## Trend and location ordering
data$Trend <- factor(data$Trend, levels = c("Increasing", "Stable", "Declining"))

data$Location <- factor(data$Location, 
                              levels = c("Eastern Hudson Bay",
                                         "North Great Slave", "Ungava Bay", 
                                         "Alaska Eastern North Slope",
                                         "Nunavut Mainland",
                                         "Seward Peninsula","Yukon North Slope",
                                         "Banks Island","East Victoria Island",
                                         "NW Victoria Island","Zackenberg"))
#### Removing outliers ####
data <- subset(data, Cu < 12)
data <- subset(data, Fe < 1050)
data <- subset(data, Ca < 2000)

# remove missing values
data <- data[!is.na(data$Cu),]
data <- data[!is.na(data$Se),]

########### A: GLMMs on age, sex and collection methods #########################
# Generalised linear mixed effects models were used to test the effect of;
# 1) collection method (harvest, capture, ground collection), and 
# 2) animal age class (calf, yearling, adult) and sex (male, female). 
# Population location and sampling year were included as nested random effects. 

# Age and sex was not known for all samples, these samples were removed from the 
# age/sex models
data2 <- subset(data, Age != "Unknown")# Subset to remove unknown age and sex
data2 <- subset(data2, Sex != "Unknown")
levels(data2$Location)

##### Sodium #######
## 1) Collection type
sodium1 <- lmer(Na ~ Collection_type + (1|Location/Year), data = data)
# Check model assumptions
plot(sodium1)
par(mfrow=c(1,2))
qqnorm(residuals(sodium1))
hist(residuals(sodium1))

# Model fit
summary(sodium1)
drop1(sodium1, test = "Chisq")
r.squaredGLMM(sodium1)
sodium1b <- lmer(Na ~(1|Location/Year), data = data)
anova(sodium1b, sodium1)

## 2) Age and sex
sodium2 <- lmer(Na ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(sodium2)
drop1(sodium2, test = "Chisq")
r.squaredGLMM(sodium2)

sodium2b <- lmer(Na ~ Sex + (1|Location/Year), data = data2)
anova(sodium2b, sodium2)
sodium2c <- lmer(Na ~ Age + (1|Location/Year), data = data2)
anova(sodium2c, sodium2)

##### Magnesium #######
## 1) Collection type
magnes1 <- lmer(Mg ~ Collection_type + (1|Location/Year), data = data)
# Check model assumptions
plot(magnes1)
par(mfrow=c(1,2))
qqnorm(residuals(magnes1))
hist(residuals(magnes1))

# Model fit
summary(magnes1)
drop1(magnes1, test = "Chisq")
r.squaredGLMM(magnes1)
magnes1b <- lmer(Mg ~(1|Location/Year), data = data)
anova(magnes1b, magnes1)

## 2) Age and sex
magnes2 <- lmer(Mg ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(magnes2)
drop1(magnes2, test = "Chisq")
r.squaredGLMM(magnes2)

magnes2b <- lmer(Mg ~ Sex + (1|Location/Year), data = data2)
anova(magnes2b, magnes2)
magnes2c <- lmer(Mg ~ Age + (1|Location/Year), data = data2)
anova(magnes2c, magnes2)


##### Calcium ######
## 1) Collection type
calcium1 <- lmer(Ca ~ Collection_type + (1|Location/Year), data = data)
# Check model assumptions
plot(calcium1)
par(mfrow=c(1,2))
qqnorm(residuals(calcium1))
hist(residuals(calcium1))

# Model fit
summary(calcium1)
drop1(calcium1, test = "Chisq")
r.squaredGLMM(calcium1)
calcium1b <- lmer(Ca ~(1|Location/Year), data = data)
anova(calcium1b, calcium1)

## 2) Age and sex
calcium2 <- lmer(Ca ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(calcium2)
drop1(calcium2, test = "Chisq")
r.squaredGLMM(calcium2)

calcium2b <- lmer(Ca ~ Sex + (1|Location/Year), data = data2)
anova(calcium2b, calcium2)
calcium2c <- lmer(Ca ~ Age + (1|Location/Year), data = data2)
anova(calcium2c, calcium2)


##### Chromium ######
## 1) Collection type
chromium1 <- lmer(Cr ~ Collection_type + (1|Location/Year), data = data)
# Check model assumptions
plot(chromium1)
par(mfrow=c(1,2))
qqnorm(residuals(chromium1))
hist(residuals(chromium1))

# Model fit
summary(chromium1)
drop1(chromium1, test = "Chisq")
r.squaredGLMM(chromium1)
chromium1b <- lmer(Cr ~(1|Location/Year), data = data)
anova(chromium1b, chromium1)

## 2) Age and sex
chromium2 <- lmer(Cr ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(chromium2)
drop1(chromium2, test = "Chisq")
r.squaredGLMM(chromium2)

chromium2b <- lmer(Cr ~ Sex + (1|Location/Year), data = data2)
anova(chromium2b, chromium2)
chromium2c <- lmer(Cr ~ Age + (1|Location/Year), data = data2)
anova(chromium2c, chromium2)


##### Manganese #######
## 1) Collection type
mangan1 <- lmer(Mn ~ Collection_type + (1|Location/Year), data = data)
# Check model assumptions
plot(mangan1)
par(mfrow=c(1,2))
qqnorm(residuals(mangan1))
hist(residuals(mangan1))

# Model fit
summary(mangan1)
drop1(mangan1, test = "Chisq")
r.squaredGLMM(mangan1)
mangan1b <- lmer(Mn ~(1|Location/Year), data = data)
anova(mangan1b, mangan1)

## 2) Age and sex
mangan2 <- lmer(Mn ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(mangan2)
drop1(mangan2, test = "Chisq")
r.squaredGLMM(mangan2)

mangan2b <- lmer(Mn ~ Sex + (1|Location/Year), data = data2)
anova(mangan2b, mangan2)
mangan2c <- lmer(Mn ~ Age + (1|Location/Year), data = data2)
anova(mangan2c, mangan2)


##### Iron #######
## 1) Collection type
iron1 <- lmer(Fe ~ Collection_type + (1|Location/Year), data = data)
# Check model assumptions
plot(iron1)
par(mfrow=c(1,2))
qqnorm(residuals(iron1))
hist(residuals(iron1))

# Model fit
summary(iron1)
drop1(iron1, test = "Chisq")
r.squaredGLMM(iron1)
iron1b <- lmer(Fe ~(1|Location/Year), data = data)
anova(iron1b, iron1)

## 2) Age and sex
iron2 <- lmer(Fe ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(iron2)
drop1(iron2, test = "Chisq")
r.squaredGLMM(iron2)

iron2b <- lmer(Fe ~ Sex + (1|Location/Year), data = data2)
anova(iron2b, iron2)
iron2c <- lmer(Fe ~ Age + (1|Location/Year), data = data2)
anova(iron2c, iron2)


##### Cobalt #######
## 1) Collection type
cobalt1 <- lmer(log(Co) ~  Collection_type +  Year +(1|Location), data = data)

# Check model assumptions
plot(cobalt1)
par(mfrow=c(1,2))
qqnorm(residuals(cobalt1))
hist(residuals(cobalt1))

# Model fit
summary(cobalt1)
drop1(cobalt1, test = "Chisq")
r.squaredGLMM(cobalt1)

## 2) Age and sex
cobalt2 <- lmer(Co ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(cobalt2)
drop1(cobalt2, test = "Chisq")
r.squaredGLMM(cobalt2)

cobalt2b <- lmer(log(Co) ~ Sex + (1|Location/Year), data = data2)
anova(cobalt2b, cobalt2)
cobalt2c <- lmer(log(Co) ~ Age + (1|Location/Year), data = data2)
anova(cobalt2c, cobalt2)


##### Copper #####
## 1) Collection type
copper1 <- lmer(Cu ~ Collection_type + (1|Location/Year), data = data)
# Check model assumptions
plot(copper1)
par(mfrow=c(1,2))
qqnorm(residuals(copper1))
hist(residuals(copper1))

# Model fit
summary(copper1)
drop1(copper1, test = "Chisq")
r.squaredGLMM(copper1)
copper1b <- lmer(Cu ~(1|Location/Year), data = data)
anova(copper1b, copper1)

## 2) Age and sex
copper2 <- lmer(Cu ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(copper2)
drop1(copper2, test = "Chisq")
r.squaredGLMM(copper2)

copper2b <- lmer(Cu ~ Sex + (1|Location/Year), data = data2)
anova(copper2b, copper2)
copper2c <- lmer(Cu ~ Age + (1|Location/Year), data = data2)
anova(copper2c, copper2)


##### Zinc #######
## 1) Collection type
zinc1 <- lmer(Zn ~ Collection_type + (1|Location/Year), data = data)
# Check model assumptions
plot(zinc1)
par(mfrow=c(1,2))
qqnorm(residuals(zinc1))
hist(residuals(zinc1))

# Model fit
summary(zinc1)
drop1(zinc1, test = "Chisq")
r.squaredGLMM(zinc1)
zinc1b <- lmer(Zn ~(1|Location/Year), data = data)
anova(zinc1b, zinc1)

## 2) Age and sex
zinc2 <- lmer(Zn ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(zinc2)
drop1(zinc2, test = "Chisq")
r.squaredGLMM(zinc2)

zinc2b <- lmer(Zn ~ Sex + (1|Location/Year), data = data2)
anova(zinc2b, zinc2)
zinc2c <- lmer(Zn ~ Age + (1|Location/Year), data = data2)
anova(zinc2c, zinc2)


##### Selenium ####
## 1) Collection type
selenium1 <- lmer(log(Se) ~  Collection_type +  Year +(1|Location), 
                  data = data)

# Check model assumptions
plot(selenium1)
par(mfrow=c(1,2))
qqnorm(residuals(selenium1))
hist(residuals(selenium1))

# Model fit
summary(selenium1)
drop1(selenium1, test = "Chisq")
r.squaredGLMM(selenium1)

## 2) Age and sex
selenium2 <- lmer(log(Se) ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(selenium2)
drop1(selenium2, test = "Chisq")
r.squaredGLMM(selenium2)

selenium2b <- lmer(log(Se) ~ Sex + (1|Location/Year), data = data2)
anova(selenium2b, selenium2)
selenium2c <- lmer(log(Se) ~ Age + (1|Location/Year), data = data2)
anova(selenium2c, selenium2)


##### Molybdenum #######
## 1) Collection type
molyb1 <- lmer(log(Mo) ~  Collection_type +  Year +(1|Location), data = data)

# Check model assumptions
plot(molyb1)
par(mfrow=c(1,2))
qqnorm(residuals(molyb1))
hist(residuals(molyb1))

# Model fit
summary(molyb1)
drop1(molyb1, test = "Chisq")
r.squaredGLMM(molyb1)

## 2) Age and sex
molyb2 <- lmer(Mo ~ Sex + Age + (1|Location/Year), data = data2)

# Model fit
summary(molyb2)
drop1(molyb2, test = "Chisq")
r.squaredGLMM(molyb2)

molyb2b <- lmer(log(Mo) ~ Sex + (1|Location/Year), data = data2)
anova(molyb2b, molyb2)
molyb2c <- lmer(log(Mo) ~ Age + (1|Location/Year), data = data2)
anova(molyb2c, molyb2)

par(mfrow=c(1,1))
########### B: MANOVA on elements and population trend #########################
# A MANOVA was used to test the association between elements and population 
# trend. 
data3 <- data[,c(1:6,9:13,15:24,26)]
data3 <- na.omit(data3)

res.man <- manova(cbind(data3$Cu, data3$Se, data3$Mo, data3$Co, data3$Zn, 
                        data3$Na, data3$Mn, data3$Mg, data3$Fe, data3$Cr, 
                        data3$Ca) ~ Trend + Collection_type, data = data3)
summary(res.man)
summary.aov(res.man)

Vars <- cbind(data3$Cu, data3$Se, data3$Mo, data3$Co, data3$Zn, data3$Na, 
              data3$Mn, data3$Mg, data3$Fe, data3$Cr, data3$Ca)
model <- lm(Vars ~ Trend, data = data3)

Manova(res.man, test.statistic = "Pillai")

Vars <- cbind(data3$Cu, data3$Se, data3$Co, data3$Zn, 
              data3$Mn, data3$Fe, data3$Ca)
model <- lm(Vars ~ Trend, data = data3)

Manova(model, test.statistic = "Pillai")
summary.aov(model)


########### C: Linear Discriminant Analysis (LDA) ##############################
# LDA is conducted to assign qiviut elements into specific dimensions defined 
# by population trend 

data3 <- data[,c(1:6,9:13,15:24,26)]
data3 <- na.omit(data3)
names(data3)
data3$Year <- as.numeric(data3$Year)

data3 <- subset(data3, Cu < 12)
data3 <- subset(data3, Fe < 1050)
data3 <- subset(data3, Ca < 2000)

## Order the locations
data3$Location <- factor(data3$Location, levels = c("Eastern Hudson Bay",
                     "North Great Slave","Ungava Bay",
                     "Alaska Eastern North Slope","Nunavut Mainland", 
                     "Seward Peninsula","Yukon North Slope","Banks Island",
                     "East Victoria Island","NW Victoria Island","Zackenberg"))

head(data3)
data3$Trend <- factor(data3$Trend, levels=c('Declining', 'Stable', "Increasing"))
data4 <- data3[,c(11,12:22)]


# Split the data into training (80%) and test set (20%)
set.seed(123)
training.samples <- data4$Trend %>%
  createDataPartition(p = 0.8, list = FALSE)
train.data <- data4[training.samples, ]
test.data <- data4[-training.samples, ]

# Normalise the data
# Estimate preprocessing parameters
preproc.param <- train.data %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(train.data)
test.transformed <- preproc.param %>% predict(test.data)
names(train.transformed)

# Fit the model
model <- lda(Trend~Cu+Se+Co+Zn+Fe+Ca+Mn, data = train.transformed)
# Make predictions
predictions <- model %>% predict(test.transformed)
# Model accuracy
mean(predictions$class==test.transformed$Trend)
table(predictions$class==test.transformed$Trend)

# Wilks Lamda criterion to assess model fit
Wilks.test(train.transformed[,-c(1,4,7,8,9,11)],factor(train.transformed$Trend))
formulaAll=Trend~Cu+Se+Mo+Co+Zn+Na+Mn+Mg+Fe+Cr+Ca
greedy.wilks(formulaAll, data = data4, niveau = 0.5) 

# View model
model
plot(model)

# Estimate the LD values
preds <- as.data.frame(cbind(as.character(predictions$class), predictions$x))
preds$V1 <- as.factor(as.character(preds$V1))
preds$LD1 <- as.numeric(preds$LD1)
preds$LD2 <- as.numeric(preds$LD2)

head(preds)
str(preds)
preds2 <- preds %>% 
  group_by(V1) %>% summarise(meanLD1 = mean(LD1), meanLD2 = mean(LD2))
preds2

#### Cross validation ####
model.cv <- lda(Trend~Cu+Se+Co+Zn+Fe+Ca+Mn, CV=TRUE, data = data4)
table(data4$Trend, model.cv$class, dnn = c('Actual Group','Predicted Group'))

confusion <- function(actual, predicted, names = NULL, printit = TRUE,
                      prior = NULL) {
  if (is.null(names))
    names <- levels(actual)
  tab <- table(actual, predicted)
  acctab <- t(apply(tab, 1, function(x) x/sum(x)))
  dimnames(acctab) <- list(Actual = names, "Predicted (cv)" = names)
  if (is.null(prior)) {
    relnum <- table(actual)
    prior <- relnum/sum(relnum)
    acc <- sum(tab[row(tab) == col(tab)])/sum(tab)
  }
  else {
    acc <- sum(prior * diag(acctab))
    names(prior) <- names
  }
  if (printit)
    print(round(c("Overall accuracy" = acc, "Prior frequency" = prior),
                4))
  if (printit) {
    cat("\nConfusion matrix", "\n")
    print(round(acctab, 4))
  }
  invisible(acctab)
}

confmat <- confusion(test.transformed$Trend, predictions$class, printit = FALSE)
fp <- (confmat[2,1]+confmat[3,1]+confmat[3,2])/3
fn <- (confmat[1,2]+confmat[1,3]+confmat[2,3])/3

#decline
tp <- confmat[1,1]
tn <- (confmat[2,2]+confmat[3,3])/2

sensiv1 <- tp/(tp+fn)
speciv1<- tn/(tn+fp)

#stable
tp <- confmat[2,2]
tn <- (confmat[1,1]+confmat[3,3])/2

sensiv2 <- tp/(tp+fn)
speciv2 <- tn/(tn+fp)

#increase
tp <- confmat[3,3]
tn <- (confmat[2,2]+confmat[1,1])/2

sensiv3 <- tp/(tp+fn)
speciv3 <- tn/(tn+fp)

# Average sensitivity and specificity
sensitivity <- (sensiv1 + sensiv2 + sensiv3 )/3
specificity <- (speciv1 + speciv2 + speciv3 )/3


########### D: Ordinal Logistic Regression ####################################
# To provide criteria to assess muskox population status from qiviut element 
# concentrations

# scale the variables
data3$Cu1 <- scale(data3$Cu)
data3$Se1 <- scale(data3$Se) 
data3$Fe1 <- scale(data3$Fe) 
data3$Zn1 <- scale(data3$Zn) 
data3$Ca1 <- scale(data3$Ca) 
data3$Mn1 <- scale(data3$Mn) 
data3$Co1 <- scale(data3$Co) 

### OLR
m <- polr(Trend ~ Cu1 + Se1 + Fe1 + Zn1 + Ca1 + Mn1 + Co1, data = data3)
Anova(m)
drop1(m, test = "Chisq")
summary(m)

ci <- confint(m, level = 0.95)
# OR and CI
round(exp(cbind(OR = coef(m), ci)),3)

m <- polr(Trend ~ Cu1 + Se1 + Fe1 + Zn1 + Ca1, data = data3)
# confidence intervals
round((Coef = coef(m)),3)
ci <- confint(m, level = 0.95)
# OR and CI
round(exp(cbind(OR = coef(m), ci)),3)


##### Copper
m <- polr(Trend ~ Cu, data = data3)
summary(m, digits = 3)
# confidence intervals
(ci <- confint(m, level = 0.95))
## OR and CI
round(exp(cbind(OR = coef(m), ci)),3)

#define new data frame that contains predictor variable
newdata <- data.frame(Cu=seq(min(data3$Cu), max(data3$Cu),len=1000))
pred <- predict(m,newdata,type = "p", conf.int=0.95)
categHi <- levels(data3$Trend)[max.col(pred)]# category with highest probability
head(categHi)

newdata <- cbind(newdata, pred, categHi)
head(newdata)
(lower <- max(newdata$Cu[newdata$categHi=="Declining"]))
(upper <- min(newdata$Cu[newdata$categHi=="Increasing"]))

#reshape the data
lnewdat <- melt(newdata, id.vars = c("Cu"),
                variable.name = "variable", value.name="value")
## view first few rows
head(lnewdat)
str(lnewdat)
lnewdat$value <- as.numeric(lnewdat$value)

## plots
ggplot(lnewdat, aes(x = Cu, y = value, colour = variable)) +
  geom_line() +theme_bw()

###
cu <- ggplot(data3, aes(x = Location, y = Cu, fill = Trend))+ 
  scale_fill_manual(values = c("#D81B60", "#FFC107","#1E88E5"))+
  geom_hline(yintercept=lower, linetype="dashed", color = "red", 
             linewidth = 0.8)+
  geom_hline(yintercept=upper, linetype="dashed", color = "blue", 
             linewidth = 0.8)+
  geom_boxplot()+
  theme_classic()+theme(panel.border =element_rect(fill=NA), 
                        text=element_text(family="serif"),
                        axis.text.x = element_blank(), 
                        axis.title.x = element_blank(),
                        legend.position = "none")+ 
  labs(y = "Copper (µg/g)", title = "a)", x = "Population trend")
cu

############## Selenium
m <- polr(Trend ~ Se, data = data3)
summary(m, digits = 3)
Anova(m)

# confidence intervals
(ci <- confint(m, level = 0.95))
## OR and CI
exp(cbind(OR = coef(m), ci))

hist(data3$Cu)

#define new data frame that contains predictor variable
newdata <- data.frame(Se=seq(min(data3$Se), max(data3$Se),len=1000))
pred <- predict(m,newdata,type = "p", conf.int=0.95)
categHi <- levels(data3$Trend)[max.col(pred)]   # category with highest value
head(categHi)

newdata <- cbind(newdata, pred, categHi)
head(newdata)
(lower <- (max(newdata$Se[newdata$categHi=="Declining"])))
(upper <- (min(newdata$Se[newdata$categHi=="Increasing"])))

#reshape the data
lnewdat <- melt(newdata, id.vars = c("Se"),
                variable.name = "variable", value.name="value")
## view first few rows
head(lnewdat)
lnewdat$value <- as.numeric(lnewdat$value)

## plots
ggplot(lnewdat, aes(x = Se, y = value, colour = variable)) +
  geom_line() 

####

se  <- ggplot(data3, aes(x = Location, y = Se, fill = Trend))+ 
  geom_hline(yintercept=lower, linetype="dashed", color = "red", 
             linewidth = 0.8)+
  geom_hline(yintercept=upper, linetype="dashed", color = "blue", 
             linewidth = 0.8)+
  geom_boxplot()+
  scale_fill_manual(values = c("#D81B60", "#FFC107","#1E88E5"))+  
  theme_classic()+theme(panel.border =element_rect(fill=NA), 
                        text=element_text(family="serif"),
                        axis.text.x = element_text(angle=55,hjust=1), 
                        legend.position = "none")+ 
  labs(y = "Selenium (µg/g)",  title = "b)", x = "Population")
se

############## Cobalt
m <- polr(Trend ~ Co, data = data3)
Anova(m)

# confidence intervals
(ci <- confint(m, level = 0.95))
## OR and CI
exp(cbind(OR = coef(m), ci))

#define new data frame that contains predictor variable
newdata <- data.frame(Co=seq(min(data3$Co), max(data3$Co),len=1000))
pred <- predict(m,newdata,type = "p", conf.int=0.95)
categHi <- levels(data3$Trend)[max.col(pred)]   # category with highest value
head(categHi)

newdata <- cbind(newdata, pred, categHi)
head(newdata)
(lower <- max(newdata$Co[newdata$categHi=="Declining"]))
(upper <- min(newdata$Co[newdata$categHi=="Increasing"]))

#reshape the data
lnewdat <- melt(newdata, id.vars = c("Co") , #"Se", "Cu", "Zn",  "Fe", "Ca"),
                variable.name = "variable", value.name="value")
## view first few rows
head(lnewdat)
lnewdat$value <- as.numeric(lnewdat$value)

## plots
ggplot(lnewdat, aes(x = (Co), y = value, colour = variable)) +
  geom_line() 

co <- ggplot(data3, aes(x = Location, y = log(Co+1), fill = Trend))+ 
  geom_boxplot()+
  scale_fill_manual(values = c("#D81B60", "#FFC107","#1E88E5"))+
  theme_classic()+theme(panel.border =element_rect(fill=NA), 
                        text=element_text(family="serif"),
                        axis.text.x = element_text(angle=55,hjust=1), 
                        legend.position = "none")+ 
  geom_text(x=11, y=18, label="P = 0.080", colour = "blue")+ 
  labs(y = "log Co (µg/g)", x = "", title = "e)") + coord_trans(y="log2")+ 
  scale_y_continuous( breaks = c(0.0009995003,0.004987542,0.009950331,0.02469261,
                                 0.09531018,0.4054651),
                      labels = c(0.001,0.005,0.01,0.025,0.1,0.5))
co

############## Zinc
m <- polr(Trend ~ Zn, data = data3)

# confidence intervals
(ci <- confint(m, level = 0.95))
## OR and CI
exp(cbind(OR = coef(m), ci))

#define new data frame that contains predictor variable
#define new data frame that contains predictor variable
newdata <- data.frame(Zn=seq(min(data3$Zn), max(data3$Zn),len=1000))
pred <- predict(m,newdata,type = "p", conf.int=0.95)
categHi <- levels(data3$Trend)[max.col(pred)]   # category with highest value
head(categHi)

newdata <- cbind(newdata, pred, categHi)
head(newdata)
lower <- min(newdata$Zn[newdata$categHi=="Declining"])
upper <- max(newdata$Zn[newdata$categHi=="Increasing"])

newdata <- data.frame(Zn=seq(min(data3$Zn), max(data3$Zn),len=100))
newdata <- cbind(newdata, predict(m, newdata, type = "probs"))
#reshape the data
lnewdat <- melt(newdata, id.vars = c( "Zn"),
                variable.name = "variable", value.name="value")
## view first few rows
head(lnewdat)

## plots
ggplot(lnewdat, aes(x = Zn, y = value, colour = variable)) +
  geom_line() 

zn <- ggplot(data3, aes(x = Location, y = Zn, fill = Trend))+ 
  geom_hline(yintercept=lower, linetype="dashed", color = "red",
             linewidth = 0.8)+
  geom_boxplot()+
  scale_fill_manual(values = c("#D81B60", "#FFC107","#1E88E5"))+
  theme_classic()+theme(panel.border =element_rect(fill=NA), 
                        text=element_text(family="serif"),
                        axis.text.x = element_text(angle=55,hjust=1), 
                        legend.position = "none")+
  scale_y_continuous(breaks = c(25, 50,75,100,125,150,175))+
  # geom_text(x=11, y=18, label="P = 0.080", colour = "blue")+ 
  labs(y = "Zn (µg/g)", x = "Population location", title = "a)")
zn


############## Manganese
m <- polr(Trend ~ Mn, data = data3)
Anova(m)

# confidence intervals
(ci <- confint(m, level = 0.95))
## OR and CI
exp(cbind(OR = coef(m), ci))

############
newdata <- data.frame(Mn=seq(min(data3$Mn), max(data3$Mn),len=100))
pred <- predict(m,newdata,type = "p", conf.int=0.95)
categHi <- levels(data3$Trend)[max.col(pred)]   # category with highest value
head(categHi)

newdata <- cbind(newdata, pred, categHi)
head(newdata)
lower <- max(newdata$Mn[newdata$categHi=="Declining"])
upper <- min(newdata$Mn[newdata$categHi=="Increasing"])

#reshape the data
lnewdat <- melt(newdata, id.vars = c("Mn"),
                variable.name = "variable", value.name="value")

lnewdat$value <- as.numeric(lnewdat$value)
## view first few rows
head(lnewdat)

## plots
ggplot(lnewdat, aes(x = Mn, y = value, colour = variable)) +
  geom_line() 

mn <- ggplot(data3, aes(x = Location, y = log(Mn+1), fill = Trend))+ 
  geom_boxplot()+
  scale_fill_manual(values = c("#D81B60", "#FFC107","#1E88E5"))+
  theme_classic()+theme(panel.border =element_rect(fill=NA), 
                        text=element_text(family="serif"),
                        axis.text.x = element_text(angle=55,hjust=1), 
                        legend.position = "none")+
  # geom_text(x=10.9, y=47, label="value = 0.41", colour = "blue")+ 
  labs(y = "log Mn (µg/g)", x = "Population location", title = "d)")+ 
  coord_trans(y="log2")+ 
  scale_y_continuous( breaks = c(0.2231436,0.4054651,0.6931472,1.252763,2.397895,
                                 3.931826),
                      labels = c(0.25,0.5,1,2.5,10,50))
mn
exp(4)-1
log(51)
log(1.1)
log(1.25)

############## Iron
m <- polr(Trend ~ Fe, data = data3)

# confidence intervals
(ci <- confint(m, level = 0.95))
## OR and CI
exp(cbind(OR = coef(m), ci))

######
newdata <- data.frame(Fe=seq(min(data3$Fe), max(data3$Fe),len=100))
pred <- predict(m,newdata,type = "p", conf.int=0.95)
categHi <- levels(data3$Trend)[max.col(pred)]   # category with highest value
head(categHi)

newdata <- cbind(newdata, pred, categHi)
head(newdata)
(lower <- min(newdata$Fe[newdata$categHi=="Declining"]))
(upper <- max(newdata$Fe[newdata$categHi=="Increasing"]))

#reshape the data
lnewdat <- melt(newdata, id.vars = c("Fe"),
                variable.name = "variable", value.name="value")
## view first few rows
head(lnewdat)
lnewdat$value <- as.numeric(lnewdat$value)

## plots
ggplot(lnewdat, aes(x = Fe, y = value, colour = variable)) +
  geom_line() 

fe <- ggplot(data3, aes(x = Location, y = log(Fe), fill = Trend))+ 
  geom_hline(yintercept=log(lower), linetype="dashed", color = "red", 
             linewidth = 0.8)+
  geom_hline(yintercept=log(upper), linetype="dashed", color = "blue", 
             linewidth = 0.8)+
  geom_boxplot()+
  scale_fill_manual(values = c("#D81B60", "#FFC107","#1E88E5"))+
  labs(y = "log Fe (µg/g)", x = "", title = "c)")+
  theme_classic()+theme(panel.border =element_rect(fill=NA), 
                        text=element_text(family="serif"),
                        axis.text.x = element_text(angle=55,hjust=1))+ 
  coord_trans(y="log2")+ 
  scale_y_continuous( breaks = c(1,1.609438,2.302585,3.218876,4.60517,6.214608, 
                                 7.600902),
                      labels = c(1, 5, 10, 25, 100, 500, 2000))

fe



############## Calcium
m <- polr(Trend ~ Ca, data = data3)

# confidence intervals
(ci <- confint(m, level = 0.95))
## OR and CI
exp(cbind(OR = coef(m), ci))

###########
newdata <- data.frame(Ca=seq(min(data3$Ca), max(data3$Ca),len=100))
pred <- predict(m,newdata,type = "p", conf.int=0.95)
categHi <- levels(data3$Trend)[max.col(pred)]   # category with highest value
head(categHi)

newdata <- cbind(newdata, pred, categHi)
head(newdata)
lower <- min(newdata$Ca[newdata$categHi=="Declining"])
upper <- min(newdata$Ca[newdata$categHi=="Increasing"])

#reshape the data
lnewdat <- melt(newdata, id.vars = c("Ca"),
                variable.name = "variable", value.name="value")

lnewdat$value <- as.numeric(lnewdat$value)
## view first few rows
head(lnewdat)

## plots
ggplot(lnewdat, aes(x = Ca, y = value, colour = variable)) +
  geom_line() 

ca <- ggplot(data3, aes(x = Location, y = Ca, fill = Trend))+ 
  geom_hline(yintercept=lower, linetype="dashed", color = "red", linewidth = 0.8)+
  geom_boxplot()+
  scale_fill_manual(values = c("#D81B60", "#FFC107","#1E88E5"))+
  theme_classic()+theme(panel.border =element_rect(fill=NA), 
                        text=element_text(family="serif"),
                        axis.text.x = element_text(angle=55,hjust=1), 
                        legend.position = "none")+
  labs(y = "Ca (µg/g)", x = "Population location", title = "b)")
ca

