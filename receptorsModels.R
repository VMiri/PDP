#receptors density maps correlations 
require(ggplot2)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(broom)
library(ggplot2)
require(MASS)

#import data
braindata <- read.csv(file.choose(), header = TRUE, sep=",")
View(braindata)
head(braindata)

# *** Rregression model + Cook's outlier detection for significantly diff regions as result of mancovas
braindata1 <- braindata[braindata$sigdiff_mancova == '1',]
head(braindata1)

model <- lm(diff1 ~ X5ht1a , data = braindata1)
model
summary(model)

# outliers detection 
cooksd <- cooks.distance(model)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, 
     labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])  # influential row numbers
head(braindata1[influential, ])  # influential observations.

braindata1.out <- braindata1[-c(), ] 
head(braindata1.out)
view(braindata1.out)

model2 <- lm(diff1 ~ receptor, data = braindata1.out)
model2
summary(model2)

view(braindata1.out)


ggplot(braindata1.out, aes(x=diff1, y= receptor)) + 
  geom_point()+
  geom_smooth(method=lm,  linetype="dashed",
              color="darkred", fill="blue")


#*** regression model for all regions 

model <- lm(diff1 ~ receptor, data = braindata)
model
summary(model)


ggplot(braindata, aes(x=diff1, y=receptor)) + 
  geom_point()+
  geom_smooth(method=lm,  linetype="dashed",
              color="darkred", fill="blue")


# outliers detection 
cooksd <- cooks.distance(model)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, 
     labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])  # influential row numbers
head(braindata[influential, ])  # influential observations.

braindata1.out <- braindata[-c(), ] 
head(braindata1.out)
view(braindata1.out)

# re-model if necessary 
model.out <- lm(diff1 ~ receptor  , data = braindata)
model.out
summary(model.out)



#model for non significant only 
view(braindata)
braindata1 <- braindata[braindata$sigdiff_mancova == '0',]
head(braindata1)
model <- lm(diff1 ~ receptor , data = braindata1)
model
summary(model)

# outliers detection 
cooksd <- cooks.distance(model)
plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, 
     labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels

influential <- as.numeric(names(cooksd)[(cooksd > 4*mean(cooksd, na.rm=T))])  # influential row numbers
head(braindata[influential, ])  # influential observations.

braindata1.out <- braindata[-c(), ] 
head(braindata1.out)
view(braindata1.out)

# re-model if necessary 
model.out <- lm(mean_diff ~ X5ht2a  , data = braindata1.out)
model.out
summary(model.out)
