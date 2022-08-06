## Introduction

#In this project, [data set](https://archive.ics.uci.edu/ml/datasets/QSAR+fish+toxicity) containing information about fish toxicity was read in, some basic exploratory data #analysis (EDA - numeric summaries and graphs) done, data transformations/feature engineering, and then some linear regression models fitted. The tidyverse library was used  to #read in, manipulate, and graph the data set.

#This dataset was used to develop quantitative regression QSAR models to predict acute aquatic toxicity towards the fish Pimephales promelas (fathead minnow) on a set of 908 #chemicals. LC50 data, which is the concentration that causes death in 50% of test fish over a test duration of 96 hours,was used as model response. The model comprised 6 #molecular descriptors: MLOGP (molecularproperties), CIC0 (information indices), GATS1i (2D autocorrelations), NdssC (atom-type counts),NdsCH ((atom-type counts), SM1_Dz (2D #matrix-based descriptors).

#The columns are set up in the order 'CIC0', 'SM1_Dz', 'GATS1', 'NdsCH', 'NdssC', 'MLOGP', 'LC50' respectively and separated by a semicolon.

 ## Dataset

#install.packages("tidyverse")
#install.packages("readr")
library(readr)

#read in the dataset and display a few rows
data <- read.csv("C:\\Users\\junhw\\Desktop\\fp\\fishtoxic.csv", col.names =list('CIC0', 'SM1_Dz', 'GATS1', 'NdsCH', 'NdssC', 'MLOGP', 'LC50'), sep = ';')
head(data,10)#view head of dataset for first 10 rows

## Exploratory Data Analysis (EDA)

### Pairplot

#install.packages("GGally")
library(GGally)
#Plot Pairplot
ggpairs(data)

### Scatterplots of LC50 vs explanatory variables

#Scatterplot of LC50 vs explanatory variables
library(tidyr)
data %>%
  gather(-LC50, key = "some_var_name", value = "some_value_name") %>%
  ggplot(aes(x = some_value_name, y = LC50)) +
  geom_point() +
  facet_wrap(~ some_var_name, scales = "free")

### Plot of Correlation Matrix

#Get correlation and plot it on heatmap
correlationmatrix <- round(cor(data),2)
library(reshape2)
melted_matrix <- melt(correlationmatrix)
library(ggplot2)

# Get lower triangle of the correlation matrix
lower_tri<-function(correlationmatrix){
  correlationmatrix[upper.tri(correlationmatrix)] <- NA
  return(correlationmatrix)
}
# Get upper triangle of the correlation matrix
upper_tri <- function(correlationmatrix){
  correlationmatrix[lower.tri(correlationmatrix)]<- NA
  return(correlationmatrix)
}

upper_tri <- upper_tri(correlationmatrix)

# Melt the correlation matrix
library(reshape2)
melted_matrix <- melt(upper_tri, na.rm = TRUE)
# Heatmap
library(ggplot2)
correlheatmap <- ggplot(data = melted_matrix, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()



correlheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

### Boxplots

boxplot(data, out.width = '25%')
#out.width = '25%'
#out.height = '50%'

### Dataset Transformations

# - To investigate relationships between the variables in terms of their medians, the dataset had to be transformed. For each variable, a new binary variable was created
# The new variable was set to **'low'** if the original variable’s value is less than or equal to the median value of that variable, otherwise it was be set to **'high'**.
# - A function called **'onehotfunc'** was created to accomplish the transformation using the **ifelse()** or **if_else()** function.
# - The **apply()** function was used to apply the **onehotfunc** function to each column of the data frame.
# - The transformed result was then turned into a data frame and saved as an R object **Transformed_df**.

medians<-apply(data, 2, median)
medians<-data.frame(medians) 

CIC0 <- data["CIC0"]>medians[1,1]
SM1_Dz <- data["SM1_Dz"]>medians[2,1]
GATS1 <- data["GATS1"]>medians[3,1]
NdsCH <- data["NdsCH"]>medians[4,1]
NdssC <- data["NdssC"]>medians[5,1]
MLOGP<- data["MLOGP"]>medians[6,1]
LC50<- data["LC50"]>medians[7,1]
dataset<-data.frame(CIC0)
dataset['SM1_Dz']<-SM1_Dz 
dataset['GATS1']<-GATS1
dataset['NdsCH']<-NdsCH
dataset['NdssC']<-NdssC
dataset['MLOGP']<-MLOGP
dataset['LC50']<-LC50



onehotfuncy <- function(x){
  if (x==TRUE) {return ('high')} 
  else {return ('low')}
}
#Create a new dataframe with values less than median labelled as 'low' and greater than median labelled 'high'

final<-data.frame(apply(dataset,c(1,2),onehotfuncy))
#Save final dataframe
save(final, file = "Transformed_df.RData")

### Creating two way tables
# - The new binary variables from the final dataframe were used to create two two-way contingency table with corresponding side-byside/stacked bar plots visual using the binary LC50 variable and two other binary variables (**CIC0** and **GATS1**).

#create two way table from data frame
data2way1 <- table(final$LC50, final$CIC0)
#data2way1

#Plot first barplot for the two way table
barplot(data2way1, legend=TRUE, beside=TRUE, main='Distribution of LC50 & CICO')
#out.width='50%'

#create second two way table from data frame
data2way2 <- table(final$LC50, final$GATS1)
#data2way2

#Plot first barplot for the two way table
barplot(data2way2, legend=TRUE, beside=TRUE, main='Distribution of LC50 & GATS1')

## Multiple Linear Regression Models

#The goal was to predict a value of LC50 using the other variables in the data set (using the original data frame, not the binary version). Four different linear regression models were fitted. 
# - At least one model included a polynomial term *Model 2* and *Model 4*. 
# - At least one model included an interaction term *Model 3* and *Model4*.

# For each model, the summary() of the fit was displayed. For *model 4*,the diagnostic plots were displayed and comments on the model fit and the normality assumption were posted below the plots. Lastly, *model 2* and *Model 4* were used to predict the LC50 value. The LC50 value at the median setting of each predictor (explanatory) variable in the selected models were predicted.

# Multiple Linear Regression model 1
fit <- lm(LC50 ~ CIC0 + SM1_Dz + GATS1 + NdsCH +  NdssC + MLOGP, data=data)
# show results for model 1
summary(fit) 

# Multiple Linear Regression model 2 with interaction terms
fit2 <- lm(LC50 ~ CIC0 + SM1_Dz + GATS1 + NdsCH + NdssC + MLOGP + SM1_Dz*GATS1 + NdsCH*NdssC, data=data)
# show results for model 2
summary(fit2) 

# Multiple Linear Regression model 3 with polynomial for variable GATS1i (2D autocorrelations) and MLOGP (molecular properties)
fit3 <- lm(LC50 ~ CIC0 + SM1_Dz + GATS1 + NdsCH + NdssC + MLOGP + I(GATS1^2) + I(MLOGP^2), data=data)
# show results for model 3
summary(fit3) 

# Multiple Linear Regression model 4 with interactions and polynomial for variable GATS1i (2D autocorrelations) and MLOGP (molecular properties)
fit4 <- lm(LC50 ~ CIC0 + SM1_Dz + GATS1 + NdsCH + NdssC + MLOGP + SM1_Dz*GATS1 + CIC0*NdssC + I(GATS1^2) + I(MLOGP^2), data=data)
# show results for model 4
summary(fit4) 

### Diagnostic Plots for Linear Regression Analysis
#Diagnostic Plots for Linear Regression Analysis
par(mfrow=c(2,2)) # Change the panel layout to 2 x 2
plot(fit4)
par(mfrow=c(1,1)) # Change back to 1 x 1





### Comments

# Residuals vs Fitted**
  # - This plot shows if residuals have non-linear patterns. In this case, the residuals are equally spread around a horizontal line at the origin without distinct patterns which is a good indication that there are no non-linear relationships.

# Normal Q-Q**
  # - This plot shows if residuals are normally distributed. In this case, most residuals follow a straight line well with a slight deviation on both extremes. There are outliers on both extremes as shown on the plot.

# Scale-Location**
  # - This plot shows if residuals are spread equally along the ranges of predictors and helps check the assumption of equal variance (homoscedasticity). In this case, there is a horizontal line with equally (randomly) spread points which means that the equal variance assumption is not violated.

### Prediction

#Model 4 used for prediction
fit4 <- lm(LC50 ~ CIC0 + SM1_Dz + GATS1 + NdsCH + NdssC + MLOGP + SM1_Dz*GATS1 + CIC0*NdssC + I(GATS1^2) + I(MLOGP^2), data=data)

#Median values of the explanatory variables used for prediction
CIC0 <- medians[1,1]
SM1_Dz <- medians[2,1]
GATS1 <- medians[3,1]
NdsCH <- medians[4,1]
NdssC <- medians[5,1]
MLOGP <- medians[6,1]

prediction4 <- 2.65716 + 0.35813*CIC0+ 0.59775*SM1_Dz -0.90876*GATS1 + 0.40120 *NdsCH - 0.03661*NdssC + 0.29949*MLOGP + 0.50538*SM1_Dz*GATS1 +  0.02893*CIC0*NdssC - 0.06832*GATS1^2 + 0.02225*MLOGP^2

prediction4

#Model 1 used for prediction
fit <- lm(LC50 ~ CIC0 + SM1_Dz + GATS1 + NdsCH +  NdssC + MLOGP, data=data)


prediction1 <- 2.17415 + 0.38574*CIC0 + 1.25580 *SM1_Dz - 0.74628*GATS1 + 0.41349 *NdsCH + 0.06436 *NdssC + 0.39000*MLOGP 
prediction1