# Data visualization using ggplot2


setwd("C:/Users/User/Documents/Bradford/R_session") # to set working directory
getwd()  # to see which is our working directory
Accessory<- read.csv("Accessory.csv", header = T, sep = ",") #reading the file
head(Accessory)#To see the first rows


data=as.matrix(Accessory[2:9]) #Selecting the rowns. as.matrix() function accepts numeric variables only.
head(data)

heatmap(data)

heatmap(data, scale="column") #scale option will scale variables.

heatmap(data, scale="column", cexCol = 0.75)

# using rcolorbrewer to change the color palette
library(RColorBrewer)
heatmap(data, scale="column", cexCol = 0.75, col= colorRampPalette(brewer.pal(8, "Reds"))(25))




