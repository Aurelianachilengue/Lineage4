getwd()
setwd("C:/Users/User/Documents/Bradford/R_session")

#Loading the file
AccessoryTable <- read.csv("AccessoryTable.csv", header = T, sep = ",")

# To vizualize the file
View(AccessoryTable)

#load reshape2 package to use melt() function
library(reshape2)

#melt AccessoryTable into long format
melt_AccessoryTable <- melt(AccessoryTable)

#add column for genomes name
melt_AccessoryTable$Genes <- rep(row.names(AccessoryTable), 8)

#To see the first rows of the table
head(AccessoryTable)

# Loading the ggplot2 to plot the heatmap
library(ggplot2)

ggplot(melt_AccessoryTable, aes(variable, Genes)) +
  geom_tile(aes(fill = value), colour = "white") +
  scale_fill_gradient(low = "white", high = "red")
theme(axis.text.y = element_text(angle = 90))

  
  

  

  
 


  
        