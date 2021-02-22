#Names
  #Tristan Wennink



setwd("C:/Users/trist/Documents/Universiteit/Econometrics3")

data <- read.csv("data_assignment1.csv", sep = ";")
df1 <- as.data.frame(data)
attach(data)

x <- 1:length(gdp)

data$gdp <- as.numeric(data$gdp)
y <- as.numeric(gdp)

y

