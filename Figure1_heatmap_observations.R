library(ggplot2)
library(dplyr)
library(ggthemes)
library(viridis)
library(lubridate)
library(gdata)

path_data <- "path/to/data/observations.txt"
data <- read.csv2(path_data, sep = "\t")

# Delimit latitude and longitude of the plot
xmin <- -30
xmax <- 60
ymin <- -35
ymax <- -28

# Make sure variables are in the correct data type
data$Decimal_latitude <- as.numeric(data$Decimal_latitude)
data$Decimal_longitude <- as.numeric(data$Decimal_longitude)
data$Month <- as.numeric(data$Month)
data$Year <- as.numeric(data$Year)
data$Day <- as.numeric(data$Day)
str(data)

# Cut the data according to coordinates
data = data[which(data$Decimal_longitude >= xmin & data$Decimal_longitude <= xmax & data$Decimal_latitude >= ymin & data$Decimal_latitude <= ymax  ),]

# Identify the numeric column you want to use for the histogram
numeric_column <- "Month"  # Replace with the actual numeric column name

# Convert the column to numeric type if necessary
data[[numeric_column]] <- as.numeric(as.character(data[[numeric_column]]))

# Specify the number of breaks or breakpoints for the histogram bins
num_breaks <- 12  # Adjust the number of breaks as per your data

# Make an histogram to visualize
obs <- as.numeric(is.na(data[[numeric_column]]))
h <- hist(data[[numeric_column]], breaks = 13, col = "grey",
          xlab = "Month", main = "6")
lines(x = density(obs), col = "red")

# Perform counts for this delimited geographic area
data1 <- data %>% 
  count(Month) %>%
  mutate(Latitude = "South African Coast (Arid / Temperate)")

############
# Repeat process for the following geographic area
data <- read.csv2(path_data, sep = "\t")

xmin <- -30
xmax <- 60
ymin <- -28
ymax <- -11


data$Decimal_latitude <- as.numeric(data$Decimal_latitude)
data$Decimal_longitude <- as.numeric(data$Decimal_longitude)
data$Month <- as.numeric(data$Month)
data$Year <- as.numeric(data$Year)
data$Day <- as.numeric(data$Day)

data = data[which(data$Decimal_longitude >= xmin & data$Decimal_longitude <= xmax & data$Decimal_latitude >= ymin & data$Decimal_latitude <= ymax  ),]
numeric_column <- "Month"  # Replace with the actual numeric column name
data[[numeric_column]] <- as.numeric(as.character(data[[numeric_column]]))
num_breaks <- 10  # Adjust the number of breaks as per your data
h <- hist(data[[numeric_column]], breaks = 10, col = "grey", 
          xlab = "Month", main = "5")

data2 <- data %>% 
  count(Month) %>%
  mutate(Latitude = "Southern Africa (Arid / Savanna)")


############
# Repeat process for the following geographic area


data <- read.csv2(path_data, sep = "\t")

xmin <- -30
xmax <- 60
ymin <- -11
ymax <- 11

data$Decimal_latitude <- as.numeric(data$Decimal_latitude)
data$Decimal_longitude <- as.numeric(data$Decimal_longitude)
data$Month <- as.numeric(data$Month)
data$Year <- as.numeric(data$Year)
data$Day <- as.numeric(data$Day)

data = data[which(data$Decimal_longitude >= xmin & data$Decimal_longitude <= xmax & data$Decimal_latitude >= ymin & data$Decimal_latitude <= ymax  ),]
numeric_column <- "Month" 
data[[numeric_column]] <- as.numeric(as.character(data[[numeric_column]]))
num_breaks <- 10  
h <- hist(data[[numeric_column]], breaks = 10, col = "grey",
          xlab = "Month", main = "4")

data3 <- data %>% 
  count(Month) %>%
  mutate(Latitude = "Equatorial Africa (Tropical Humid)")


############
# Repeat process for the following geographic area


data <- read.csv2(path_data, sep = "\t")

xmin <- -30
xmax <- 60
ymin <- 11
ymax <- 21

data$Decimal_latitude <- as.numeric(data$Decimal_latitude)
data$Decimal_longitude <- as.numeric(data$Decimal_longitude)
data$Month <- as.numeric(data$Month)
data$Year <- as.numeric(data$Year)
data$Day <- as.numeric(data$Day)

data = data[which(data$Decimal_longitude >= xmin & data$Decimal_longitude <= xmax & data$Decimal_latitude >= ymin & data$Decimal_latitude <= ymax  ),]
numeric_column <- "Month"  
data[[numeric_column]] <- as.numeric(as.character(data[[numeric_column]]))
num_breaks <- 10  
h <- hist(data[[numeric_column]], breaks = 10, col = "grey",
          xlab = "Month", main = "3")

data4 <- data %>% 
  count(Month) %>%
  mutate(Latitude = "Sahel and Sahara (Desert / Arid)")


############
# Repeat process for the following geographic area

data <- read.csv2(path_data, sep = "\t")

xmin <- -30
xmax <- 60
ymin <- 21
ymax <- 45

data$Decimal_latitude <- as.numeric(data$Decimal_latitude)
data$Decimal_longitude <- as.numeric(data$Decimal_longitude)
data$Month <- as.numeric(data$Month)
data$Year <- as.numeric(data$Year)
data$Day <- as.numeric(data$Day)

data = data[which(data$Decimal_longitude >= xmin & data$Decimal_longitude <= xmax & data$Decimal_latitude >= ymin & data$Decimal_latitude <= ymax  ),]
numeric_column <- "Month" 
data[[numeric_column]] <- as.numeric(as.character(data[[numeric_column]]))
num_breaks <- 10  
h <- hist(data[[numeric_column]], breaks = 10, col = "grey",
          xlab = "Month", main = "2")

data5 <- data %>% 
  count(Month) %>%
  mutate(Latitude = "Southern Europe – Northern Africa (Temperate / Desert)")

############
# Repeat process for the following geographic area

data <- read.csv2(path_data, sep = "\t")

xmin <- -30
xmax <- 60
ymin <- 45
ymax <- 70

data$Decimal_latitude <- as.numeric(data$Decimal_latitude)
data$Decimal_longitude <- as.numeric(data$Decimal_longitude)
data$Month <- as.numeric(data$Month)
data$Year <- as.numeric(data$Year)
data$Day <- as.numeric(data$Day)

data = data[which(data$Decimal_longitude >= xmin & data$Decimal_longitude <= xmax & data$Decimal_latitude >= ymin & data$Decimal_latitude <= ymax  ),]
numeric_column <- "Month"  # Replace with the actual numeric column name
data[[numeric_column]] <- as.numeric(as.character(data[[numeric_column]]))
num_breaks <- 10  
h <- hist(data[[numeric_column]], breaks = 10, col = "grey",
          xlab = "Month", main = "1")

data6 <- data %>% 
  count(Month) %>%
  mutate(Latitude = "Northern Europe (Polar / Temperate)")


################ 
# CREATE HEATMAP
################

data_heatmap <- combine(data1, data2, data3, data4, data5, data6)
data_heatmap.nona <- na.omit(data_heatmap)


# Assign color variables
col1 = "#f6feee" #d8e1cf
col2 = "#438484"
data_heatmap.nona %>%
  count(n)

# create frequencies for each latitud
obs <- data_heatmap.nona %>%
  group_by(Latitude) %>%
  mutate(freq = n / sum(n))

# check levels if we have to reorder
obs$Latitude <- as.factor(obs$Latitude)
levels(obs$Latitude)


# reorder the levels
obs$Latitude <- factor(obs$Latitude, levels = levels(obs$Latitude)[c(4,5,1,3,6,2)])
obs$Month <- as.factor(obs$Month)

# do the heatmap plot
ggplot(data = obs, aes(x = Month, y = Latitude)) + 
  geom_tile(aes(fill = freq)) + scale_fill_gradient(low = col1, high = col2) + theme_tufte(base_family="Helvetica")
