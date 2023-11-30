
library(spData)
library(sp)
library(sf)
library(spdep)
library("INLA")
library("spatstat")
library("sp")
library("maptools")
library("latticeExtra")
library("gridExtra")
library("spdep")
library("rgdal")
library("rnaturalearth")
library("wbstats")
library("leaflet")
library("DT")
library("ggplot2")
library("sp")
library("RColorBrewer")
library(ggplot2)
library(pander)
library(SpatialEpi)
library(dplyr)
library("INLA")
library("spatstat")
library("sp")
library("maptools")
library("latticeExtra")
library("gridExtra")
library("spdep")
library("rgdal")
library(spdep)

library(INLA)

# read data 
setwd(getwd())
md= read.csv("~/G6PD.csv", header = TRUE) 

# select few columns

md$number_males[is.na(md$number_males)]<-0

md$number_females[is.na(md$number_females)]<-0
md$number_males_deficient[is.na(md$number_males_deficient)]<-0
md$number_females_deficient[is.na(md$number_females_deficient)]<-0

md1 = data.frame(country = md$country)


# The total number of people
md1$population =md$number_females+md$number_males


# The total number of the cases
md1$cases = md$number_females_deficient+md$number_males_deficient



md2<-md1
# aggregate data
md3 <- aggregate(
  x = md2$cases  ,
  by = list(country = md2$country),
  FUN = sum
)


md4 <- aggregate(
  x = md2$population  ,
  by = list(country = md2$country),
  FUN = sum
)

# the expected value
md3$pop <- md4$x
overall_rate = sum(md3$x)/sum(md3$pop)

md3$E <- overall_rate*md3$pop
#md4 = subset(md3, md3$E!=0)
head(md3)


# Malaria disease data

# read data 
mdd= read.csv("~/public_pf.csv", header = TRUE) 

mdd$examined[is.na(mdd$examined)]<-0
mdd$pf_pos[is.na(mdd$pf_pos)]<-0


mdd1 = data.frame(country = mdd$country,population2 = mdd$examined,cases2 = mdd$pf_pos)

mdd2<-mdd1

# aggregate data
mdd3 <- aggregate(
  x = mdd2$cases2  ,
  by = list(country = mdd2$country),
  FUN = sum
)

mdd4 <- aggregate(
  x = mdd2$population2  ,
  by = list(country = mdd2$country),
  FUN = sum
)

# the expected value
mdd3$pop <- mdd4$x

overall_rate = sum(mdd3$x)/sum(mdd3$pop)


mdd3$E <- overall_rate*mdd3$pop

head(mdd3)

# Merge data 

mddd3 = merge(mdd3, md3, by = "country")

mddd3 = mddd3[order(mddd3$country),]


mddd3 = mddd3[c(-1,-3, -4,-5,-7,-8,-9,-11,-14,-19,-20,-21,-25,-27,-29,-30,-32,
                -33,-34,-35,-36,-38,-39,-41,-42,-44,-46,-47),] 
Common.countries = mddd3$country
length(Common.countries)

mddd3 = mddd3[order(mddd3$country),]

merg.data = data.frame(country = mddd3$country , M.cases = mddd3$x.x ,M.pop = mddd3$pop.x ,M.E = mddd3$E.x,
                       D.cases = mddd3$x.y , D.pop = mddd3$pop.y , D.E = mddd3$E.y ,
                       b1 = 1:length(mddd3$x.y), b2 = 1:length(mddd3$x.y), ID = 1:length(mddd3$x.y))



merg.data$SMRm = merg.data$M.cases / merg.data$M.E

merg.data$SMRg = merg.data$D.cases / merg.data$D.E
merg.data <- merg.data[order(merg.data$country),]

# The map 
map <- ne_countries(returnclass = "sp")
names(map)[names(map) == "iso_a3"] <- "ISO3"
names(map)[names(map) == "name"] <- "NAME"
map <- map[order(map$name_long),]
map$name_long[157] <- "Gambia"
map$name_long[131] <- "Congo"

map <- map[order(map$name_long),]

# Select only the common countries from the map

sub_map <- match(map$name_long , Common.countries)
map$cc <- sub_map
map_2 <- subset(map,map$cc!="NA" )
map_2 <- map_2[order(map_2$name_long),]

Our_index <- match(map$name_long , Common.countries)

map_2$D.E = merg.data$D.E
map_2$M.E = merg.data$M.E
map_2$M.cases = merg.data$M.cases
map_2$D.cases = merg.data$D.cases
map_2$M.pop = merg.data$M.pop
map_2$D.pop = merg.data$D.pop


library(cleangeo)
rr <- clgeo_CollectionReport(map_2)
summary(rr)
issues <- rr[rr$type == NA,]
map_2_c <- clgeo_Clean(map_2)


