##########################################################################
################CREATE PIE CHARTS OF GENETIC BELONGING MAP################

library(shapefiles)#Creates shapefile needed for map
library(mapplots)#creates map

#Pull Shapefiles
setwd("D:/2023/Paper/map/Guate")
Guate<-read.shapefile(shape.name = "GTM_adm1")
setwd("D:/2023/Paper/map/ElSal")
ElSal<-read.shapefile(shape.name = "SLV_adm1")
setwd("D:/2023/Paper/map/Hon")
Hon<-read.shapefile(shape.name = "HND_adm1")

#Create baseline map
xlim <- c(-90.5,-88.5)
ylim <- c(13.6,15)
basemap(xlim, ylim, main="Genetic Composition Map")
col=rainbow(3,alpha=0.6)
draw.shape(Guate, col="cornsilk")
draw.shape(ElSal, col="bisque")
draw.shape(Hon, col="darkgray")

#Add pie plots

#Carrizal
add.pie(z=c(9,90,1), x=-89.98881, y=14.37851, radius=0.12, col=col, labels=c("","",""), label.dist=1.2)
#text(-89.99,14.235,"ECa",cex=1.1) #Lets me add text at given coordinates

#Chaperno
add.pie(z=c(40,30,30), x=-89.79210, y=14.35872, radius=0.1, col=col, labels="")
#text(-89.79,14.237,"ECh",cex=1.1)

#LLano Santa Maria
add.pie(z=c(48,44,8), x=-89.81717, y=14.26666, radius=0.02, col=col, labels="")

#La Prensa
add.pie(z=c(46,49,5), x=-89.26602, y=14.72538, radius=0.04, col=col, labels=c("","",""),label.dist=2)
#text(-89.25,14.795,"LPr",cex=1.1)

#El Cerron
add.pie(z=c(42,57,1), x=-89.28631, y=14.69, radius=0.03, col=col, labels=c("","",""),label.dist=1.2)
#text(-89.32,14.65,"ECe",cex=1.1)

#Azacualpa
add.pie(z=c(81,17,2), x=-89.52917, y=13.75968, radius=0.09, col=col, labels="", label.dist=1.3)

#Zacatal
add.pie(z=c(94,4,2), x=-89.48223, y=13.89918, radius=0.07, col=col, labels="",cex=0.9,label.dist=1.7)
#text(-89.4,13.84,"F",cex=1.1) #Lets me add text at given coordinates

#Dona Maria
add.pie(z=c(74,25,1), x=-89.53540, y=13.95435, radius=0.04, col=col, labels=c("","",""),label.dist=1.2)

#Cerro Alto
add.pie(z=c(81,18,1), x=-89.65000, y=13.68323, radius=0.04, col=col, labels=c("","",""),label.dist=1.4)

#Chilcuyo
add.pie(z=c(82,15,3), x=-89.52790, y=14.06611, radius=0.05, col=col, labels="",label.dist=1.4)

#El Jute
add.pie(z=c(77,12,11), x=-89.64132, y=14.11384, radius=0.05, col=col, labels=c("","",""),label.dist=1.2)

#La Primavera
add.pie(z=c(98,1,1), x=-89.82592, y=13.94329, radius=0.04, col=col, labels="", label.dist=1.5)
, y=14.75, radius=0.15, col=col, labels=c("Group 1","Group 2","Group 3"),label.dist=1.2)


#########################################################################
###########CREATE AMERICA MAP############################################


library(shapefiles)#Creates shapefile needed for map
library(mapplots)#creates map

#Pull Shapefiles
setwd("D:/2023/Paper/map/Guate")
Guate<-read.shapefile(shape.name = "GTM_adm1")
setwd("D:/2023/Paper/map/ElSal")
ElSal<-read.shapefile(shape.name = "SLV_adm1")
setwd("D:/2023/Paper/map/Hon")
Hon<-read.shapefile(shape.name = "HND_adm1")
setwd("D:/2023/Paper/map/CA")
Mex<-read.shapefile(shape.name = "MEX_adm0")
setwd("D:/2023/Paper/map/Belize")
Bel<-read.shapefile(shape.name = "BLZ_adm0")
setwd("D:/2023/Paper/map/Nic")
Nic<-read.shapefile(shape.name = "NIC_adm0")


#Create baseline map
xlim <- c(-92,-85)
ylim <- c(12,22)
basemap(xlim, ylim, main="Genetic Composition Map")
draw.shape(Guate, col="cornsilk")
draw.shape(ElSal, col="bisque")
draw.shape(Hon, col="darkgray")
draw.shape(Mex, col="darkgray")
draw.shape(Bel, col="darkgray")
draw.shape(Nic, col="darkgray")
