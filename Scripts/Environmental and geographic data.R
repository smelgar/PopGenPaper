#https://gis.stackexchange.com/questions/227585/using-r-to-extract-data-from-worldclim
library(raster)
library(sp)
r <- getData("worldclim",var="bio",res=2.5)
#https://community.rstudio.com/t/raster-getdata-function-alternative/155898/3
#library(geodata)
#climVarGT <- worldclim_country(country= "GT", var="bio", res=10, path= "./Output", version = "2.1")
#climVarSV <- worldclim_country(country= "SV", var="bio", res=10, path= "./Output", version = "2.1")

#Bio 1 and Bio12 are mean anual temperature and anual precipitation:
r <- r[[c(1,12)]]
names(r) <- c("Temp","Prec")
#create random points as example, in your case use coordinates to create a SpatialPoint object.

coords <- data.frame(row.names = c("Azacualpa","Cerro Alto ","Chilcuyo",
                                   "El Cerron","El Jute","El Zacatal","El_Carrizal","El_Chaperno",
                                   "La Prensa","La Primavera","Llano de Dona Maria",
                                   "Llano Santa María"),
           x = c(-89.52917,-89.65,-89.5279,-89.28631,
                           -89.64132,-89.48223,-89.98881,-89.7921,-89.26602,-89.5354,
                           -89.82592,-89.81717),
           y = c(13.75968,13.68323,14.06611,14.70905,
           14.11384,13.89918,14.37851,14.35872,14.72538,13.95435,
           13.94329,14.26666)
)



points <- SpatialPoints(coords, proj4string = r@crs)
#points <- spsample(as(r@extent, 'SpatialPolygons'),n=10, type="random")    
#Finally, use extract. With cbind.data.frame and coordinates you will get the desire data.frame.

values <- extract(r,points)

df <- cbind.data.frame(coordinates(points),values)
rownames(df) <- rownames(coords)
#I used random points, so I got a lot of NA. It is to be expected.

head(df)
#temperature
plot(r[[1]],ext=c(-90.5,-89,13.6,14.8))
plot(points,add=T)
#annual pp
plot(r[[2]],ext=c(-90.5,-89,13.6,14.8))
plot(points,add=T)

#https://stackoverflow.com/questions/32363998/function-to-calculate-geospatial-distance-between-two-points-lat-long-using-r
library(geosphere)
distances <- distm(coords, fun = distHaversine) #distances in meters
#add location names to rows and columns
rownames(distances) <- rownames(coords) 
colnames(distances) <- rownames(coords)



#rownames(FstVillages)=colnames(FstVillages)=c("El_Chaperno",	"Chilcuyo",	"Llano de Dona Maria",	"La Primavera",	"Cerro Alto", 	"Azacualpa",	"El Zacatal",	"La Prensa",	"El Jute",	"El_Carrizal")

geodist <- data.frame(
  El_Chaperno = c(0,43288.04,46389.32,
                  52858.37,72477.52,61121.83,69852.41,31745.75,
                  21326.77),
  Chilcuyo = c(43288.04,0,34972.6,
               12467.41,34111.91,19226.26,78634.29,13348.99,
               60687.82),
  Llano.de.Dona.Maria = c(46389.32,34972.6,0,
                          31411,38032.49,37458.84,105954.39,27530.38,51540.01),
  La.Primavera = c(52858.375,12467.411,
                   31410.999,0,21681.024,8409.606,90614.359,21120.318,
                   68003.134),
  Azacualpa = c(72477.52,34111.91,
                38032.49,21681.02,0,16336.95,111187.56,41244.93,
                84904.92),
  El.Zacatal = c(61121.832,19226.264,
                 37458.842,8409.606,16336.954,0,94882.843,29432.562,
                 76403.32),
  La.Prensa = c(69852.41,78634.29,
                105954.39,90614.36,111187.56,94882.84,0,79193.11,
                86926.27),
  El.Jute = c(31745.75,13348.99,
              27530.38,21120.32,41244.93,29432.56,79193.11,0,
              47684.09),
  El_Carrizal = c(21326.77,60687.82,
                  51540.01,68003.13,84904.92,76403.32,86926.27,47684.09,
                  0)
)

FstDist <- data.frame(
  El_Chaperno = c(0,0.0361,0.2188,0.0876,
                  0.0773,0.0847,0.1171,0.021,0.0798),
  Chilcuyo = c(0.0361,0,0.2261,0.0519,
               0.0059,0.0161,0.0736,-0.004,0.067),
  Llano.de.Dona.Maria = c(0.2188,0.2261,0,0.2002,
                          0.1967,0.2131,0.351,0.1955,0.2553),
  La.Primavera = c(0.0876,0.0519,0.2002,0,
                   0.0394,0.0743,0.1505,0.0242,0.1104),
  Azacualpa = c(0.0773,0.0059,0.1967,
                0.0394,0,0.0644,0.1154,0.0326,0.1006),
  El.Zacatal = c(0.0847,0.0161,0.2131,
                 0.0743,0.0644,0,0.132,0.0443,0.1174),
  La.Prensa = c(0.1171,0.0736,0.351,
                0.1505,0.1154,0.132,0,0.0738,0.0969),
  El.Jute = c(0.021,-0.004,0.1955,
              0.0242,0.0326,0.0443,0.0738,0,0.0572),
  El_Carrizal = c(0.0798,0.067,0.2553,
                  0.1104,0.1006,0.1174,0.0969,0.0572,0)
)

rownames(FstDist) <- colnames(FstDist)
rownames(geodist) <- colnames(geodist)
FstDist = FstDist[-8,-8] # erase "El Jute" wich has Fst= -0.04 from Chilcuyo.
geodist = geodist[-8,-8]
FstDist <- as.dist(FstDist)
geodist <- as.dist(geodist)
#https://stats.oarc.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/
library(ade4)
mantel.rtest(FstDist, geodist, nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(FstDist,geodist,xlab="Fst",ylab="meters",main="Genetic vs Geographic distances")

# distance matrix with climatic data
distclim <- dist(df[,c(3,4)])

# Clúster jerárquico https://r-charts.com/es/parte-todo/hclust/
hc <- hclust(distclim)

# Dendrograma
plot(hc)

#Mantel test for all variables against Fst distances
#first find projection of all coordinates in order to use tem for calculating distances
#######UTM projection https://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm
#library(sp)
library(rgdal)

#Function
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}


UTM = LongLatToUTM(coords$x,coords$y,16)
rownames(UTM) = rownames(coords)
#next define dfUTM with projection data instead of coordinates.
dfUTM = df[,c(3,4)]
dfUTM$X = UTM$X
dfUTM$Y = UTM$Y
#next standardize all variables by column https://www.marsja.se/how-to-standardize-data-in-r-numeric-only/#:~:text=or%20other%20factors.-,What%20does%20it%20mean%20to%20standardize%20variables%20in%20R%3F,%E2%80%9Cstandardization%20to%20unit%20variance%E2%80%9D.
dfUTM_std = scale(dfUTM)
dfUTM_std = dfUTM_std[c(8,3,11,10,1,6,9,7),]
#to compare dfUTM_std villages with FstDist villages:
#rownames(as.matrix(FstDist))
#rownames(dfUTM_std[c(8,3,12,10,1,6,9,7),])
#now find a distance matrix using geographic and environmental variables
distallvars = dist(dfUTM_std) #this uses euclidean distances
#next do mantel test
mantel.rtest(FstDist, distallvars, nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(FstDist,distallvars,xlab="Fst",ylab="geographic and environmental distances")

# Altitude and genetic distances:
dfUTM$Altitude <- c(2011L,963L,1521L,4060L,2296L,2650L,
                         5053L,3344L,3048L,2396L,2316L,3499L)
#Mantel test with just altitude and Fst
mantel.rtest(FstDist, dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(FstDist,dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]),xlab="Fst",
     ylab="Altitude differences")

# altitude and temperature
mantel.rtest(dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]), dist(dfUTM$Temp[c(8,3,11,10,1,6,9,7)]), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]),dist(dfUTM$Temp[c(8,3,11,10,1,6,9,7)]),
     xlab="Altitude",ylab="Temperature differences")

#altitude and precipitation
mantel.rtest(dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]), dist(dfUTM$Prec[c(8,3,11,10,1,6,9,7)]), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]),
     dist(dfUTM$Prec[c(8,3,11,10,1,6,9,7)]),xlab="Altitude",ylab="Precipitation differences")

# this is a data frame taken from genind analysis using Hierftat
VillageHoHsFis <- data.frame(
  stringsAsFactors = FALSE,
  Village = c("Guatemala_Jutiapa_LlanoSantaMaría","Guatemala_Jutiapa_ElChaperno",
              "ElSalvador_SonsonateAhuachapan_","ElSalvador_SantaAna_Chilguyo",
              "ElSalvador_SonsonateAhuachapan_CaserioSanRaymundo",
              "ElSalvador_SantaAna_CaserioSaritaLoma",
              "ElSalvador_SonsonateAhuachapan_CantónLasFlores",
              "ElSalvador_SonsonateAhuachapan_ElAguacate","ElSalvador_SantaAna_ElZacatal",
              "ElSalvador_SantaAna_SantaAna","Guatemala_Chiquimula_ElCerron",
              "Guatemala_Chiquimula_LaPrensa","ElSalvador_SantaAna_ElJute",
              "Guatemala_Jutiapa_ElCarrizal"),
  Hs = c("NaN","0.13081269","NaN",
         "0.12917018","0.09292035","0.13067553","0.13752914",
         "0.12478989","0.14075258","0.14318828","NaN","0.11404729",
         "0.14762532","0.12969333"),
  Ho = c(0.0813449,0.08078505,
         0.07373272,0.07544387,0.09301075,0.07435345,0.08655914,
         0.06934161,0.09801075,0.09497753,0.08948546,0.06291028,
         0.12096344,0.07975882),
  Fis = c("NaN","0.31220881","NaN",
          "0.2879435","-0.15602837","0.24861878","0.20994475",
          "0.35233414","0.22677875","0.22262357","NaN","0.28225806",
          "0.09175505","0.31221695")
)

dfUTM$Village2 <- rownames(VillageHsHoFis[c(8,7,4,11,13,9,14,2,12,6,5,1),])
dfUTM$Hs <- VillageHsHoFis$Hs[c(8,7,4,11,13,9,14,2,12,6,5,1)]
dfUTM$Ho <- VillageHsHoFis$Ho[c(8,7,4,11,13,9,14,2,12,6,5,1)]
dfUTM$Fis <- VillageHsHoFis$Fis[c(8,7,4,11,13,9,14,2,12,6,5,1)]

#Altitude and Ho
mantel.rtest(dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]), 
             dist(dfUTM$Ho[c(8,3,11,10,1,6,9,7)]), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]),dist(dfUTM$Ho[c(8,3,11,10,1,6,9,7)]),
     xlab="Altitude",ylab="Ho differences")

#Altitude and Hs
mantel.rtest(dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]), 
             dist(dfUTM$Hs[c(8,3,11,10,1,6,9,7)]), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]),dist(dfUTM$Hs[c(8,3,11,10,1,6,9,7)]),
     xlab="Altitude",ylab="Hs differences")

#Altitude and Fis
mantel.rtest(dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]), 
             dist(dfUTM$Fis[c(8,3,11,10,1,6,9,7)]), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(dist(dfUTM$Altitude[c(8,3,11,10,1,6,9,7)]),
     dist(dfUTM$Fis[c(8,3,11,10,1,6,9,7)]),xlab="Altitude",ylab="Fis differences")



# Temperature and genetic distances:

mantel.rtest(FstDist, dist(dfUTM$Temp[c(8,3,11,10,1,6,9,7)]), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(FstDist,dist(dfUTM$Temp[c(8,3,11,10,1,6,9,7)]),xlab="Fst",ylab="Temperature differences")

# Precipitation and genetic distances:

mantel.rtest(FstDist, dist(dfUTM$Prec[c(8,3,11,10,1,6,9,7)]), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(FstDist,dist(dfUTM$Prec[c(8,3,11,10,1,6,9,7)]),xlab="Fst",ylab="Precipitation differences")

# geographic and genetic distances:

mantel.rtest(FstDist, dist(dfUTM[c(8,3,11,10,1,6,9,7),3:4]), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(FstDist,dist(dfUTM[c(8,3,11,10,1,6,9,7),3:4]),xlab="Fst",ylab="Geographic difference")

#Mantel test with all variables including altitude:
dfUTM_std = scale(dfUTM)
dfUTM_std = dfUTM_std[c(8,3,11,10,1,6,9,7),]
#to compare dfUTM_std villages with FstDist villages:
#rownames(as.matrix(FstDist))
#rownames(dfUTM_std[c(8,3,12,10,1,6,9,7),])
#now find a distance matrix using geographic and environmental variables
distallvars = dist(dfUTM_std) #this uses euclidean distances
#next do mantel test
mantel.rtest(FstDist, distallvars, nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(FstDist,distallvars,xlab="Fst",ylab="geographic and environmental distances")

#Mantel test for individuales
load("GeneticDistance47.Robj")
indivdf <- data.frame(
  stringsAsFactors = FALSE,
  sample = c("A9944","CHJ0007","CHJ0011",
             "CHJ0083","CHJ0155","CHJ0188","S236","S327","S328",
             "S333","S336","S358","S359","S443","S445","S446",
             "S447","S448","S451","S455","S476","S478","S486",
             "S507a","S507b","TPG1017","TPG154","TPG715","TPS0053",
             "TPS0068","TPS0180","TPS0225","TPS0330","CHJ1016",
             "CHJ0168","CHJ0169","CHJ0172","CHJ0360","CHJ0362",
             "CHJ0363","CHJ0364","CHJ0366","CHJ0367","CHJ0376",
             "CHJ0378","CHJ0386"),
  Village = c("Llano Santa María",
              "El Chaperno","El Chaperno","El Chaperno","El Chaperno",
              "El Chaperno","Chilcuyo","Llano de Dona Maria",
              "Llano de Dona Maria","La Primavera","La Primavera","Cerro Alto",
              "Cerro Alto","Azacualpa","Azacualpa","Azacualpa",
              "Azacualpa","Azacualpa","Azacualpa","Azacualpa",
              "El Zacatal","El Zacatal","El Zacatal","El Zacatal","El Zacatal",
              "El Cerron","La Prensa","La Prensa","Chilcuyo",
              "Chilcuyo","El Jute","El Jute","El Jute","El Carrizal",
              "El Chaperno","El Chaperno","El Chaperno","El Carrizal",
              "El Carrizal","El Carrizal","El Carrizal",
              "El Carrizal","El Carrizal","El Carrizal","El Carrizal",
              "El Carrizal"),
  Temp = c(22.5,22.6,22.6,22.6,22.6,
           22.6,24,23.6,23.6,22.9,22.9,24.4,24.4,23,23,23,
           23,23,23,23,22.4,22.4,22.4,22.4,22.4,22.2,22.2,
           22.2,24,24,23.3,23.3,23.3,17.8,22.6,22.6,22.6,
           17.8,17.8,17.8,17.8,17.8,17.8,17.8,17.8,17.8),
  Prec = c(134.4,118.1,118.1,118.1,
           118.1,118.1,173.7,175.3,175.3,183.8,183.8,198.6,
           198.6,181.3,181.3,181.3,181.3,181.3,181.3,181.3,184.5,
           184.5,184.5,184.5,184.5,129,129,129,173.7,173.7,
           157.7,157.7,157.7,169,118.1,118.1,118.1,169,169,
           169,169,169,169,169,169,169),
  X = c(19601.73,19884.73,19884.73,
        19884.73,19884.73,19884.73,22700.97,19464.18,19464.18,
        22606.67,22606.67,21334.54,21334.54,22651.22,
        22651.22,22651.22,22651.22,22651.22,22651.22,22651.22,
        23175.09,23175.09,23175.09,23175.09,23175.09,25381.81,
        25602.2,25602.2,22700.97,22700.97,21481.32,21481.32,
        21481.32,17764.23,19884.73,19884.73,19884.73,
        17764.23,17764.23,17764.23,17764.23,17764.23,17764.23,
        17764.23,17764.23,17764.23),
  Y = c(157906.1,158922,158922,
        158922,158922,158922,155650.2,154327,154327,154414,
        154414,151426.3,151426.3,152258.5,152258.5,152258.5,
        152258.5,152258.5,152258.5,152258.5,153797.3,153797.3,
        153797.3,153797.3,153797.3,162739.4,162917.9,
        162917.9,155650.2,155650.2,156192,156192,156192,159167.8,
        158922,158922,158922,159167.8,159167.8,159167.8,
        159167.8,159167.8,159167.8,159167.8,159167.8,159167.8),
  Elevation = c(349.9,334.4,334.4,334.4,
                334.4,334.4,152.1,231.6,231.6,239.6,239.6,96.3,96.3,
                201.1,201.1,201.1,201.1,201.1,201.1,201.1,265,
                265,265,265,265,406,304.8,304.8,152.1,152.1,229.6,
                229.6,229.6,505.3,334.4,334.4,334.4,505.3,505.3,
                505.3,505.3,505.3,505.3,505.3,505.3,505.3)
)

mantel.rtest(GeneticDistance47, dist(indivdf$Elevation), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(GeneticDistance47, dist(indivdf$Elevation),xlab="Genetic distance between individuals",ylab="elevation differences")

mantel.rtest(GeneticDistance47, dist(indivdf$Temp), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(GeneticDistance47, dist(indivdf$Temp),xlab="Genetic distance between individuals",ylab="temperature differences")

mantel.rtest(GeneticDistance47, dist(indivdf$Prec), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(GeneticDistance47, dist(indivdf$Prec),xlab="Genetic distance between individuals",ylab="precipitation differences")

mantel.rtest(GeneticDistance47, dist(indivdf[,5:6]), nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(GeneticDistance47, dist(indivdf[,5:6]),xlab="Genetic distance between individuals",ylab="geographic distances")

#####
#next standardize all variables by column https://www.marsja.se/how-to-standardize-data-in-r-numeric-only/#:~:text=or%20other%20factors.-,What%20does%20it%20mean%20to%20standardize%20variables%20in%20R%3F,%E2%80%9Cstandardization%20to%20unit%20variance%E2%80%9D.
indivdf_std = scale(indivdf[3:7])

#now find a distance matrix using geographic and environmental variables for individuals
distallvarsind = dist(indivdf_std) #this uses euclidean distances
#next do mantel test
mantel.rtest(GeneticDistance47, distallvarsind, nrepet = 9999)
par(mar = c(4, 4, 4, 4))
plot(GeneticDistance47,distallvarsind,xlab="Fst",ylab="geographic and environmental distances",main="Pairs of individuals")


# BIO1 = Annual Mean Temperature
# 
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# 
# BIO3 = Isothermality (BIO2/BIO7) (×100)
# 
# BIO4 = Temperature Seasonality (standard deviation ×100)
# 
# BIO5 = Max Temperature of Warmest Month
# 
# BIO6 = Min Temperature of Coldest Month
# 
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# 
# BIO8 = Mean Temperature of Wettest Quarter
# 
# BIO9 = Mean Temperature of Driest Quarter
# 
# BIO10 = Mean Temperature of Warmest Quarter
# 
# BIO11 = Mean Temperature of Coldest Quarter
# 
# BIO12 = Annual Precipitation
# 
# BIO13 = Precipitation of Wettest Month
# 
# BIO14 = Precipitation of Driest Month
# 
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# 
# BIO16 = Precipitation of Wettest Quarter
# 
# BIO17 = Precipitation of Driest Quarter
# 
# BIO18 = Precipitation of Warmest Quarter
# 
# BIO19 = Precipitation of Coldest Quarter




