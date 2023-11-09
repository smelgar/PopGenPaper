#Relatedness analysis 20230324
#if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("gdsfmt")
# BiocManager::install("SNPRelate")
library(SNPRelate)
library(here)
library(dendextend)
library(circlize)


#Here the analysis was done with VCF filtered for missingness and only 50 individuals, MAF=0.01

setwd("/home/sergio/OneDrive/DocumentosData 1/Proyectos/Proyecto Salvador")
vcffile =  "batch_1NombreCorto_(1)_filteredMAF0.01.vcf"
gdsfile = "snpsNCMAF0.01.gds"
snpgdsVCF2GDS(vcffile, gdsfile, verbose=FALSE)
genofile <- snpgdsOpen(gdsfile)
#when I repeated this, I needed to use showfile.gds(closeall=TRUE) to close the gbs file, from https://github.com/zhengxwen/gdsfmt/issues/17
#to see the sample ids: read.gdsn(index.gdsn(genofile, "sample.id")) taken from https://www.bioconductor.org/packages/devel/bioc/vignettes/SNPRelate/inst/doc/SNPRelate.html
RV <- snpgdsIBS(genofile)
colors <- read.csv("colors.csv")
m <- 1 - RV$ibs
#colnames(m) <- rownames(m) <- RV$sample.id
colnames(m) <- rownames(m) <- colors$SampleID
#select only the samples with low missingness taken with match(rownames(SNPSmsSnpInd),vcfsamplenames)-1
samples50 <- c("A9944", "CHJ0007", "CHJ0011", "CHJ0083", "CHJ0155", "CHJ0188", "S174", "S236", "S327",
               "S328", "S333", "S336", "S358", "S359", "S443", "S445", "S446", "S447",
               "S448", "S451", "S455", "S476", "S478", "S486", "S507a", "S507b", "S649",
               "TEX0010", "TEX0023", "TPG1017", "TPG154", "TPG715", "TPS0053", "TPS0068", "TPS0180", "TPS0225",
               "TPS0330", "CHJ1016", "CHJ0168", "CHJ0169", "CHJ0172", "CHJ0360", "CHJ0362", "CHJ0363", "CHJ0364",
               "CHJ0366", "CHJ0367", "CHJ0376", "CHJ0378", "CHJ0386")
match(samples50,rownames(m))

GeneticDistance <- as.dist(m[match(samples50,rownames(m)),match(samples50,rownames(m))]) #these numbers have to be changed to adjust to the number of samples
HC <- hclust(GeneticDistance, "ave")
#labels(GeneticDistance) #to see the labels
#graphic HC with font size 0.6
HCD <- as.dendrogram(HC)
#exclude samples from colors object
colors <- colors[match(samples50,rownames(m)),] #these numbers have to be changed to adjust to the number of samples

labels_colors(HCD) <- "black"#colors$HouseColor[order.dendrogram(HCD)]
  


plot(HC, cex = 0.8)
par(cex=0.7)
plot(HCD, horiz = TRUE)
#house
#par(mfrow=c(2,2))
par(xpd=TRUE)
par(cex=0.4)
par(mar = c(4,4,1,18),xpd=TRUE)
HCD <- set(HCD, "labels_cex", 0.5)
plot(HCD, horiz =TRUE)
colored_bars(cbind(colors[,c(10,8,6,4)]), HCD, rowLabels = c("Country","Department","Village","House"),horiz = TRUE)
legend("topright", inset = c(-0.2,0), fill = c("green","blue"), legend = levels(factor(colors$Country)),cex = 0.5, inset = c(-0.63,0), title = "Country", xjust = 0)
legend("topright", fill = unique(colors$DepartmentColor), legend = unique(colors$Department),cex = 0.5, inset=c(-0.64,0.11), title = "Department", xjust = 0)
legend("topright", fill = unique(colors$VillageColor)[-3], legend = unique(colors$Village)[-3], cex = 0.5, inset=c(-1.0,0.3), xpd=TRUE, title = "Village", xjust = 0)
rect.dendrogram(HCD,h=0.07,horiz=TRUE,which = 59)
rect.dendrogram(HCD,h=0.08,horiz=TRUE,which = 28)
rect.dendrogram(HCD,h=0.1,horiz=TRUE,which = 20)
#colores: https://stackoverflow.com/questions/18802519/label-and-color-leaf-dendrogram
#y https://stackoverflow.com/questions/34539746/color-side-bar-dendrogram-plot/
#legend outside plot: https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics

#https://r-charts.com/es/parte-todo/dendrograma-circular/
#color branches and other formats https://cran.r-project.org/web/packages/dendextend/vignettes/dendextend.html#the-set-function

#House with village.
par(mar=c(3,0,3,10),xpd=TRUE)
colors$HouseColor[colors$HouseColor=="white"] <- 'gray70'
  labels_colors(HCD) <- colors$HouseColor[order.dendrogram(HCD)]
  HCD <- set(HCD, "labels_cex", 0.7)
  HCD <- set(HCD, "branches_lwd", 3)
  #HCD <- set(HCD, "branches_col", c(rep(1,33),"aquamarine4",2,1,"aquamarine4",rep(1,100)))
  HCD %>% set("branches_col", c(1,1,1,1,1,"aquamarine",1,1,rep(1,26),"black","aquamarine4","aquamarine4",
                                rep(1,3),rep("chocolate2",2),1,1,"chocolate2",rep(1,4),
                                "green","green",1,"chartreuse4","chartreuse",1,1,1,1,1,1,1,1,
                                "chartreuse2","darkolivegreen2",
                                rep(1,27),"aquamarine",rep(1,19),"aquamarine",
                                1,1,1,1,1,1,1,1,1,"antiquewhite4","antiquewhite3",rep(1,100))) %>% 
    set("leaves_pch", c(rep(19,2),8,
                        rep(19,13),17,17,
                        19,rep(17,3),19,17,17,
                        9,9,rep(19,3),
                        9,9,rep(19,13),8,
                        rep(19,9),8,rep(19,6))) %>% 
    set("leaves_col", colors$VillageColor[order.dendrogram(HCD)]) %>%
    circlize_dendrogram()  
  #plot(horiz=TRUE)
  legend("topright", fill = unique(colors$VillageColor)[-3], legend = unique(colors$Village)[-3],cex = 0.5, 
         inset=c(-0.5,0), 
         title = "Village", xjust = 0)
  legend("topright",  legend = c("Same house clustered","Same village","Same house dispersed"),
         cex = 0.5, 
         inset=c(-0.24,0.8), 
         title = "Symbols", xjust = 0,pch = c(17, 9, 8))
  
  #Village, Department, Country
  par(mar=c(3,0,3,10),expand=TRUE)
  labels_colors(HCD) <- colors$CountryColor[order.dendrogram(HCD)]
  HCD %>% set("branches_col", c(1,"gray70",1,1,"darkorange","cornsilk3",1,"deeppink",
                                rep(1,2),"bisque",1,"brown3",1,"bisque","cornsilk3",
                                1,"coral1",1,1,"deeppink","cornsilk3",1,"burlywood",1,1,
                                "burlywood",1,"burlywood","cornsilk3",1,
                                "black",1,1,1,"cornsilk3","cornsilk3",
                                1,"burlywood",1,"cyan3","cyan3",1,
                                1,"cyan3",1,1,"cornflowerblue",1,
                                "cornflowerblue",
                                "cornflowerblue",1,"cornflowerblue",
                                "cornflowerblue",1,"deeppink",1,1,
                                "cornflowerblue",1,"cornflowerblue",
                                1,"cornflowerblue","cornflowerblue",1,
                                "brown1",1,1,1,"cornsilk3",1,"bisque2",
                                1,1,"orange","orange",1,1,"brown1","brown1",
                                1,"coral1",1,1,"bisque2",1,rep("bisque2",2),
                                1,1,"brown1","cornsilk3",1,"darkorange",
                                1,"brown1",1,
                                "brown1",1,"brown1",1,"bisque2",1,"gray",
                                1,"black",1,"bisque",1,1,"coral1","cornsilk3",
                                1,1,"coral1","coral1",1,"bisque2",1,1,1,"aquamarine",
                                "bisque",
                                rep(1,100))) %>% 
    set("leaves_pch", 19) %>% 
    set("leaves_col", colors$DepartmentColor[order.dendrogram(HCD)]) %>%
    
    #plot(horiz = TRUE)
    circlize_dendrogram()
  legend("topright", fill = c("green","blue"), legend = levels(factor(colors$Country)),cex = 0.5, inset=c(-0.15,-0.15), title = "Country (labels)", xjust = 0)
  legend("topright", fill = unique(colors$DepartmentColor), legend = unique(colors$Department),cex = 0.5, inset=c(-0.18,0.06), title = "Department (leaves)", xjust = 0)
  legend("topright", fill = unique(colors$VillageColor)[-3], legend = unique(colors$Village)[-3], cex = 0.5, inset=c(-0.51,0.4), xpd=TRUE, title = "Village (branches)", xjust = 0)
  
  
  #HCD <- set(HCD, "branches_col", c(rep(1,4),2,rep(1,4),3))
  
  #house
  #dev.new(width=10, height=10)
  
  #par(mfrow=c(2,2))
  par(mar=c(4,4,4,4),expand=TRUE)
  
  labels_colors(HCD) <- colors$HouseColor[order.dendrogram(HCD)]
  #plot(HCD, horiz =TRUE)
  circlize_dendrogram(HCD)
  
  
  #village
  
  
  par(mar=c(4,4,4,4),expand=TRUE)
  labels_colors(HCD) <- colors$VillageColor[order.dendrogram(HCD)]
  HCD %>% set("branches_col", c(1,"gray70",1,1,"darkorange","cornsilk3",1,"deeppink",
                                rep(1,2),"bisque",1,"brown3",1,"bisque","cornsilk3",
                                1,"coral1",1,1,"deeppink","cornsilk3",1,"burlywood",1,1,
                                "burlywood",1,"burlywood","cornsilk3",1,
                                "black",1,1,1,"cornsilk3","cornsilk3",
                                1,"burlywood",1,"cyan3","cyan3",1,
                                1,"cyan3",1,1,"cornflowerblue",1,
                                "cornflowerblue",
                                "cornflowerblue",1,"cornflowerblue",
                                "cornflowerblue",1,"deeppink",1,1,
                                "cornflowerblue",1,"cornflowerblue",
                                1,"cornflowerblue","cornflowerblue",1,
                                "brown1",1,1,1,"cornsilk3",1,"bisque2",
                                1,1,"orange","orange",1,1,"brown1","brown1",
                                1,"coral1",1,1,"bisque2",1,rep("bisque2",2),
                                1,1,"brown1","cornsilk3",1,"darkorange",
                                1,"brown1",1,
                                "brown1",1,"brown1",1,"bisque2",1,"gray",
                                1,"black",1,"bisque",1,1,"coral1","cornsilk3",
                                1,1,"coral1","coral1",1,"bisque2",1,1,1,"aquamarine",
                                "bisque",
                                rep(1,100))) %>% 
    #plot(horiz = TRUE)
    circlize_dendrogram()
  legend("topright", fill = c("green","blue"), legend = levels(factor(colors$Country)),cex = 0.5, inset=c(-0.63,0), title = "Country", xjust = 0)
  legend("topright", fill = unique(colors$DepartmentColor), legend = unique(colors$Department),cex = 0.5, inset=c(-0.64,0.11), title = "Department", xjust = 0)
  legend("topright", fill = unique(colors$VillageColor)[-3], legend = unique(colors$Village)[-3], cex = 0.5, inset=c(-1.0,0.3), xpd=TRUE, title = "Village", xjust = 0)
  
  
  #Department
  par(mar=c(4,4,4,4),expand=TRUE)
  labels_colors(HCD) <- colors$DepartmentColor[order.dendrogram(HCD)]
  HCD %>% set("branches_col", c(1,"chocolate2",1,1,"yellow","chocolate2",1,"yellow",
                                rep(1,2),"red",1,"red",1,"red","chocolate2",
                                1,"yellow",1,1,"yellow","chocolate2",1,"yellow",1,1,
                                "yellow",1,"yellow","chocolate2",1,
                                "green",1,1,1,"chocolate2","chocolate2",
                                1,"yellow",1,"chocolate2","chocolate2",1,
                                1,"chocolate2",1,1,"chocolate2",1,
                                "chocolate2",
                                "chocolate2",1,"chocolate2",
                                "chocolate2",1,"yellow",1,1,
                                "chocolate2",1,"chocolate2",
                                1,"chocolate2","chocolate2",1,
                                "green",1,1,"yellow","chocolate2",1,"yellow",
                                1,1,"aquamarine","aquamarine",1,1,"green","green",
                                1,"yellow",1,1,"yellow",1,rep("yellow",2),
                                1,1,"green","chocolate2",1,"yellow",
                                1,"green",1,
                                "green",1,"green",1,"yellow",1,"green",
                                1,"green",1,"red",1,1,"yellow","chocolate2",
                                1,1,"yellow","yellow",1,"yellow",1,"green",1,"red",
                                "red",
                                rep(1,100))) %>% 
    #plot(horiz = TRUE)
    circlize_dendrogram()
  
  #Country
  par(mar=c(4,4,4,4),expand=TRUE)
  labels_colors(HCD) <- colors$CountryColor[order.dendrogram(HCD)]
  HCD %>% set("branches_col", c(1,"blue",1,1,"green","blue",1,"green",
                                rep(1,2),"blue",1,"blue",1,"blue","blue",
                                1,"green",1,1,"green","blue",1,"green",1,1,
                                "green",1,"green","blue",1,
                                "green",1,1,1,"blue","blue",
                                1,"green",1,"blue","blue",1,
                                1,"blue",1,1,"blue",1,
                                "blue",
                                "blue",1,"blue",
                                "blue",1,"green",1,1,
                                "blue",1,"blue",
                                1,"blue","blue",1,
                                "green",1,1,"green","blue",1,"green",
                                1,1,"green","green",1,1,"green","green",
                                1,"green",1,1,"green",1,rep("green",2),
                                1,1,"green","blue",1,"green",
                                1,"green",1,
                                "green",1,"green",1,"green",1,"green",
                                1,"green",1,"blue",1,1,"green","blue",
                                1,1,"green","green",1,"green",1,1,1,"blue",
                                "blue",
                                rep(1,100))) %>% 
    #plot(horiz = TRUE)
    circlize_dendrogram()
  
  # Patricias' proposal
  #Country and Dept. 
  #labels of the dendrogram as sample ID
  HC$labels <- colors$SampleID
  HCD <- as.dendrogram(HC)

  HCD <- set(HCD, "labels_cex", 0.9)
  HCD <- set(HCD, "branches_lwd", 3)
  colors$CountryColor2 = ifelse(colors$Country == "Guatemala", "darkblue", "brown3")
  colors$DepartmentColor2 = ifelse(colors$Department == "Jutiapa", "blue2",
                                   ifelse(colors$Department == "Chiquimula", "dodgerblue3",
                                          ifelse(colors$Department == "Sonsonate", "magenta",
                                                 ifelse(colors$Department == "Ahuachapan","orange4",
                                                        ifelse(colors$Department == "Santa_Ana","deeppink3",NA
                                                        )))))
  
  par(mar=c(4,2,4,4),expand=TRUE)
  labels_colors(HCD) <- colors$DepartmentColor2[order.dendrogram(HCD)]
 HCD %>% set("branches_col", c("brown3",rep("darkblue",3),rep("brown3",7),rep("darkblue",5),
                                rep("brown3",12),
                                rep("brown3",8),rep("darkblue",3),"brown3",rep("darkblue",22),
                                rep("brown3",3),rep("darkblue",7),rep("brown3",6),"darkblue",
                                rep("brown3",100))) %>%
    set("leaves_pch", 19) %>% 
    set("leaves_col", colors$CountryColor2[order.dendrogram(HCD)]) %>%
    #plot(horiz = TRUE)
    circlize_dendrogram()
  legend("topright", fill = unique(colors$CountryColor2), legend = rev(levels(factor(colors$Country))),cex = 1, inset=c(-0.2,0), title = "\n Country \n (branches and leaves)", xjust = 0)
  legend("topright", fill = unique(colors$DepartmentColor2)[c(5,1,4,3,2)], legend = unique(colors$Department)[c(5,1,4,3,2)],cex = 1, inset=c(-0.2,0.5), title = "\n Department \n (labels)", xjust = 0)
  

  
  #Villages
  #https://stackoverflow.com/questions/66923282/create-new-column-based-on-another-variable
  colors$VillageColor2 = ifelse(colors$Village == "Santa_Ana", "orange",
                         ifelse(colors$Village == "El Zacatal, Amatepec", "hotpink",
                         ifelse(colors$Village == "El_Jute", "orange3",
                         ifelse(colors$Village == "El Zacatal, Amatepec", "orange4",
                         ifelse(colors$Village == "El_Jute", "magenta",
                         ifelse(colors$Village == "El Aguacate, Azacualpa, Municipio de Armenia", "red",
                         ifelse(colors$Village == "Chilguyo, Texistepeque", "coral1",
                         ifelse(colors$Village == "Caserio Sarita Loma, Aldea La Primavera", "deeppink",
                         ifelse(colors$Village == "Caserio San Raymundo, Canton Llano", "coral3",
                         ifelse(colors$Village == "Cantón Las Flores Cerro Alto, Municipio Caluco ", "coral4",
                         ifelse(colors$Village == "Llano_Santa_María", "cyan3",
                         ifelse(colors$Village == "El_Chaperno", "green",
                         ifelse(colors$Village == "El_Cerrón,_Olopa","blue",
                         ifelse(colors$Village == "La_Prensa,_Olopa","blueviolet",
                         ifelse(colors$Village == "El_Carrizal","dodgerblue",
                                "gray62"
                                                            )))))))))))))))
  
  
  #change house colors:
  colors$HouseColor = ifelse(colors$House == 'JUT-D-FM', 'deeppink2',
                      ifelse(colors$House == 'CHAP-D-65', 'deepskyblue4',
                      ifelse(colors$House == 'CHAP-D-125', 'cyan',
                      ifelse(colors$House == 'CHAP-D-155', 'cyan4',
                      ifelse(colors$House == 'CHAP-D-114-115', 'dodgerblue',
                      ifelse(colors$House == 'CHAP-D-54', 'dodgerblue4',
                      ifelse(colors$House == 'CHAP-D-64', 'deepskyblue',
                      ifelse(colors$House == '', 'gray62',
                      ifelse(colors$House == 'CHI-D137A', 'green',
                      ifelse(colors$House == 'CHI-D204', 'green3',
                      ifelse(colors$House == 'CHI-D48', 'green4',
                      ifelse(colors$House == 'CHI-D13', 'darkolivegreen',
                      ifelse(colors$House == 'CARR-D-39', 'deeppink2',
                      ifelse(colors$House == 'CHAP-D-81', 'purple',
                      ifelse(colors$House == 'CHAP-D-156', 'darkturquoise',
                      ifelse(colors$House == 'CARR-D-63', 'deeppink',
                      ifelse(colors$House == 'CARR-D-58', 'deeppink4',
                      ifelse(colors$House == 'CARR-D-57', 'darkorchid',
                      ifelse(colors$House == 'CARR-D-95', 'darkorchid4',
                      ifelse(colors$House == 'CARR-D-93', 'firebrick1',
                      ifelse(colors$House == 'CARR-D-18A', 'firebrick4',
                      ifelse(colors$House == 'CHI-D19', 'royalblue4',
                      NA))))))))))))))))))))))
  
  colors$VillageIDs = ifelse(colors$Village == "Santa_Ana", "SAn",
                                ifelse(colors$Village == "El Zacatal, Amatepec", "EZa",
                                ifelse(colors$Village == "El_Jute", "EJu",
                                ifelse(colors$Village == "El Aguacate, Azacualpa, Municipio de Armenia", "EAg",
                                ifelse(colors$Village == "Chilguyo, Texistepeque", "ChT",
                                ifelse(colors$Village == "Caserio Sarita Loma, Aldea La Primavera", "CSL",
                                ifelse(colors$Village == "Caserio San Raymundo, Canton Llano", "CSR",
                                ifelse(colors$Village == "Cantón Las Flores Cerro Alto, Municipio Caluco ", "CLF",
                                ifelse(colors$Village == "Llano_Santa_María", "LSM",
                                ifelse(colors$Village == "El_Chaperno", "ECh",
                                ifelse(colors$Village == "El_Cerrón,_Olopa","ECe",
                                ifelse(colors$Village == "La_Prensa,_Olopa","LPr",
                                ifelse(colors$Village == "El_Carrizal","ECa",""
                                )))))))))))))
  #create new variable using "unknown" for villages from El Salvador and keeping
  #Guatemalan villages.
  colors$Village2 = colors$Village
  colors$Village2 = ifelse(colors$Country == "El_Salvador","Unknown",
                           colors$Village2)
  #
  #change the labels for HC with village ids
  HC$labels <- paste0(colors$SampleID," ", colors$VillageIDs)
  HCD <- as.dendrogram(HC)
  
  par(mar=c(5,0,4,13),xpd=TRUE)
  colors$HouseColor[colors$HouseColor=="white"] <- 'gray62'
    labels_colors(HCD) <- colors$HouseColor[order.dendrogram(HCD)]
    HCD <- set(HCD, "labels_cex", 1)
    HCD <- set(HCD, "branches_lwd", 4)
    #HCD <- set(HCD, "branches_col", c(rep(1,33),"aquamarine4",2,1,"aquamarine4",rep(1,100)))
    HCD %>% 
      set("branches_col", c(rep(1,7),1,rep(1,3),rep("purple",3),rep(1,12),1,1,
                            rep(1,10),rep(1,7),rep(1,4),
                            1,1,rep("firebrick1",3),rep(1,1),rep("firebrick4",4),
                            rep(1,100)
      )) %>%
      set("leaves_cex", 2) %>% 
      set("leaves_bg" , colors$VillageColor2[order.dendrogram(HCD)]) %>% 
      #in order to know the houses with more than 1 bug I used:
      #table(colors$House)
      #then I used:
      #colors[which(colors$House=="CARR-D-18A"),]
      #colors[which(colors$House=="CARR-D-93"),]
      #colors[which(colors$House=="CHAP-D-81"),]
      
      set("leaves_pch", c(rep(19,3),rep(13,2),rep(8,2),rep(19,3),
                          rep(13,4),rep(19,6),rep(13,3),
                          rep(8,6), rep(19,3),rep(13,4),rep(19,4),rep(13,6),
                          rep(19,4)
      )) %>%
      
      #set("leaves_pch", 19) %>% 
      set("leaves_col", colors$VillageColor2[order.dendrogram(HCD)]) %>%
      circlize_dendrogram()  
    #plot(horiz=TRUE)
    legend("topright",fill = unique(colors$VillageColor2)[c(7, 5, 6, 4, 8, 13, 9, 10, 14, 11, 2, 12, 1)], 
           legend = paste(unique(colors$Village)[c(7, 5, 6, 4, 8, 13, 9, 10, 14, 11, 2, 12, 1)],"(",unique(colors$VillageIDs)[c(7, 5, 6, 4, 8, 13, 9, 10, 14, 11, 2, 12, 1)],")"),cex = 0.8, 
           inset=c(-0.6,-0.2), 
           title = "\n Village \n (leaves)", xjust = 0)
    
    legend("topright", fill = unique(colors$HouseColor)[-7], legend = unique(colors$House)[-7],cex = 0.8, 
           inset=c(-0.4,0.39), 
           title = "\n House (branches \n and labels)", xjust = 0)
    
    
    legend("topright",  legend = c("Same house genetically similar","Same village genetically similar"),
           cex = 0.7, 
           inset=c(-0.1,1.005), 
           title = "Symbols", xjust = 0,pch = c(8, 13))
    
    # add table in legend: https://stackoverflow.com/questions/15406969/add-table-aligned-text-blocks-to-plot-in-r
    # with colors: https://stackoverflow.com/questions/61288600/how-to-put-a-table-below-a-ggplot-and-colour-the-rows-by-the-same-grouping-facto
    # https://stackoverflow.com/questions/52491085/plot-the-table-with-colorized-cells-in-r
    
    
    
    
    #exportar a un archivo, hay que transformar primero clase "dist" a "data.frame"
    #write.table(as.matrix(GeneticDistance), "/media/sergio/OS/Users/SergioAlejandro/Documents/Proyectos/Lenap 2015/Proyecto SNP/Parentesco dentro de casas/distanciasNC.txt", sep="\t")
    #save an svg image of the plot 
    #dev.copy(svg, "dendrogram20220823.svg")
    #create the plot and then:
    #dev.off()
    #save an png image of the plot
    #dev.copy(png, "dendrogram.png")
    #dev.off()
    

