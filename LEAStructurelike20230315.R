library(LEA)
setwd("/home/sergio/OneDrive/DocumentosData 1/Proyectos/Proyecto Salvador")
vcffile =  "batch_1NombreCorto_(1)_filteredMAF0.01ind50(copy).vcf"
genofile = "batch_1NombreCorto_(1)_filteredMAF0.01ind50(copy).geno"


######################
#borrar <- read.vcf("batch_1NombreCorto_(1)_filteredMAF0.01ind50.vcf")
#write.vcf(borrar,"batch_1NombreCorto_(1)_filteredMAF0.01ind50(copy).vcf")
######################

vcf2geno(vcffile, output.file = genofile, force = TRUE)

# cross entropy according to Lori's script
set.seed(1)
krange = 1:10
nrep = 300
project_Dept = snmf(genofile, K = krange, 
                    entropy = TRUE, repetitions = nrep, ploidy = 2, project = "new")
par(mar=(c(4,4,4,4)))

plot(project_Dept, col = "blue", pch = 18, cex = 2, 
     main = "Cross-entropy")
allbest = NULL
for (x in krange) {
  allbest[x] = min(cross.entropy(project_Dept,K=x))
}
 plot(allbest, col="blue", type="b", Main = "Best entropy values for each k", ylab = "entropy (min)",
      xlab = "k", pch = 18, cex = 2)
 allbest = NULL
 for (x in krange) {
   allbest[x] = median(cross.entropy(project_Dept,K=x))
 }
 plot(allbest, col="blue", type="b", Main = "Best entropy values for each k", ylab = "entropy (median)",
      xlab = "k", pch = 18, cex = 2)
ap <- which.min(allbest)
#select the best run for K = ap 
best = which.min(cross.entropy(project_Dept, K = ap))
#make a STRUCTURE diagram of the ancentry proportion for each individual
# barchart(project_Dept, K = ap, run = best, border = NA, space = 0, 
#          sort.by.Q = FALSE, col = (c("gray", "orange", "green","red")), xlab = "", 
#          ylab = "Ancestry proportions", main = c("Village", ap)) -> bp
# #label the axis with the breed or village
# axis(1, at = 1:length(bp$order), labels = Village_IDs$One_village_Breed, 
#      las=2, cex.axis = .5) 
##########################################

#from: https://bookdown.org/hhwagner1/LandGenCourse_book/WE_9.html
#snmf2 <- LEA::snmf(paste0(here::here(), "/data/stickleback.geno"), 
#                   K=1:8, ploidy=2, entropy=T, alpha=100, project="new")
# snmf2 <- LEA::snmf(genofile, K=1:8, ploidy=2, entropy=T, 
#                    alpha=100, project="new")
# par(mfrow=c(1,1))
# plot(snmf2, col="blue4", cex=1.4, pch=19)
# 
# K=5
# snmf = LEA::snmf(genofile, 
#                  K = K, alpha = 100, project = "new")

qmatrix = LEA::Q(project_Dept, K = ap, run= best)
#save(qmatrix,file="qmatix20230214c.Robj")
#50 samples ID:
#datapasta::vector_paste(rownames(SNPSmsSnpInd))
#50 samples ID
ids = c("A9944", "CHJ0007", "CHJ011", "CHJ083", "CHJ155", "CHJ188", "S174", "S236", 
        "S327", "S328", "S333", "S336", "S358", "S359", "S443", "S445", "S446", 
        "S447", "S448", "S451", "S455", "S476", "S478", "S486", "S507a", "S507b", 
        "S649", "TEX0010", "TEX0023", "TPG1017", "TPG154", "TPG715", "TPS0053", "TPS0068", 
        "TPS0180", "TPS0225", "TPS0330", "CHJ1016", "CHJ168", "CHJ169", "CHJ172", "CHJ360", 
        "CHJ362", "CHJ363", "CHJ364", "CHJ366", "CHJ367", "CHJ376", "CHJ378", "CHJ386")
##50 samples dept.
#Dept and country
#colors <- read.csv("colors.csv")
# Dept <- colors[match(ids,colors$SampleIDC),]$Department
# datapasta::vector_paste(Dept)
Dept <- c("Jutiapa", "Jutiapa", "Jutiapa", "Jutiapa", "Jutiapa", "Jutiapa", "Sonsonate", 
          "Santa_Ana", "Ahuachapan", "Ahuachapan", "Santa_Ana", "Santa_Ana", "Sonsonate", 
          "Sonsonate", "Sonsonate", "Sonsonate", "Sonsonate", "Sonsonate", "Sonsonate", 
          "Sonsonate", "Sonsonate", "Santa_Ana", "Santa_Ana", "Santa_Ana", "Santa_Ana", 
          "Santa_Ana", "Santa_Ana", "Santa_Ana", "Santa_Ana", "Chiquimula", "Chiquimula", 
          "Chiquimula", "Santa_Ana", "Santa_Ana", "Santa_Ana", "Santa_Ana", "Santa_Ana", 
          "Jutiapa", "Jutiapa", "Jutiapa", "Jutiapa", "Jutiapa", "Jutiapa", "Jutiapa", 
          "Jutiapa", "Jutiapa", "Jutiapa", "Jutiapa", "Jutiapa", "Jutiapa")
#to join Sonsonate and Ahuachapan:
Dept[Dept=="Ahuachapan"] <- "Sonsonate-Ahuachapan"
Dept[Dept =="Sonsonate"] <- "Sonsonate-Ahuachapan"

# (Country <- colors[match(ids,colors$SampleIDC),]$Country)
# datapasta::vector_paste(Country)
Country <- c("Guatemala", "Guatemala", "Guatemala", "Guatemala", "Guatemala", "Guatemala", 
             "El_Salvador", "El_Salvador", "El_Salvador", "El_Salvador", "El_Salvador", 
             "El_Salvador", "El_Salvador", "El_Salvador", "El_Salvador", "El_Salvador", 
             "El_Salvador", "El_Salvador", "El_Salvador", "El_Salvador", "El_Salvador", 
             "El_Salvador", "El_Salvador", "El_Salvador", "El_Salvador", "El_Salvador", 
             "El_Salvador", "El_Salvador", "El_Salvador", "Guatemala", "Guatemala", 
             "Guatemala", "El_Salvador", "El_Salvador", "El_Salvador", "El_Salvador", 
             "El_Salvador", "Guatemala", "Guatemala", "Guatemala", "Guatemala", "Guatemala", 
             "Guatemala", "Guatemala", "Guatemala", "Guatemala", "Guatemala", "Guatemala", 
             "Guatemala", "Guatemala")
#gsub use taken from: https://www.digitalocean.com/community/tutorials/sub-and-gsub-function-r
#Village = colors[match(ids,gsub('_.','',colors$SampleIDC)),]$Village
# datapasta::vector_paste(Village)
Village = c("Llano_Santa_María", "El_Chaperno", "El_Chaperno", "El_Chaperno", "El_Chaperno", 
            "El_Chaperno", NA, "Chilguyo, Texistepeque", "Caserio San Raymundo, Canton Llano", 
            "Caserio San Raymundo, Canton Llano", "Caserio Sarita Loma, Aldea La Primavera", 
            "Caserio Sarita Loma, Aldea La Primavera", 
            "Cantón Las Flores Cerro Alto, Municipio Caluco", 
            "Cantón Las Flores Cerro Alto, Municipio Caluco", 
            "El Aguacate, Azacualpa, Municipio de Armenia", 
            "El Aguacate, Azacualpa, Municipio de Armenia", 
            "El Aguacate, Azacualpa, Municipio de Armenia", 
            "El Aguacate, Azacualpa, Municipio de Armenia", 
            "El Aguacate, Azacualpa, Municipio de Armenia", 
            "El Aguacate, Azacualpa, Municipio de Armenia", 
            "El Aguacate, Azacualpa, Municipio de Armenia", 
            "El Zacatal, Amatepec", "El Zacatal, Amatepec", "El Zacatal, Amatepec", 
            "El Zacatal, Amatepec", "El Zacatal, Amatepec", "Santa_Ana", "Santa_Ana", 
            "Santa_Ana", "El_Cerrón,_Olopa", "La_Prensa,_Olopa", "La_Prensa,_Olopa", 
            "Chilguyo, Texistepeque", "Chilguyo, Texistepeque", "El_Jute", "El_Jute", 
            "El_Jute", "El_Carrizal", "El_Chaperno", "El_Chaperno", "El_Chaperno", 
            "El_Carrizal", "El_Carrizal", "El_Carrizal", "El_Carrizal", "El_Carrizal", "El_Carrizal", 
            "El_Carrizal", "El_Carrizal", "El_Carrizal")


rownames(qmatrix) = ids
qmatrix = as.data.frame(qmatrix)
qmatrix$Dept = Dept
qmatrix$Country = Country
qmatrix$Village = Village
#oder from https://stackoverflow.com/questions/1296646/sort-order-data-frame-rows-by-multiple-columns
qmatrix = qmatrix[with(qmatrix,order(Country,Dept,Village)),]



#vertical labels from:https://stackoverflow.com/questions/10286473/rotating-x-axis-labels-in-r-for-barplot
par(mar=c(10,8,0.5,0.5))
par(cex=1.5)
barplot(t(qmatrix[,1:ap]), col=RColorBrewer::brewer.pal(9,"Paired"), 
        #border=NA, 
        space=0, 
        #xlab="Individuals", 
        ylab="Admixture coefficients", las=2, cex.lab = 2 ) #ylabel size: cex.lab
#Add population labels to the axis:
#https://stackoverflow.com/questions/62408799/is-there-a-way-to-add-a-x-axis-below-the-primary-x-axis-in-r
#the following code is for text wrapping so text size can be increased
qmatrix$Dept2 = ifelse(qmatrix$Dept == "Santa_Ana", "Santa\nAna",
                ifelse(qmatrix$Dept == "Sonsonate-Ahuachapan", "Sonsonate-\nAhuachapan",
                       qmatrix$Dept
                 ))
for (i in 1:length(unique(qmatrix$Dept2))){
  axis(1, at=median(which(qmatrix$Dept2==unique(qmatrix$Dept2)[i]))-0.5, 
       labels=unique(qmatrix$Dept2)[i], line=6, 
       col="black", col.axis = "black", cex.lab = 3, 
       padj=0.5,cex.axis = 1.5)} #padj is to separate text from tick marks when I added text wrapping

#the following axis functions add the El Carrizal and El Chaperno labels
axis(1,at=median(which(qmatrix$Village=="El_Carrizal"))-0.5,labels="El Carrizal",
    line=3,col="black",col.axis="black",cex.lab=3,padj=0.5,cex.axis=1,lwd.tick=0)
axis(1,at=median(which(qmatrix$Village=="El_Chaperno"))-0.5,labels="El Chaperno",
     line=3,col="black",col.axis="black",cex.lab=3,padj=0.5,cex.axis=1,lwd.tick=0)

mtext("Individuals", 1, line=1, at=-4)
mtext("Villages",1,line=5,at=-5)
mtext("Dept", 1, line=8, at=-6)
# the following abline sometimes draws long lines, so I changed it to "segments()" https://r-charts.com/es/r-base/segments/
# for (i in 1:length(unique(qmatrix$Dept))){
#   abline(v = max(which(qmatrix$Dept==unique(qmatrix$Dept)[i])), col="red",lwd=4)}
# abline(v = 0,col="red", lwd=4)
for (i in 1:length(unique(qmatrix$Village))){
  #  abline(v = max(which(qmatrix$Dept==unique(qmatrix$Dept)[i])), col="red",lwd=4)
  segments(x0 = max(which(qmatrix$Village==unique(qmatrix$Village)[i])), 
           x1=max(which(qmatrix$Village==unique(qmatrix$Village)[i])),
           y0=-0.02, y1=1, col="blue", lwd=5,lty = 1)
}

for (i in 1:length(unique(qmatrix$Dept))){
#  abline(v = max(which(qmatrix$Dept==unique(qmatrix$Dept)[i])), col="red",lwd=4)
  segments(x0 = max(which(qmatrix$Dept==unique(qmatrix$Dept)[i])), 
           x1=max(which(qmatrix$Dept==unique(qmatrix$Dept)[i])),
                 y0=-0.02, y1=1, col="red", lwd=6)
  }

segments(x0=0, y0=-0.02,x1=0,y1=1,col="red",lwd=6)
#I exported this using width of 1800

#In order to get the data tables and make this reproducible I used the following:
#to get the frequencies of villages in the correct order:
# a = as.data.frame(table(qmatrix$Village))
# a[match(unique(qmatrix$Village),a$Var1),]
# #to get the data for the cross entropy graphs:
# allbest
# #to get the table of data for the graphs:
# qmatrix
#save(project_Dept,file = "project_Dept300rep150323.Robj") to save the file


