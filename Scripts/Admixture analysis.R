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
krange = 1:10
nrep = 200
project_Dept = snmf(genofile, K = krange, 
                    entropy = TRUE, repetitions = nrep, ploidy = 2, project = "new")
par(mar=(c(4,4,4,4)))
plot(project_Dept, col = "blue", pch = 18, cex = 2, 
     main = "Cross-entropy")
allbest = NULL
for (x in krange) {
  allbest[x] = min(cross.entropy(project_Dept,K=x))
  
}
# plot(allbest, col="blue", type="b", Main = "Best entropy values for each k", ylab = "entropy",
#      xlab = "k", pch = 18, cex = 2)
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

rownames(qmatrix) = ids
qmatrix = as.data.frame(qmatrix)
qmatrix$Dept = Dept
qmatrix$Country = Country
#oder from https://stackoverflow.com/questions/1296646/sort-order-data-frame-rows-by-multiple-columns
qmatrix = qmatrix[with(qmatrix,order(Country,Dept)),]



#vertical labels from:https://stackoverflow.com/questions/10286473/rotating-x-axis-labels-in-r-for-barplot
par(mar=c(10,8,0.5,0.5))
par(cex=1.5)
barplot(t(qmatrix[,1:ap]), col=RColorBrewer::brewer.pal(9,"Paired"), 
        #border=NA, 
        space=0, 
        #xlab="Individuals", 
        ylab="Admixture coefficients", las=2 )
#Add population labels to the axis:
#https://stackoverflow.com/questions/62408799/is-there-a-way-to-add-a-x-axis-below-the-primary-x-axis-in-r
for (i in 1:length(unique(qmatrix$Dept))){
  axis(1, at=median(which(qmatrix$Dept==unique(qmatrix$Dept)[i]))-0.5, labels=unique(qmatrix$Dept)[i], line=5, 
       col="black", col.axis = "black", cex.lab = 3)}
mtext("Individuals", 1, line=1, at=-4)
mtext("Dept", 1, line=5, at=-6)
# the following abline sometimes draws long lines, so I changed it to "segments()" https://r-charts.com/es/r-base/segments/
# for (i in 1:length(unique(qmatrix$Dept))){
#   abline(v = max(which(qmatrix$Dept==unique(qmatrix$Dept)[i])), col="red",lwd=4)}
# abline(v = 0,col="red", lwd=4)

for (i in 1:length(unique(qmatrix$Dept))){
#  abline(v = max(which(qmatrix$Dept==unique(qmatrix$Dept)[i])), col="red",lwd=4)
  segments(x0 = max(which(qmatrix$Dept==unique(qmatrix$Dept)[i])), x1=max(which(qmatrix$Dept==unique(qmatrix$Dept)[i])),
                 y0=-0.02, y1=1, col="red", lwd=4)
  }

segments(x0=0, y0=-0.02,x1=0,y1=1,col="red",lwd=4)
#I exported this using width of 2000

