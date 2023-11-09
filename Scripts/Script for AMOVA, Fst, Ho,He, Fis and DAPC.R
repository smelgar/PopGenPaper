library(vcfR)
library(adegenet)
library(poppr)
library(mmod)
library(ape)
library(magrittr)
library(ggplot2)

setwd("/home/sergio/OneDrive/DocumentosData 1/Proyectos/Proyecto Salvador")
#vcf <- read.vcfR("batch_1NombreCorto (1).vcf", verbose = FALSE )
#vcf <- read.vcfR("batch_1NombreCorto_(1)_filtered (1).vcf", verbose = FALSE)
vcf <- read.vcfR("batch_1NombreCorto_(1)_filteredMAF0.01.vcf", verbose = FALSE)
vcfsamplenames <- colnames(vcf@gt)
#the object SNPSmsSnpInd has the list of selected individuals with missingness less than 0.2
#taken from the genetic filtering analysis using "rownames(SNPSmsSnpInd)":
samples50 <- c("A9944", "CHJ0007_L", "CHJ011_A", "CHJ083_A", "CHJ155_A", "CHJ188_A", "S174_A", "S236_A", 
               "S327_A", "S328_A", "S333_A", "S336_A", "S358_A", "S359_A", "S443_A", "S445_A", 
               "S446_A", "S447_A", "S448_A", "S451_A", "S455_A", "S476_A", "S478_A", "S486_A", 
               "S507a_A", "S507b_A", "S649_A", "TEX0010", "TEX0023", "TPG1017_A", "TPG154", "TPG715", 
               "TPS0053", "TPS0068", "TPS0180", "TPS0225", "TPS0330", "CHJ1016_L", "CHJ168_A", "CHJ169_A", 
               "CHJ172_A", "CHJ360_A", "CHJ362_A", "CHJ363_A", "CHJ364_A", "CHJ366_A", "CHJ367_A", "CHJ376_A", 
               "CHJ378_A", "CHJ386_A")
#these are the excluded samples:
setdiff(vcfsamplenames, samples50)



#vcf2 has here extraction of the individuals that are listed in the samples of SNPSmsSnpInd, and preserving
#the first column named "format"

vcf2 <- vcf[,c(1,match(samples50,vcfsamplenames))]
#write.vcf(vcf2, file = "batch_1NombreCorto_(1)_filteredMAF0.01ind50.vcf.gz")
#rownames(SNPSmsSnpInd)
#colnames(vcf@gt)[-1]
strata <-as.data.frame(read.csv2("GeneindStrata.csv",sep = ","))


strata2 <- strata[match(colnames(vcf2@gt)[-1],strata$Sample),]

#to correct the error I needed to change the vcf ID column to something unique https://www.biostars.org/p/305137/
#this was taken from filtering script: y$ID <- paste(y$X.CHROM,"_",y$POS,"_",y$ID,sep="")
vcf2@fix[,"ID"] = paste0(vcf@fix[,"CHROM"],"_",vcf@fix[,"POS"],vcf@fix[,"ID"])
#I used this before, and also works vcf2@fix[,"ID"] = paste0(vcf@fix[,"CHROM"],"_",vcf@fix[,"POS"],"_",1:2)
#head(vcf)
#vcf@fix[,"ID"]

my_genind <- vcfR2genind(vcf2)

# to see if there are duplicates:https://stackoverflow.com/questions/16905425/find-duplicate-values-in-r
n_occur = data.frame(table(vcf2@fix[,"ID"]))
n_occur[n_occur$Freq > 1,]
#Here the ID variable have only unique values.  Error corrected.

my_genind <- vcfR2genind(vcf2)
#names(as.data.frame(vcf@gt)) #to get sample names, but better use:
indNames(my_genind)
#to extract individuals see: https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2015-December/001388.html
#I got sample names as an object and the exported as csv, I got the GeneindStrata.csv from Google Docs and
#exported "GeneindStrata.csv".  Then I imported it.

#strata(my_genind) = as.data.frame(strata$VillageDepartmentCountry)
#splitStrata(my_genind, ~Country/Department/Village, sep = "_")

strata(my_genind) = data.frame(Country = strata2$Country,Department = strata2$DepartmentSA,Village =strata2$Village.1)

setPop(my_genind) <- ~Country/Department
#from https://github.com/mossmatters/Structure-Pipeline/blob/master/structureplot.R
#poppr(my_genind)
#To save my_genind: save(my_genind,file="genindfiltered50ind.robj"), taken from: https://lists.r-forge.r-project.org/pipermail/adegenet-forum/2015-August/001227.html
#to open the object: load("genindfiltered50ind.robj")
#https://grunwaldlab.github.io/Population_Genetics_in_R/Pop_Structure.html
#my_genind_dist <- provesti.dist(my_genind)
#par(mar = c(1,1,1,1),xpd=TRUE)
# theTree <- my_genind_dist %>%
#   nj() %>%    # calculate neighbor-joining tree
#   ladderize() # organize branches by clade
# plot(theTree)
# add.scale.bar(length = 0.05) # add a scale bar showing 5% difference.

# Analysis from https://grunwaldlab.github.io/Population_Genetics_in_R/Pop_Structure.html
par(mar = c(1,1,3,1),xpd=TRUE)
set.seed(999)
my_genind %>%
  genind2genpop(pop = ~Country/Department) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = nei.dist)


my_genind %>%
  genind2genpop(pop = ~Country/Department/Village) %>%
  aboot(cutoff = 50, quiet = TRUE, sample = 1000, distance = nei.dist)

#AMOVA 
My_Genind_amova <- poppr.amova(my_genind, ~Country/Department)
My_Genind_amova
My_Genind_amova$results
My_Genind_amova$componentsofcovariance
#from https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html:
My_Genind_amova_test    <- randtest(My_Genind_amova, nrepet = 999)
My_Genind_amova_test

#https://rdrr.io/cran/hierfstat/man/basic.stats.html
#calculate Fst between populations
library("hierfstat")
mybasicstats = basic.stats(my_genind)
mybasicstats$overall #to see the global Ho, Hs, HT, Dst, etc.
t(summary(mybasicstats$Ho)[4,])
t(summary(mybasicstats$Hs)[4,])
t(summary(mybasicstats$Fis)[4,])

#Mann Whitney U for Ho != Hs
for (i in 1:4) {
print(colnames(mybasicstats$Ho)[i]);
print(wilcox.test(mybasicstats$Ho[,i], mybasicstats$Hs[,i]))
}

#from: https://statsandr.com/blog/one-sample-wilcoxon-test-in-r/
for (i in 1:4) {
  print(colnames(mybasicstats$Fis)[i]);
  print(wilcox.test(mybasicstats$Fis[,i]
                    ))
}


#https://www.geeksforgeeks.org/how-to-find-confidence-intervals-in-r/
# R program to find the confidence interval
# Calculate the mean and standard error

model <- lm(mybasicstats$Ho ~ 1)
# Find the confidence interval of Ho
confint(model, level=0.95)

model <- lm(mybasicstats$Hs ~ 1)
# Find the confidence interval of Hs
confint(model, level=0.95)

model <- lm(mybasicstats$Fis ~ 1)
# Find the confidence interval of Fis
confint(model, level=0.95)
model

#calculate number of non NaN elements in Fis column https://stackoverflow.com/questions/24027605/determine-the-number-of-na-values-in-a-column
length(mybasicstats$Fis) - sum(is.na(mybasicstats$Fis))
length(mybasicstats$Ho) - sum(is.na(mybasicstats$Ho))
length(mybasicstats$Hs) - sum(is.na(mybasicstats$Hs))

par(mar=c(2, 2, 2, 2))
Ho = as.data.frame(mybasicstats$Ho)
par(mfrow = c(2, 2))
for(i in 1:length(Ho)) {
  hist(as.numeric(unlist(Ho[i]) ) ,main=paste0("Histogram of Ho of ",names(Ho)[i]),
       ylim = c(0, 800))
}
Hs = as.data.frame(mybasicstats$Hs)
for(i in 1:length(Hs)) {
  hist(as.numeric(unlist(Hs[i]) ) ,main=paste0("Histogram of Hs of ",names(Hs)[i]),
       ylim = c(0, 800))
}

#The overlapping Ho and He in the same plots.
#https://www.dataanalytics.org.uk/plot-two-overlapping-histograms-on-one-chart-in-r/#:~:text=To%20make%20sure%20that%20both,combined%20range%20of%20the%20samples.
ax <- pretty(0:1,n=12)
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
#reorder columns:
Ho <- Ho[,c(4,1,2,3)]
Hs <- Hs[,c(4,1,2,3)]
departmentsforgraph <- gsub("_"," ",colnames(Ho))
departmentsforgraph <- gsub("lS","l S",departmentsforgraph,ignore.case = FALSE)
departmentsforgraph <- gsub("eA","e-A",departmentsforgraph,ignore.case = FALSE)
departmentsforgraph <- gsub("Salvador","Salvador,",departmentsforgraph,ignore.case = FALSE)
departmentsforgraph <- gsub("Guatemala","Guatemala,",departmentsforgraph,ignore.case = FALSE)

#Histograms of Ho and Hs
par(mfrow = c(1, 1))
# #for common legend: https://stackoverflow.com/questions/10389967/common-legend-for-multiple-plots-in-r
# m <- matrix(c(1,2,3,1,4,5,1,6,6),nrow = 3,ncol = 3,byrow = TRUE)
# layout(mat = m,heights = c(0.4,0.4,0.2),widths = c(0.1,0.45,0.45))
# plot(1, type = "n", axes=N, xlab="", ylab="SNPs number", bty = 'n')
# text(x=1,y=1.1,"Number of SNPs",srt=90,cex=2)
# #legend(x="left", legend="SNPs number",cex=2, bty = "n",seg.len=1, xpd = TRUE, horiz = F)
# 
# for(i in 1:length(Ho)) {
#   histHo <- hist(as.numeric(unlist(Ho[i])), breaks = ax, plot = FALSE) #first histogram
#   histHs <- hist(as.numeric(unlist(Hs[i])), breaks = ax, plot = FALSE) #second histogram
#     plot(histHo, col = "yellow",ylim = c(0,800), main = paste0("\n\n\n",departmentsforgraph[i]),xlab="", ylab="") #https://stackoverflow.com/questions/33480500/is-there-a-way-to-show-overlapping-histograms-in-r-without-adjusting-transparenc
#   plot(histHs, col = "red",ylim = c(0,800), add = TRUE, xlab="", ylab="")
#   plot(histHo, col = "yellow",ylim = c(0,800), add = TRUE, xlab="", ylab="")  
#   plot(histHs, col = rgb(1,0,0,0.5),ylim = c(0,800), add = TRUE, xlab="", ylab="")
# 
# }
# 
# plot(1, type = "n", axes=FALSE, xlab="", ylab="", bty = 'n')
# plot_colors <- c("yellow","red")
# legend(x="top", y=NULL, legend = c("Ho","Hs"), fill = c("yellow","red"), 
#        cex = 1.8,horiz = TRUE,inset = -0,bty = "n" )
# legend(x="bottom", legend="Heterozygocities",cex=2, bty = "n",seg.len=1, xpd = TRUE)

#for common legend: https://stackoverflow.com/questions/10389967/common-legend-for-multiple-plots-in-r
m <- matrix(c(1,2,3,1,4,5,1,6,6),nrow = 3,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.4,0.4,0.2),widths = c(0.1,0.45,0.45))
plot(1, type = "n", axes=FALSE, xlab="", ylab="SNPs number", bty = 'n')
text(x=1,y=1.1,"Number of SNPs",srt=90,cex=2)
#legend(x="left", legend="SNPs number",cex=2, bty = "n",seg.len=1, xpd = TRUE, horiz = F)

for(i in 1:length(Ho)) {
  histHo <- hist(as.numeric(unlist(Ho[i])), breaks = ax, plot = FALSE) #first histogram
  histHs <- hist(as.numeric(unlist(Hs[i])), breaks = ax, plot = FALSE) #second histogram
  plot(histHo, col = "yellow",ylim = c(0,800), main = paste0("\n",departmentsforgraph[i]),
       xlab="", ylab="",cex.main=1.5) #https://stackoverflow.com/questions/33480500/is-there-a-way-to-show-overlapping-histograms-in-r-without-adjusting-transparenc
  plot(histHs, col = "red",ylim = c(0,800), add = TRUE, xlab="", ylab="")
  plot(histHo, col = "yellow",ylim = c(0,800), add = TRUE, xlab="", ylab="")  
  plot(histHs, col = rgb(1,0,0,0.5),ylim = c(0,800), add = TRUE, xlab="", ylab="")
  # if(i==2) {
  #   legend(x="topright",legend=c("Ho","Hs"),fill=c("yellow","red"),cex=1.8)
  # }
  
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="", bty = 'n')
plot_colors <- c("yellow","red")
legend(x="top", y=NULL, legend = c("Ho","Hs"), fill = c("yellow","red"), 
       cex = 1.8,horiz = TRUE,inset = -0.5,bty = "n" )
legend(x="top", legend="Heterozygocities",cex=2, bty = "n",seg.len=1, xpd = TRUE)












# #for common legend: https://stackoverflow.com/questions/10389967/common-legend-for-multiple-plots-in-r
# m <- matrix(c(1,2,3,1,4,5,1,6,6),nrow = 3,ncol = 3,byrow = TRUE)
# layout(mat = m,heights = c(0.4,0.4,0.2),widths = c(0.1,0.45,0.45))
# plot(1, type = "n", axes=N, xlab="", ylab="SNPs number", bty = 'n')
# text(x=1,y=1.1,"Number of SNPs",srt=90,cex=2)
# #legend(x="left", legend="SNPs number",cex=2, bty = "n",seg.len=1, xpd = TRUE, horiz = F)
# 
# for(i in 1:length(Ho)) {
#   histHo <- hist(as.numeric(unlist(Ho[i])), breaks = ax, plot = FALSE) #first histogram
#   histHs <- hist(as.numeric(unlist(Hs[i])), breaks = ax, plot = FALSE) #second histogram
#   plot(histHo, col = "yellow",ylim = c(0,800), main = paste0("\n\n\n",departmentsforgraph[i]),xlab="", ylab="") #https://stackoverflow.com/questions/33480500/is-there-a-way-to-show-overlapping-histograms-in-r-without-adjusting-transparenc
#   plot(histHs, col = "red",ylim = c(0,800), add = TRUE, xlab="", ylab="")
#   plot(histHo, col = "yellow",ylim = c(0,800), add = TRUE, xlab="", ylab="")  
#   plot(histHs, col = rgb(1,0,0,0.5),ylim = c(0,800), add = TRUE, xlab="", ylab="")
#   if(i==2) {
#     legend(x="topright",legend=c("Ho","Hs"),fill=c("yellow","red"),cex=1.8)
#   }
#   
# }
# 
# plot(1, type = "n", axes=FALSE, xlab="", ylab="", bty = 'n')
# #plot_colors <- c("yellow","red")
# #legend(x="top", y=NULL, legend = c("Ho","Hs"), fill = c("yellow","red"), 
# #       cex = 1.8,horiz = TRUE,inset = -0,bty = "n" )
# legend(x="top", legend="Heterozygocities",cex=2, bty = "n",seg.len=1, xpd = TRUE)



par(mfrow = c(1, 1))
FstDepart  <- pairwise.neifst(my_genind)
plot(hclust(as.dist(FstDepart)), main="Fst dendrogram between departments")
# to get confidence intervals for pairwise Fst:
boot.ppfst(my_genind)  #https://rdrr.io/cran/hierfstat/man/boot.ppfst.html


setPop(my_genind) <- ~Country/Department/Village
FstVillage  <- pairwise.neifst(my_genind) #write.csv(FstVillage,"FstVillage.csv")
#plot(hclust(as.dist(FstVillage))) #now it gives an error which is corrected below
#erase columns and rows with NAs https://stackoverflow.com/questions/54153670/remove-na-values-from-distance-matrix-in-r
FstVillageNA = FstVillage[rowSums(is.na(FstVillage)) == 1, colSums(is.na(FstVillage)) == 1, drop = FALSE]
plot(hclust(as.dist(FstVillageNA)), main="Fst dendrogram between villages")
#write.csv(FstVillageNA,"FstVillageNA.csv")

#there was an error since there are some NaNs in the matrix 
#https://github.com/jokergoo/ComplexHeatmap/issues/155
#I found this function: 
#' Hclust cannot handle matrices in which for some pairs of rows and columns,
#' only 1 or fewer shared values are non-NA. This function recurrently
#' identifies the most aggravating column/row, excludes that column/row and checks
#' whether more columns/rows need to be excluded
#'
#' @param mat Matrix to investigate
#' @param min_shared_fields Minimum number of positions that are not NA in both
#' vectors in order not to flag the vector pair as problematic
#'
identify_problematic_combs <- function(mat, min_shared_fields = 1) {
  exclude_rows <- NULL
  exclude_cols <- NULL
  stopifnot(is.matrix(mat))
  
  ## Loop over candidate removals
  for (k in 1:nrow(mat)) {
    candidate_rows <- setdiff(1:nrow(mat), exclude_rows)
    problem_row_combs <- NULL
    for (i in candidate_rows) {
      i_idx <- which(candidate_rows == i)
      for (j in candidate_rows[i_idx:length(candidate_rows)]) {
        if (sum(!is.na(mat[i, ]) & !is.na(mat[j, ])) <= min_shared_fields) {
          problem_row_combs <- rbind(problem_row_combs, c(i, j))
        }
      }
    }
    if (is.null(problem_row_combs)) break
    exclude_rows <- c(exclude_rows,
                      as.integer(names(which.max(table(problem_row_combs)))))
  }
  
  for (k in 1:ncol(mat)) {
    candidate_cols <- setdiff(1:ncol(mat), exclude_cols)
    problem_col_combs <- NULL
    for (i in candidate_cols) {
      i_idx <- which(candidate_cols == i)
      for (j in candidate_cols[i_idx:length(candidate_cols)]) {
        if (sum(!is.na(mat[, i]) & !is.na(mat[, j])) <= min_shared_fields) {
          problem_col_combs <- rbind(problem_col_combs, c(i, j))
        }
      }
    }
    if (is.null(problem_col_combs)) break
    exclude_cols <- c(exclude_cols,
                      as.integer(names(which.max(table(problem_col_combs)))))
  }
  
  return(list('row' = exclude_rows, 'column' = exclude_cols))
}


remove_problematic_combs <- function() {
  problematic_combs <- identify_problematic_combs(
    mat = mat, min_shared_fields = min_shared_fields)
  if (!is.null(problematic_combs$row)) {
    mat <- mat[-problematic_combs$row, ]
  }
  if (!is.null(problematic_combs$column)) {
    mat <- mat[, -problematic_combs$column]
  }
  return(mat)
}
formals(remove_problematic_combs) <- formals(identify_problematic_combs)
#use the functions to eliminate the problematic rows(columns)
identify_problematic_combs(FstVillage,10)
par(mar = c(3,3,3,3),xpd=TRUE)
FstVillage <- remove_problematic_combs(FstVillage,10)

#plot without three problematic villages
plot(hclust(as.dist(FstVillage)))

# #https://popgen.nescent.org/DifferentiationSNP.html
# #unsupervised clustering
# # using Kmeans and DAPC in adegenet 
# set.seed(20160308) # Setting a seed for a consistent result
# grp <- find.clusters(my_genind, max.n.clust = 10, n.pca = 20, choose.n.clust = FALSE) 
# names(grp)
# grp$grp
# par(mar=c(0,0,0,0))
# dapc1 <- dapc(my_genind, grp$grp, n.pca = 20, n.da = 6) 
# scatter(dapc1) # plot of the group

#https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
#DAPC using the Departments
#dapc <- dapc(my_genind, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(my_genind) - 1)
#scatter(dapc, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2, posi.da = "topright")

setPop(my_genind) <- ~Country/Department
dapc <- dapc(my_genind, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = nPop(my_genind) - 1)
scatter(dapc, cell = 0, pch = 18:23, cstar = 0, mstree = TRUE, lwd = 2, lty = 2, posi.da = "topright")

#from https://popgen.nescent.org/startMicrosatellite.html
div = summary(my_genind)
plot(div$Hobs, xlab="Loci number", ylab="Observed Heterozygosity", 
     main="Observed heterozygosity per locus")
plot(div$Hobs, div$Hexp, xlab="Observed Heterozygosity", ylab="Expected Heterozygosity", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus")

basicstat <- basic.stats(my_genind,diploid=TRUE,digits=4)
names(basicstat)
Fis = summary(basicstat$Fis)
Fis = as.data.frame(Fis)
names(Fis) = c("No", "Dept", "Fis")
Fis
Fis = Fis[grep("Mean",Fis$Fis),]
Fis$Fis = gsub("Mean   : ", "", Fis$Fis)

#to get the mean and the sd in two tables:
summaryHoHeFismean = matrix(, nrow = 4, ncol = 3)
summaryHoHeFisSD = matrix(, nrow = 4, ncol = 3)
for (j in 3:5) {
  for (i in 1:4) {
    print(paste(colnames(basicstat[[j]])[i],names(basicstat)[j]))
    print(paste("mean = ",mean(basicstat[[j]][,i],na.rm=T)))
    print(paste("sd = ",sd(basicstat[[j]][,i],na.rm=T)))
    summaryHoHeFismean[i,j-2]=mean(basicstat[[j]][,i],na.rm=T)
    summaryHoHeFisSD[i,j-2]=sd(basicstat[[j]][,i],na.rm=T)
  }}
colnames(summaryHoHeFismean) = names(basicstat[3:5])
rownames(summaryHoHeFismean) = colnames(basicstat[[3]])
summaryHoHeFismean =summaryHoHeFismean[c(4,1,3,2),]

colnames(summaryHoHeFisSD) = names(basicstat[3:5])
rownames(summaryHoHeFisSD) = colnames(basicstat[[3]])
summaryHoHeFisSD =summaryHoHeFisSD[c(4,1,3,2),]
summaryHoHeFismean = round(summaryHoHeFismean,4)
summaryHoHeFisSD = round(summaryHoHeFisSD,4)
summaryHoHeFismean
summaryHoHeFisSD


#DAPC taken from the tutorial chrome-extension://oemmndcbldboiebfnladdacbdfmadadm/https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
setPop(my_genind) <- ~Country/Department
grp <- find.clusters(my_genind, max.n.clust=40)
# I chose number of groups as 2 and 4.
table(pop(my_genind), grp$grp)
#table.value(table(pop(my_genind),grp$grp), col.lab=paste("inf", 1:4),row.lab=paste("ori", 1:4))
table.value(table(pop(my_genind),grp$grp), col.lab=paste("inf", 1:4),row.lab=unique(pop(my_genind)))

#erase strata data from genind object and add the pop$pop as the population
strata(my_genind) = NULL
pop(my_genind) = grp$grp


#To see number of PC that should be retained:
par(mar = c(5,5,1,1),xpd=TRUE)


dapc1 <- dapc(my_genind, grp = grp$grp,n.pca = 10, n.da = length(unique(grp$grp)-1))
#nevertheless, there is not a clear limit under 62.
#it is easier to use about 10 pcs retained and 4 DF retained.
#scatter(dapc1, posi.da = "topright")
#scatter(dapc1,1,1, col=myCol, bg="white",
#        scree.da=FALSE, legend=TRUE, solid=.4)
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:20, col = c("blue", "red", "darkgreen", "black"),
        cex=1,clab=0, leg=TRUE,txt.leg=paste("Cluster",1:4), cell=0, cstar=0)

#scatter(dapc(my_genind, var.contrib = TRUE, scale = FALSE, n.pca = 30, n.da = length(unique(grp$grp)-1)),
#        posi.da = "topright", grp = grp$grp)

#For different number of groups, in this case 10 groups:
#reset strata and pop data for genind object
strata(my_genind) = data.frame(Country = strata$Country,Department = strata$DepartmentSA,Village =strata$Village.1)

setPop(my_genind) <- ~Country/Department
#obtain 10 groups
grp <- find.clusters(my_genind, max.n.clust=40,n.pca = 62, n.clust = 10)

table(pop(my_genind), grp$grp)

table.value(table(pop(my_genind),grp$grp), col.lab=paste("inf", 1:10),row.lab=unique(pop(my_genind)))

#erase strata data from genind object and add the pop$pop as the population
strata(my_genind) = NULL
pop(my_genind) = grp$grp


#To see number of PC that should be retained:
par(mar = c(5,5,1,1),xpd=TRUE)

dapc2 <- dapc(my_genind, grp = grp$grp,n.pca = 10, n.da = (length(unique(grp$grp))-1))
scatter(dapc2, posi.da="bottomright", bg="white", pch =15:20 , col = c("blue", "red", "darkgreen", "black","pink","orange",
                                                                             "magenta","brown","brown3","purple"),
                                                                             cex=1,clab=0, leg=TRUE,txt.leg=paste("Cluster",unique(grp$grp)), cell=0, cstar=0)

scatter(dapc2, posi.da="bottomright", bg="white", pch =15:20 , col = c("blue", "red", "darkgreen", "black","pink","orange",
                                                                             "magenta","brown","brown3","purple"),
                                                                             cex=1,clab=0.5, leg=TRUE,txt.leg=paste("Cluster",unique(grp$grp)), cell=0, cstar=0)

#groups from 2 to 10

set.seed(100)

for (x in 3:10) {
  set.seed(100)
  grp <- find.clusters(my_genind, max.n.clust=40,n.pca = dim(dataDAPC)[1], n.clust = x)
  strata(my_genind) = data.frame(Country = strata2$Country,Department = strata2$DepartmentSA,Village =strata2$Village.1)
  setPop(my_genind) <- ~Country/Department
  table(pop(my_genind), grp$grp)
  
  table.value(table(pop(my_genind),grp$grp), col.lab=paste("inf", 1:10),row.lab=unique(pop(my_genind)))
  
  #erase strata data from genind object and add the pop$pop as the population
  strata(my_genind) = NULL
  pop(my_genind) = grp$grp
  
  
  par(mar = c(5,5,1,1),xpd=TRUE)
  
  dapc2 <- dapc(my_genind, grp = grp$grp,n.pca = 10, n.da = (length(unique(grp$grp))-1))
  
  # scatter(dapc2, posi.da="bottomright", bg="white", pch =15:20 , col = c("blue", "red", "darkgreen", "black","pink","orange",
  #                                                                        "magenta","brown","brown3","purple"),
  #         cex=1,clab=0, leg=TRUE,txt.leg=paste("Cluster",unique(grp$grp)), cell=0, cstar=0)
  # 
  scatter(dapc2, posi.da="bottomleft", bg="white", pch =15:20 , col = c("blue", "red", "darkgreen", "black","pink","orange",
                                                                              "magenta","brown","brown3","purple"),
                                                                              cex=1,clab=0.5, leg=TRUE,txt.leg=paste("Cluster",unique(grp$grp)), cell=0, cstar=0)
  # scatter(dapc2, posi.da="bottomright", bg="white", pch =15:20 , col = c("blue", "red", "darkgreen", "black","pink","orange",
  #                                                                       "magenta","brown","brown3","purple"),
  #         cex=1,clab=0.5, leg=TRUE,txt.leg=paste("Cluster",unique(grp$grp)), cell=0, cstar=0, xlim = c(-0.1,0.1),ylim = c(-5,5))
}
#######################################33333
x=3
set.seed(100)
grp <- find.clusters(my_genind, max.n.clust=40,n.pca = dim(dataDAPC)[1], n.clust = x)
strata(my_genind) = data.frame(Country = strata2$Country,Department = strata2$DepartmentSA,Village =strata2$Village.1)
setPop(my_genind) <- ~Country/Department
table(pop(my_genind), grp$grp)

table.value(table(pop(my_genind),grp$grp), col.lab=paste("inf", 1:10),row.lab=unique(pop(my_genind)))

#erase strata data from genind object and add the pop$pop as the population
strata(my_genind) = NULL
pop(my_genind) = grp$grp


par(mar = c(5,5,1,1),xpd=TRUE)

dapc2 <- dapc(my_genind, grp = grp$grp,n.pca = 10, n.da = (length(unique(grp$grp))-1))
############################*****************


dataDAPC <- dapc2$ind.coord
dataDAPC <- as.data.frame(dataDAPC)
dataDAPC$Dept <- as.factor(strata2$DepartmentSA)
#dataDAPC$DeptCol <- colors$DepartmentColor
dataDAPC$DeptCol <- rainbow(dim(dataDAPC)[1])
#plot(dataDAPC$LD1, dataDAPC$LD2, col=dataDAPC$DeptColor, pch=16)
# https://www.geeksforgeeks.org/how-to-color-scatter-plot-points-in-r/
#in a for loop, ggplot needs "print()": https://stackoverflow.com/questions/15678261/ggplot-does-not-work-if-it-is-inside-a-for-loop-although-it-works-outside-of-it
print(ggplot(dataDAPC, aes(x=LD1, y=LD2,group=Dept))+ 
        geom_point(aes(color=Dept)) + ggtitle(paste0("K = ",x)))# +xlim(-25,75)+ylim(-25,75)

#################################################33
#This is the DAPC graph with asterix pointing to El Carrizal
dataDAPC <- dapc2$ind.coord
dataDAPC <- as.data.frame(dataDAPC)
dataDAPC$Village <- as.factor(strata2$Village )
dataDAPC$Dept <- strata2$Department
dataDAPC$Country <- strata2$Country
dataDAPC$Dept[dataDAPC$Dept== "Santa_Ana"] <- "Santa Ana"
dataDAPC$DeptCol <-   ifelse(dataDAPC$Dept == "Jutiapa", "blue2",
                                                ifelse(dataDAPC$Dept == "Chiquimula", "dodgerblue3",
                                                ifelse(dataDAPC$Dept == "Sonsonate", "magenta",
                                                ifelse(dataDAPC$Dept == "Ahuachapan","orange4",
                                                ifelse(dataDAPC$Dept == "Santa Ana","deeppink3",NA)))))

dataDAPC$VillageIDs <- ifelse(dataDAPC$Village == "Santa_Ana", "",
                      ifelse(dataDAPC$Village == "El Zacatal, Amatepec", "EZa",
                      ifelse(dataDAPC$Village == "El Jute", "EJu",
                      ifelse(dataDAPC$Village == "El Aguacate, Azacualpa, Municipio de Armenia", "Aza",
                      ifelse(dataDAPC$Village == "Chilguyo, Texistepeque", "Chi",
                      ifelse(dataDAPC$Village == "Caserio Sarita Loma, Aldea La Primavera", "LPi",
                      ifelse(dataDAPC$Village == "Caserio San Raymundo, Canton Llano", "LDM",
                      ifelse(dataDAPC$Village == "Cantón Las Flores Cerro Alto, Municipio Caluco ", "CAl",
                      ifelse(dataDAPC$Village == "Llano Santa María", "LSM",
                      ifelse(dataDAPC$Village == "El_Chaperno", "ECh",
                      ifelse(dataDAPC$Village == "El Cerron, Olopa","ECe",
                      ifelse(dataDAPC$Village == "La Prensa, Olopa","LPr",
                      ifelse(dataDAPC$Village == "El_Carrizal","ECa",""
                      )))))))))))))
dataDAPC$group <- c("A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "A", "B", "A", "A", "B", "A", "B", "B", "A", "B", "B", "B", "B", "B", "B", "B", "B", "A", "A", "C", "C", "C", "A", "B", "A", "A", "A", "C", "A", "A", "A", "C", "C", "C", "C", "C", "C", "C", "C", "C")


dataDAPC$VillageCol = ifelse(dataDAPC$Village == "Santa Ana", "gray62",
                      ifelse(dataDAPC$Village == "El Zacatal, Amatepec", "hotpink",
                      ifelse(dataDAPC$Village == "El Jute", "orange3",
                      ifelse(dataDAPC$Village == "El Zacatal, Amatepec", "orange4",
                      ifelse(dataDAPC$Village == "El Jute", "magenta",
                      ifelse(dataDAPC$Village == "El Aguacate, Azacualpa, Municipio de Armenia", "red",
                      ifelse(dataDAPC$Village == "Chilguyo, Texistepeque", "coral1",
                      ifelse(dataDAPC$Village == "Caserio Sarita Loma, Aldea La Primavera", "deeppink",
                      ifelse(dataDAPC$Village == "Caserio San Raymundo, Canton Llano", "coral3",
                      ifelse(dataDAPC$Village == "Cantón Las Flores Cerro Alto, Municipio Caluco ", "coral4",
                      ifelse(dataDAPC$Village == "Llano Santa María", "cyan3",
                      ifelse(dataDAPC$Village == "El_Chaperno", "green",
                      ifelse(dataDAPC$Village == "El Cerron, Olopa","blue",
                      ifelse(dataDAPC$Village == "La Prensa, Olopa","blueviolet",
                      ifelse(dataDAPC$Village == "El_Carrizal","dodgerblue",
                      "gray62"
                      )))))))))))))))
dataDAPC$VillageShort <- ifelse(dataDAPC$Village == "Santa Ana", "",
                        ifelse(dataDAPC$Village == "El Zacatal, Amatepec", "El Zacatal",
                        ifelse(dataDAPC$Village == "El Jute", "El Jute",
                        ifelse(dataDAPC$Village == "El Aguacate, Azacualpa, Municipio de Armenia", "Azacualpa",
                        ifelse(dataDAPC$Village == "Chilguyo, Texistepeque", "Chilcuyo",
                        ifelse(dataDAPC$Village == "Caserio Sarita Loma, Aldea La Primavera", "La Primavera",
                        ifelse(dataDAPC$Village == "Caserio San Raymundo, Canton Llano", "Llano de Dona María",
                        ifelse(dataDAPC$Village == "Cantón Las Flores Cerro Alto, Municipio Caluco ", "Cerro Alto",
                        ifelse(dataDAPC$Village == "Llano Santa María", "Llano Santa María",
                        ifelse(dataDAPC$Village == "El_Chaperno", "El Chaperno",
                        ifelse(dataDAPC$Village == "El Cerron, Olopa","El Cerron",
                        ifelse(dataDAPC$Village == "La Prensa, Olopa","La Prensa",
                        ifelse(dataDAPC$Village == "El_Carrizal","El Carrizal",""
                        )))))))))))))
dataDAPC$VillageCol 
#change the order of departmentshttps://www.datanovia.com/en/blog/how-to-change-ggplot-legend-order/
dataDAPC$Dept <- factor(dataDAPC$Dept, levels = c("Chiquimula" ,"Jutiapa" ,"Ahuachapan" ,   "Santa Ana" , "Sonsonate" ))
dataDAPC$DeptCol <- factor(dataDAPC$DeptCol, levels = c("dodgerblue" ,"blue2" ,"orange4" ,   "deeppink3" , "magenta" ))

dataDAPC$VillageShort <- factor(dataDAPC$VillageShort, levels = (c("El Cerron",
                                                                   "La Prensa",
                                                                 "El Carrizal",
                                                                 "El Chaperno",
                                                                 "Llano Santa María",
                                                                 "Llano de Dona María",
                                                                 "Chilcuyo",
                                                                 "El Jute",
                                                                 "El Zacatal",
                                                                 "La Primavera",
                                                                 "Azacualpa",
                                                                 "Cerro Alto",
                                                                 NA
                                                                 )))
dataDAPC$legend <- paste(dataDAPC$VillageShort, dataDAPC$VillageIDs)
#to change NA in dataDAPC$legend to "El Salvador":
dataDAPC$legend[is.na(dataDAPC$legend)]<-"El Salvador"
dataDAPC$legend[dataDAPC$legend=="NA "]<-"El Salvador"

dataDAPC$legend <- factor(dataDAPC$legend, levels = (c("El Cerron ECe",
                                                                   "La Prensa LPr",
                                                                   "El Carrizal ECa",
                                                                   "El Chaperno ECh",
                                                                   "Llano Santa María LSM",
                                                                   "Llano de Dona María LDM",
                                                                   "Chilcuyo Chi",
                                                                   "El Jute EJu",
                                                                   "El Zacatal EZa",
                                                                   "La Primavera LPi",
                                                                   "Azacualpa Aza",
                                                                   "Cerro Alto CAl",
                                                                   "El Salvador"
)))

dataDAPC$VillageCol <- factor(dataDAPC$VillageCol, levels = (c("blue",
                                                               "blueviolet",
                                                               "dodgerblue",
                                                               "green",
                                                               "cyan3",
                                                               "coral3",
                                                               "coral1",
                                                               "orange3",
                                                               "hotpink",
                                                               "deeppink",
                                                               "red",
                                                               "coral4",
                                                               "gray"
)))

#https://www.statology.org/r-add-column-to-data-frame-based-on-other-columns/
dataDAPC$Villagepch <- with(dataDAPC,ifelse(dataDAPC$Village=="El_Carrizal",8,19))
#dataDAPC$DeptCol <- colors$DepartmentColor
#dataDAPC$DeptCol <- rainbow(dim(dataDAPC)[1])
#plot(dataDAPC$LD1, dataDAPC$LD2, col=dataDAPC$DeptColor, pch=16)
# https://www.geeksforgeeks.org/how-to-color-scatter-plot-points-in-r/
#in a for loop, ggplot needs "print()": https://stackoverflow.com/questions/15678261/ggplot-does-not-work-if-it-is-inside-a-for-loop-although-it-works-outside-of-it
library(ggrepel)
library(ggforce)
image = ggplot(dataDAPC, aes(x=LD1, y=LD2,group=Dept,label=dataDAPC$VillageIDs))+ 
#  ggplot(dataDAPC, aes(x=LD1, y=LD2,group=Dept,label=dataDAPC$group))+ 
        geom_point(aes(color=legend), #shape=dataDAPC$Villagepch, 
                   cex=2.5, stroke = 1) + #stroke: https://stackoverflow.com/questions/51644418/how-to-change-the-linewidth-of-point-using-the-geom-point
        scale_color_manual(values = levels(unique(dataDAPC$VillageCol)))+ #https://stackoverflow.com/questions/15965870/fill-and-border-colour-in-geom-point-scale-colour-manual-in-ggplot
        ggtitle(paste0("K = ",x)) + # +xlim(-25,75)+ylim(-25,75) 
  labs(colour = "Village")+
    theme(text = element_text(size = 15))
#print(image)
image + geom_text_repel() + annotate(geom="text", x=-2, y=3, label="A", cex=10,color="black")+
              annotate(geom="text", x=-5, y=0.5, label="B", cex=10,color="black")+
              annotate(geom="text", x=6, y=1, label="C", cex=10,color="black") +
              stat_ellipse(mapping = aes(group=group),color="darkgray",level = 0.9) +  #https://r-charts.com/correlation/scatter-plot-ellipses-ggplot2/
              xlab("DAPC1 (74.56%)")+ ylab("DAPC2 (25.44%)")
# to get Ho and Hs from a genind object by village:
strata(my_genind) = data.frame(Country = strata2$Country,Department = strata2$DepartmentSA,Village =strata2$Village.1)
strata(my_genind)
setPop(my_genind) = ~Country/Department/Village
pop(my_genind)
Hs(my_genind)
Ho(my_genind)
VillageHsHo = data.frame(Hs = Hs(my_genind), Ho = Ho(my_genind))
stats = basic.stats(my_genind)
VillageHsHoFis = data.frame(Hs=colMeans(stats$Hs,na.rm = TRUE), #https://stackoverflow.com/questions/21807987/calculate-the-mean-for-each-column-of-a-matrix-in-r
                            Ho=colMeans(stats$Ho,na.rm = TRUE),
                            Fis=colMeans(stats$Fis,na.rm = TRUE))
VillageHsHoFis #this is for Mantel test.
# library(datapasta)
# df_paste()

#  https://stackoverflow.com/questions/52171700/r-draw-ellipse-in-scatterplot-not-stat-ellipse
#  geom_ellipse(aes(x0=0, y0=0, a=1, b=2,angle=30))#+geom_polygon()

  
# image + geom_label_repel(box.padding   = 0.35, #https://stackoverflow.com/questions/15624656/label-points-in-geom-point
#                          point.padding = 0.5,
#                          segment.color = 'grey50')
#image + geom_text(check_overlap = TRUE)
#ggsave(file="DAPC20230510.svg", plot=image, width=8, height=6)
#################################################


# strata(my_genind) = data.frame(Country = strata2$Country,Department = strata2$DepartmentSA,Village =strata2$Village.1)
# setPop(my_genind) <- ~Country/Department
# 
# 
# #read csv "colors" to get information about houses  and departments:
# colors <- read.csv("colors.csv")
# par(mar = c(5,5,5,5),xpd=TRUE)
# 
# #Graphs with departmente based coloring
# for (x in 3:10) {
#  
# grp <- find.clusters(my_genind, max.n.clust=40,n.pca = dim(dataDAPC)[1], n.clust = x )
# #erase strata data from genind object and add the pop$pop as the population
# strata(my_genind) = NULL
# pop(my_genind) = grp$grp
# table(pop(my_genind), grp$grp)
# dapc2 <- dapc(my_genind, grp = grp$grp,n.pca = 10, n.da = (length(unique(grp$grp))-1))
# dataDAPC <- dapc2$ind.coord
# dataDAPC <- as.data.frame(dataDAPC)
# dataDAPC$Dept <- as.factor(strata2$DepartmentSA)
# #dataDAPC$DeptCol <- colors$DepartmentColor
# dataDAPC$DeptCol <- rainbow(dim(dataDAPC)[1])
# #plot(dataDAPC$LD1, dataDAPC$LD2, col=dataDAPC$DeptColor, pch=16)
# # https://www.geeksforgeeks.org/how-to-color-scatter-plot-points-in-r/
# #in a for loop, ggplot needs "print()": https://stackoverflow.com/questions/15678261/ggplot-does-not-work-if-it-is-inside-a-for-loop-although-it-works-outside-of-it
# print(ggplot(dataDAPC, aes(x=LD1, y=LD2,group=Dept))+ geom_point(aes(color=Dept)) + ggtitle(paste0("K = ",x)))# +xlim(-25,75)+ylim(-25,75)
# #print(ggplot(dataDAPC, aes(x=LD1, y=LD2,group=Dept))+ geom_point(aes(color=Dept)) +xlim(-5,5)+ylim(-5,5)+ ggtitle(paste0("K = ",x)))
# }
# #warnings are due to some points taken out of the graphs when zooming in.

#things to try: confidence intervals https://rdrr.io/cran/hierfstat/man/boot.vc.html,
#https://grunwaldlab.github.io/Population_Genetics_in_R/Pop_Structure.html
#https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html

