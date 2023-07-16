############################################################################################
##########################  Filtering Genotype Matrix
############################################################################################

## Download Required Packages

library(raster) # to extract coordinates from our shapefile
library(vcfR) # to open VCF in R. 
setwd("/home/sergio/OneDrive/DocumentosData 1/Proyectos/Proyecto Salvador")
## Load Data
Archivo <- "batch_1NombreCorto (1).vcf"
#Archivo <- "batch_1NombreCorto_(1)_filtered.vcf"
#Archivo <- "batch_3.vcf"
#Archivo <- "batch_1NombreCorto.vcf"
IDunicos <- paste0("IDunicos",Archivo)
#Esto es para colocar la columna "ID" con valores únicos, si no, no corre.
y <- read.csv(Archivo, sep = "\t",skip=9)
y$ID <- paste(y$X.CHROM,"_",y$POS,"_",y$ID,sep="")

#luego se agregan las primeras líneas y se sustituye la línea del encabezado para que quede el archivo "vcf" listo. Eso lo hice manualmente.

encabezado <- readLines(Archivo,10)
write(encabezado,IDunicos)
write.table(y,IDunicos, col.names=FALSE, row.names=FALSE,append=TRUE,quote=FALSE,sep="\t")


#write.table(y,"IDunicos3.vcf",sep = "\t",quote = FALSE,row.names = F)

VCF <- read.vcfR(IDunicos) # import the VCF in R

#luego manualmente junté las primeras 9 líneas y eliminé la primera columna con ^\d{1,4}\s


### Create a Genotype Matrix

## The VCF object contains several additional information concerning genetic variant quality that are not needed in this exercise. 
## We want to obtain a Genotype Matrix out of the VCF object, we can do it as shown here below:  

GT <- vcfR2loci(VCF)

GT[1:10,1:10] # The GT matrix shows the genotype for each individual (row) and genetic variant (column)

dim(GT) # 62 individuals,  genetic markers 

## For our calculations, we prefer having the genotypes coded as 0,1,2 (homozygous 0/0, heterozygous, homozygous 1/1).
## We can perform the conversion as shown here below.

SNPS <- as.matrix(GT)

SNPS[SNPS=='0/0'] = 0
SNPS[SNPS=='0/1'] = 1
SNPS[SNPS=='1/0'] = 1
SNPS[SNPS=='1/1'] = 2

SNPS <- apply(SNPS,2,as.numeric) # since the SNPS matrix used to be a character matrix, we have to convert it. 
rownames(SNPS) <- rownames(GT) # we also assign rownames as those in the GT matrix. 

## The SNPS matrix has the genotype values coded as 0,1,2:

SNPS[1:10,1:10]

### Check for genetic diversity metrics

# Before continuing, we use the save() function to register the SNPS object as a file. 
#save(SNPS, file='SNPS.robj')

# In case you want to start again the filtering procedure, you can regenerate the SNPS object calling the load function. 
#load('SNPS.robj') # this avoids you re-processing the vcf object. 

## There are R packages that allow to compute this metrics directly with built-in function.
## In this exercise we use the hand-written formula to better show the meaning of each procedure.

### Missingness by SNP
## For each SNP, it computes the % of missing observation across the individuals. 
#para grabar SNPS como archivo de excel:
#library(openxlsx)
#write.xlsx(SNPS,"genotipos.xlsx")


MNsnp <- apply(SNPS, 2, function(x) { sum(is.na(x))/length(x) })
hist(MNsnp, breaks=30,ylab="Frequency", xlab="Proportion of missing observations by SNP" )
#hist(MNsnp, breaks=10) # We can observe that most of the SNPS have a Missingness below 0.01 (i.e. ~ 2 individual ot of 160).
SNPS <- SNPS[,MNsnp<0.2] # We apply a 20% threshold

### Missingness by Individual
## For each individual, it computes the % of missing observation across the SNPs. 
MNind <- apply(SNPS, 1, function(x) { sum(is.na(x))/length(x) })
barplot(MNind,las=2)
hist(MNind, breaks=100,xlab = "Missigness per individual after filtering SNPs with missingness above 0.2") # We can observe that all the missingness by Individual are below 1% (i.e. ~ 350 snps out of 35000)

SNPS <- SNPS[MNind<0.2,] # We apply a 20% threshold
dim(SNPS)


# from here on, they are filtering for minor allele frequencies less than 0.05 and major genotype frequency less than 0.95

### Minor Allele Frequency
## For each SNP, it computes the frequency of the most rare allele across all individuals. 

MAF <- apply(SNPS, 2, function(x) {
  a=(sum(x=='1', na.rm=T)+sum(x=='0', na.rm=T)*2)/(sum(is.na(x)==F)*2)
  A=(sum(x=='1', na.rm=T)+sum(x=='2', na.rm=T)*2)/(sum(is.na(x)==F)*2)
  return(min(a,A))
})
hist(MAF, breaks=100) # The distribution has a left skew. RQ1) What genotypes would you filter out? Why?
SNPS <- SNPS[,MAF>0.05] # We apply a 5% threshold. 

### Major Genotype Frequency
## For each SNP, it computes the frequency of the most frequent genotype across all individuals. 

MGF <- apply(SNPS, 2, function(x) {
  aa=(sum(x=='0', na.rm=T))/(sum(is.na(x)==F))
  aA=(sum(x=='1', na.rm=T))/(sum(is.na(x)==F))
  AA=(sum(x=='2', na.rm=T))/(sum(is.na(x)==F))
  return(max(aa,aA,AA))
})
hist(MGF, breaks=100) # RQ2) What genotypes would you filter out? Why?
SNPS <- SNPS[,MGF<0.95] # We apply a 95% threshold 


dim(SNPS) # The filtering reduced the number of SNPS of ~2500, while no individual was filtered out. 

# We register the filtered genetic matricx as an rObject file so that we can use it in the next steps of the study. 
fSNPS=SNPS
save(fSNPS, file='fSNPS.robj')

