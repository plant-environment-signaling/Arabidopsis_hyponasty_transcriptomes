### Based on script by Basten Snoek, edited by Jesse Küpers

#DEG Analyses used in paper are #1 Treatment effect per timepoint, per tissue (Model 1) (Figure 1, Data S1)
                                #2 Treatment * Tissue effect over relevant timepoints (100-300 minutes) in the petiole samples (Model 2) (Figure 2, Data S2). Used for GO enrichment on subsetted genes
                                #3 Tissue effect in WL only between petiole samples over all timepoints (Model 4) (Data S3)

# Open data and metadata
mydata <- read.delim(file="C:\\Kupers et al 2022 RNA seq data\\output_readCounts_normalized.txt",row.names = 1)
mydata[1:5,]
my.info <- read.delim(file="C:\\Kupers et al 2022 RNA seq data\\MetaData.txt")
my.info[1:5,]
rownames(my.info) <- my.info$SampleName
my.info <- my.info[colnames(mydata),]

#1) treatment effect per timepoint per tissue ####

# select genes
meanreads <- apply(mydata,1,mean)  
ndat <- mydata[meanreads>1,]       #select genes that have average read number >1
odat <- ndat[((apply(ndat>0,1,sum))>(length(ndat)/6)),] #Selects genes with reads in > one-sixth (one tissue in one treatment) of the samples
pdat <- (odat+1)    # Data adjustment to deal with 0 read samples

#Create lists per timepoint (ptp) of p-value & mean expression in WL vs FRtip for the three tissues independently
#The resulting files are used to generate the heatmap in Figure 1 as well as Data S1.
ptp <- NULL
ptp$Tip <- NULL
ptp$Abax <- NULL
ptp$Adax <- NULL

uni.tis <- as.character(unique(my.info$Tissue))     
uni.tp <- unique(my.info$Timepoint)

#Leaf Tip
use.tis <- uni.tis[1]
for ( k in 1:9){
  use.tp <- uni.tp[k]                           
  selc <- my.info$Timepoint == use.tp & my.info$Tissue == use.tis
  
  trt <- my.info$Treatment[selc]
  pmat <- matrix(NA,ncol = 4,nrow=nrow(pdat[,selc]))
  for ( i in 1:nrow(pdat)){
    gexp <- as.numeric(pdat[i,selc])
    tmp <- anova(lm(gexp~trt))
    pmat[i,1] <- tmp$`Pr(>F)`[1]
  }
  pmat[,2] <- apply(pdat[,selc & my.info$Treatment == "WL"],1,mean)
  pmat[,3] <- apply(pdat[,selc & my.info$Treatment == "Frtip"],1,mean)
  pmat[,4] <- log2(pmat[,3]/pmat[,2])
  ptp$Tip[[k]] <- pmat 
  png(filename = paste(use.tis,use.tp, ".png", sep=""), units="in", width=6, height=6, res=300)
  plot(ptp$Tip[[k]][,4],-log10(ptp$Tip[[k]][,1]),main = paste(use.tis,use.tp), abline(2,0, col="red"), xlab = "log2 Fold Change", ylab = "-log10 p-value")
  dev.off()
}


#Abaxial petiole
use.tis <- uni.tis[2]
for ( k in 2:9){        #During library prep 40 min petiole samples were contaminated, therefore not sequenced.
  use.tp <- uni.tp[k]                           
  selc <- my.info$Timepoint == use.tp & my.info$Tissue == use.tis
  
  trt <- my.info$Treatment[selc]
  pmat <- matrix(NA,ncol = 4,nrow=nrow(pdat[,selc]))
  for ( i in 1:nrow(pdat)){
    gexp <- as.numeric(pdat[i,selc])
    tmp <- anova(lm(gexp~trt))
    pmat[i,1] <- tmp$`Pr(>F)`[1]
  }
  pmat[,2] <- apply(pdat[,selc & my.info$Treatment == "WL"],1,mean)
  pmat[,3] <- apply(pdat[,selc & my.info$Treatment == "Frtip"],1,mean)
  pmat[,4] <- log2(pmat[,3]/pmat[,2])
  ptp$Abax[[k]] <- pmat 
  png(filename = paste(use.tis,use.tp, ".png", sep=""), units="in", width=6, height=6, res=300)
  plot(ptp$Abax[[k]][,4],-log10(ptp$Abax[[k]][,1]),main = paste(use.tis,use.tp), abline(2,0, col="red"), xlab = "log2 Fold Change", ylab = "-log10 p-value")
  dev.off()
}


#Adaxial petiole
use.tis <- uni.tis[3]
for ( k in 2:9){        #During library prep 40 min petiole samples were contaminated, therefore not sequenced.
  use.tp <- uni.tp[k]                           
  selc <- my.info$Timepoint == use.tp & my.info$Tissue == use.tis
  
  trt <- my.info$Treatment[selc]
  pmat <- matrix(NA,ncol = 4,nrow=nrow(pdat[,selc]))
  for ( i in 1:nrow(pdat)){
    gexp <- as.numeric(pdat[i,selc])
    tmp <- anova(lm(gexp~trt))
    pmat[i,1] <- tmp$`Pr(>F)`[1]
  }
  pmat[,2] <- apply(pdat[,selc & my.info$Treatment == "WL"],1,mean)
  pmat[,3] <- apply(pdat[,selc & my.info$Treatment == "Frtip"],1,mean)
  pmat[,4] <- log2(pmat[,3]/pmat[,2])
  ptp$Adax[[k]] <- pmat 
  png(filename = paste(use.tis,use.tp, ".png", sep=""), units="in", width=6, height=6, res=300)
  plot(ptp$Adax[[k]][,4],-log10(ptp$Adax[[k]][,1]),main = paste(use.tis,use.tp), abline(2,0, col="red"), xlab = "log2 Fold Change", ylab = "-log10 p-value")
  dev.off()
}


#Create ".out" files that contain data for all timepoints and can be reloaded when necessary
Tip.out <- ptp$Tip
Abax.out <- ptp$Abax
Adax.out <- ptp$Adax



#2) model testing the effect of treatment, tissue and time on selected timepoints and selected tissues for Figure 2#### 
use.set <- pdat

#Select tissues and timepoints 
selc <- my.info$Tissue %in% c("Abax","Adax")
uni.time <- unique(my.info$Timepoint)
selc2 <- my.info$Timepoint %in% uni.time[4:9]  #4:9 = 100 - 300 minutes, this removes 40 - 80 where there is no petiole FR effect

#Get p-values for 3-way interaction per gene
mtis <- my.info$Tissue[selc & selc2] 
mtime <- as.factor(my.info$Timepoint[selc & selc2])
mtreat <- as.factor(my.info$Treatment[selc & selc2])
xm <- as.numeric(use.set[1,selc & selc2])
tmp <- anova(lm(xm~mtis*mtime*mtreat))

#Make get.p function to get p-values
get.p <- function(xm,mtis,mtime,mtreat)  {  
  tmp <- anova(lm(xm~mtis*mtime*mtreat))
  return(tmp$`Pr(>F)`)
}

all.int.p <- t(apply(use.set[,selc & selc2],1,get.p,mtis,mtime,mtreat))
all.int.p <- as.data.frame(all.int.p) 
colnames(all.int.p) <- rownames(tmp)
all.int.p[is.na(all.int.p)] <- 1

#Make a large datafile "tistreat" that has mean expression value for abax & adax in 2 treatments with the full ANOVA results
#This files is used to generate Figure 2 and Data S2 
wlabax <- as.data.frame(apply(pdat[,my.info$Tissue == "Abax" & selc2 & my.info$Treatment == "WL"],1,mean))
wladax <- as.data.frame(apply(pdat[,my.info$Tissue == "Adax" & selc2 & my.info$Treatment == "WL"],1,mean))
frabax <- as.data.frame(apply(pdat[,my.info$Tissue == "Abax" & selc2 & my.info$Treatment == "Frtip"],1,mean))
fradax <- as.data.frame(apply(pdat[,my.info$Tissue == "Adax" & selc2 & my.info$Treatment == "Frtip"],1,mean))

wlabax$gene <- rownames(wlabax)
wladax$gene <-rownames(wladax)
frabax$gene <- rownames(frabax)
fradax$gene <-rownames(fradax)

allexpr<- merge(wlabax, wladax, all=TRUE, by="gene")
allexpr<- merge(allexpr, frabax, all=TRUE, by="gene")
allexpr<- merge(allexpr, fradax, all=TRUE, by="gene")
colnames(allexpr) <- c("gene", "wlabax", "wladax", "frabax", "fradax")

all.int.p$gene <- rownames(all.int.p)
tistreat<- merge(allexpr, all.int.p, all= TRUE, by= "gene")
rownames(tistreat) <- tistreat$gene
tistreat<- tistreat[,2:12]



#3) Get abaxial vs adaxial differences in WL for Data S3####
selc <- my.info$Tissue %in% c("Abax","Adax")
selc2 <-my.info$Treatment %in% c("WL")
mtis <- my.info$Tissue[selc & selc2] 
mtime <- as.factor(my.info$Timepoint)[selc & selc2]
mtreat <- as.factor(my.info$Treatment)[selc]
xm <- as.numeric(pdat[1,selc&selc2])
tmp <- anova(lm(xm~mtis*mtime))

# Get p-values with edited get.p function
get.p <- function(xm,mtis,mtime)  {  
  tmp <- anova(lm(xm~mtis*mtime))
  return(tmp$`Pr(>F)`)
}

all.p <- apply(pdat[,selc & selc2],1,get.p,mtis,mtime)
all.p <- t(all.p)

#Adjust p-values
all.p.adj <- all.p
for (i in 1:4) {
  all.p.adj[,i]<- p.adjust(all.p[,i])
}

# Turn p-values of all NA's into 1
all.p.adj[is.na(all.p.adj)] <- 1                    

colnames(all.p.adj) <- rownames(tmp)
all.p.adj <- as.data.frame(all.p.adj)
View(all.p.adj)