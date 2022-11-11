### Based on script by Basten Snoek, edited by Jesse Küpers

#load metadata files
load(file="C:\\Kupers et al 2022 RNA seq data\\obj_enrGO.out")

load(file="C:\\Kupers et al 2022 RNA seq data\\obj_GOnames.out")
rownames(GOnames) <- GOnames[,1]

mydata <- read.delim(file="C:\\Kupers et al 2022 RNA seq data\\output_readCounts_normalized.txt",row.names = 1)
all.genes <- unique(rownames(mydata))


#Make required functions for enrichment
get.all.count <- function(all.genes){ 
     all.count <- table(enrGO$GO[enrGO$geneid %in% all.genes])
     return(all.count)
}


get.go.list <- function(all.count){
     GO.id <- names(all.count)
     GO.desk <- matrix(NA,length(GO.id),3)
     for ( i in 1:length(GO.id)){ GO.desk[i,] <- c(GOnames[GOnames[,1] == GO.id[i],],NA,NA,NA)[1:3]}
     return(GO.desk)
}


go.test <- function(selected.genes,enrGO,all.genes,all.count){
     selc.count <- table(enrGO$GO[enrGO$geneid %in% selected.genes])
     tot <- length(all.genes)
     drawn <- length(selected.genes)
     p.hyp <- phyper(selc.count,all.count,tot-all.count,drawn,lower.tail = F)
     GOgenes <- as.numeric(all.count)
     Setgenes <- as.numeric(selc.count)
     pval <- as.numeric(p.hyp)
     GO.id <- names(all.count)
     GO.type <- GO.desk[,2] ; GO.name <- GO.desk[,3] 
     enr.test <- data.frame(GO.id,GO.type,GO.name,GOgenes,Setgenes,tot,drawn,pval)
     return(enr.test)
}


signif.go <- function(enr.test,p.thr,set.thr){
     these.go <- enr.test[enr.test$Setgenes>3 & enr.test$pval<0.001,]
     these.go <- these.go[order(these.go$pval),]
     return(these.go)
}


all.count <- get.all.count(all.genes)  


GO.desk <- get.go.list(all.count)       


#Perform GO analyses

# Use a list of Arabidopsis Gene Identifiers (e.g. "AT5G18010")      
selected.genes <-   c("AT5G43190", "AT4G38410", "AT5G47370", "AT5G18050", "AT1G15580", "AT2G22810", "AT3G50340", "AT3G03830", "AT4G37770", "AT1G52830",
"AT1G65920", "AT5G18030", "AT3G25900", "AT1G29490", "AT4G38850", "AT5G18010", "AT2G18010", "AT3G50350", "AT5G66580", "AT1G26945",
"AT3G60630", "AT5G18020", "AT3G18200", "AT5G43700", "AT5G04190", "AT5G18080", "AT4G03140", "AT5G39860", "AT1G67900", "AT5G65800",
"AT1G78100", "AT3G59900", "AT3G12670", "AT5G18060", "AT4G34770", "AT1G30190", "AT3G03840", "AT3G03850", "AT1G29500", "AT1G04610") 

enr.test <- go.test(selected.genes,enrGO,all.genes,all.count)
enriched.GO<- signif.go(enr.test)

#Select only biological processes
enriched.GO <- dplyr::filter(enriched.GO, grepl("biological_process",GO.type)) 
View(enriched.GO)
