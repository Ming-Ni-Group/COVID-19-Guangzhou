# Raxml tree and global iSNV plot 
rm(list=ls())
library(stringr)
library(pheatmap)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(treeio)
library(dplyr)
library(aplot)
# Loading
options(stringsAsFactors = F)
snppositionfile <- "./result/Samples_all_iSNVs.csv"
positionfile <- "./data/Position.txt"
treefile <- "./result/RAxML_bestTree.raxmltree"
sampleinfofile <- "./data/High_coverage_genome_sequences_information.csv"
snpannofile <- "./data/Mutation_annotion.txt"

## position loading
positionlist <- read.table(positionfile)[, 1]
##SNP loding
snpdf <- read.csv(snppositionfile, sep=",", header = T)
## smaple info loading
sampleinfo <- read.csv(sampleinfofile, sep=",", header = T)
row.names(sampleinfo) <- sampleinfo$SID

### sample info group
#### location group
sampleinfo$Location <- sampleinfo$Region_exprosure
attach(sampleinfo)
sampleinfo[which(Location_exprosure != "Guangzhou" & 
                 Country_exprosure == "China"), "Location" ] <- "China" 
sampleinfo[which(Location_exprosure == "Guangzhou"), "Location" ] <- "Local" 
detach(sampleinfo)

##### group list
Locationgroup <- lapply(unique(sampleinfo$Country_exposure), function(x, sampleinfo){
  sampleinfo[which(sampleinfo$Country_exposure== x), "SID"]
}, sampleinfo=sampleinfo) 
names(Locationgroup) <- unique(sampleinfo$Country_exposure)

## load snp type
snpannotation <- read.csv(snpannofile, sep="\t", header = T)
## load tree
tree <- read.tree(treefile)
tree$tip.label <- str_extract(tree$tip.label, pattern="ID[0-9]+")
tree$tip.label[which(is.na(tree$tip.label))] <- "EPI_ISL_402125"
## Adding group info
tree <- groupOTU(tree, Locationgroup, "Location")

### sample names
snames <- colnames(snpdf)[2: ncol(snpdf) ]
Treesnames <- tree$tip.label[ str_detect(tree$tip.label, "ID") ]
sIDs <- str_extract(snames, pattern="ID[0-9]+")
treeIDs <- str_extract(Treesnames, pattern="ID[0-9]+")

# SNP code in the select region
## not_snp = 0; snp = 1; N equas to not_snp;
pindex <- str_c("p.", c(1:nrow(snpdf)))
row.names(snpdf) <- pindex
colnames(snpdf) <- c("EPI_ISL_402125", sIDs)
## Base mutation type
basem_type <- snpdf
BaseMutationType <- function(x){
  x <- as.character(x)
  refbase <- x[1] # refbase
  samep <- which(x==refbase)
  np <- which(x=="N") # Nbase
  MergeBase <- c( which(x == "W"), which(x == "M") )
  mtype <- paste0(refbase, ">", x)
  mtype[c(samep, np, MergeBase) ] <- "0"
  return (mtype)
}
## base mutation
basem_type <- apply(snpdf, 1, BaseMutationType )
basem_type <- as.data.frame( t(basem_type) )
colnames( basem_type ) <- colnames(snpdf)

## transform snp to annotation type
snpdf_type <- snpdf # whole length
XtypeMatch <- function(ix, snpdf, snpannotation){
  xlist <- as.character( snpdf[ix, ] )
  p_anno <- snpannotation[snpannotation$pos==ix, ]
  xtypes <- c()
  # xbase match
  xbaseMatch <- function(xbase, p_anno){
    ptype <- which(p_anno$Alt == xbase) # match alteration
    refbase <- p_anno$Ref[1]
    xtype <- if (xbase == refbase) {
      "0" 
    } else if (xbase == "N") {
       "N"
    } else if (length(ptype) >=1 ) {
      p_anno[ptype, "eff"]
    } else if (xbase == "-" ) {
      "del" 
    } else  { "Unknown"} 
    return(xtype) 
  }
  ## each row
  xtypes <- sapply(xlist, xbaseMatch, p_anno = p_anno)
  return (xtypes)
}
### loop to search the snp type
snpdf_type  <- sapply( c(1:nrow(snpdf)), XtypeMatch, 
                        snpdf=snpdf,
                        snpannotation=snpannotation)
snpdf_type <- as.data.frame( t(snpdf_type) )
rownames(snpdf_type) <- rownames(snpdf)
colnames(snpdf_type) <- colnames(snpdf)

## transform snp value
snpvalueDF <- apply(snpdf_type, 2, function(x){
                    xval <- x
                    xval[which(x != "0")] <- 4 # Other and Unknown N 
                    xval[which(x == "0")] <- 0 # not snp
                    xval[which(x == "Unknown")] <- 4
                    xval[which(str_detect(x, pattern="stream") )] <- 1 # UTR variants
                    xval[which(x == "missense_variant")] <- 2 
                    xval[which(x == "synonymous_variant") ] <- 3 
                    xval[which(x == "del")] <- 5 # - del 
                    xval <- as.numeric(xval)
                    return (xval)})
row.names(snpvalueDF) <- rownames(snpdf_type)
snpvalueDF <- as.data.frame(snpvalueDF)

# ggtree plot -------------------------------------------------------------
## ggtree plot
gtree1 <- ggtree(tree ) + geom_tree() +
          geom_tiplab(aes(color=Local_External, label=""),align=T , linesize = 0.5) #+
gtree2 <- gtree1 + geom_tiplab(aes(label=paste0(label,"/", Location,"/", Lineage) ) , size=1.5)
# gtree with anotation
annodf <- sampleinfo[tree$tip.label,c("Location", "Location") ] 
treehmap <- gheatmap(gtree1, annodf, font.size=3, width=0.1, colnames = F) 
treehmap_lab <- gheatmap(gtree2, annodf, font.size=3, width=0.1, colnames = F) 

### tree data and tree order
treedata <- treehmap$data[which(treehmap$data$isTip), ]
treeIDorders <- treedata[order(treedata$y), "label"]$label

snpdfexport <- t( snpdf[, treeIDorders] )

# Sample Whole Genome SNP plot --------------------------------------------
# Genome SNP dot plot
ReshapMelt <- function(datadf){  
  datadf$position <- as.numeric(  str_remove( rownames(datadf), "p.") )
  
  datadf_melt <- melt(datadf, id.vars=c("position"), 
                      variable.name = "SampleName", value.name = "Value" )
  datadf_melt$MergeID <- paste0(as.character(datadf_melt$SampleName), 
                                "_", 
                                datadf$position)
  return(datadf_melt)
}
allregion <- c(86:29300)
snpvalue_genome <- snpvalueDF[allregion, treeIDorders] ### region select and reorder
snpgenome_melt <- ReshapMelt(snpvalue_genome)
snpgenome_melt  <- snpgenome_melt [, c(1,2,3)]
colnames(snpgenome_melt ) <- c("Position", "SampleName", "Mutation_Effect")
### remove 0
snpgenome_melt <- snpgenome_melt[ which(snpgenome_melt$Mutation_Effect!=0 |
                                        snpgenome_melt$Position %in% c(86,29300) ),  ]
xticks <- c(c(86,29300),c(5654, 8371, 11083, 16846, 22949) ,
            positionlist) %>% sort()
xticks <- c(241, 1059, 2416, 3037, 6312, 8782,12436, 12439, 
            13730, 14408, 15324, 19524, 23403, 23929, 25563, 28881, 28882, 28883,
            28311, 28144) %>% sort()
### plot
genomeplot <- ggplot(snpgenome_melt, aes(Position, SampleName,
                     color = factor(Mutation_Effect)) ) + geom_point(size=1) +
              scale_color_manual(name="Mutation Effect", 
                     labels = c("NonMutation" , "UTR variants", "Missense variant", "Synonymous variant", "N or other" ,"Del"),
                     values = c("0"= "white","1"="#009e73","2"="#d55e00","3"="#56b4e9", "4" = "#D3D3D3" , "5"="#000000") )  +
              labs(x = "", y = "", fontsize=6) + 
              scale_x_continuous(breaks = xticks, labels = xticks ) + 
              scale_y_discrete(breaks = as.character(snpgenome_melt$SampleName) ) + 
              theme(panel.background = element_blank(),
                    panel.grid.major.x = element_line(size=0.5, color="gray"),
                    axis.text.x = element_text(angle=90, size=5),
                    axis.text.y = element_text(size=4),
                    panel.grid.major.y = element_line(size=0.25, color="#D2D3D3"),
                    legend.position = "right")
## histogram
snpgenome_melt$Exposure <- sampleinfo[snpgenome_melt$SampleName, "Exposure_type"]
snpgenome_melt$Exposure_Detailed <- sampleinfo[snpgenome_melt$SampleName, "Exposure_Detailed_classfication"]
snpgenome_meltbar <- snpgenome_melt[ which(snpgenome_melt$Mutation_Effect!=0 &
                                        snpgenome_melt$Mutation_Effect!=4),  ]
gbar <- ggplot(snpgenome_meltbar, aes(Position) ) + 
  geom_bar(aes(fill = factor(Exposure)), width=50) +
  scale_y_continuous(breaks=seq(0,100,50), labels=c(0,50,100) ) +
  labs(x="", y="Count") +
  theme( #panel.background = element_rect(fill = "transparent",colour = "black"),
        panel.grid.major = element_blank(),
        legend.background = element_blank(),
        axis.text.x = element_text(angle=90, size=3) )

## Merge plot and save
genomepdf_merge <- paste0( treefile, "Genome_mutation_merge.pdf" )
gplot <-  genomeplot %>% 
          insert_left(treehmap, width=.5) %>%
          insert_bottom(gbar, height=.2) 
ggsave(genomepdf_merge, gplot, width=20, height=18, useDingbats=FALSE)