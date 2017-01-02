
# -------------------------------------
# specify paths and load functions
# -------------------------------------
DATA_DIR <- paste(ROOT_DIR, "/Analysis/MLR2/WS", sep="")      # SPECIFY HERE
PROG_DIR <- paste(ROOT_DIR, "/Analysis/MLR2/prog", sep="")      # SPECIFY HERE
RES_DIR  <- paste(ROOT_DIR, "/Analysis/MLR2/res/", key, sep="")      # SPECIFY HERE
source(file.path(PROG_DIR,'SFfunc.R'))
dir.create(file.path(RES_DIR), showWarnings=F)
source(file.path(PROG_DIR,'SFplotting.R'))
source(file.path(PROG_DIR,'support_functions.R'))
source(file.path(PROG_DIR,'color.R'))
colorfile<-file.path(PROG_DIR, "CBSafe15.csv")
colors<-readColorFile(colorfile)
colors<-as.character(colors)
#knitr::opts_knit$set(root.dir = ROOT_DIR)
library(colorout)
library(Matrix)
library(monocle)
library(stringr)
library(slam)
library(pheatmap)
library(matrixStats)
#library(plyr)
library(dplyr)
library(reshape2)
library(piano)
library(DDRTree)
library(gridExtra)
library(XLConnect)
library(tsne)
library(Rtsne)
library(e1071)
library(RColorBrewer)
#library(densityClust)
library(devtools)
load_all(file.path(ROOT_DIR, "monocle-dev"))
load_all(file.path(ROOT_DIR, "densityClust"))
#load_all(file.path(ROOT_DIR, "fstree"))
# Set global ggplot2 properties for making print-scaled PDF panels
SFtheme<-theme_bw(base_size=14) + 
	theme(panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA))
theme_set(SFtheme)



set.seed(0)
print("done")
mix<-readRDS(file.path(DATA_DIR, "CDS2_Final.RDS"))
genes2<-read.table(file.path(DATA_DIR, "genes2.tsv"))

# Filter super high mRNA cells, which are probably not singletons:
mRNA_thresh <- 10000
removedhigh<-mix[,which(pData(mix)$Total_mRNAs>mRNA_thresh)]
mix <- mix[,pData(mix)$Total_mRNAs < 10000]

# Remove low mRNA cells: (There are none)
keep <- detectGenes(mix, min_expr = 0.1)
which(!rownames(pData(mix)) %in% rownames(pData(keep)))
mix<-keep
rm(keep)
cthSVM <- newCellTypeHierarchy()
cthSVM <- addCellType(cthSVM, "CD3s", classify_func=function(x) {x["CD3D",] > 0})
cthSVM <- addCellType(cthSVM, "CD4s", classify_func=function(x) {x["CD4",] > 0}, parent_cell_type_name = "CD3s")
cthSVM <- addCellType(cthSVM, "CD8s", classify_func=function(x) {x["CD8A",] > 0 | x["CD8B",] > 0 }, parent_cell_type_name = "CD3s")
cthSVM <- addCellType(cthSVM, "Bcells", classify_func=function(x) 
  {x["MS4A1",] > 0})
cthSVM <- addCellType(cthSVM, "Monos", classify_func=function(x) 
  {x["CD14",] > 0  })
cthSVM <- addCellType(cthSVM, "NKs", classify_func=function(x) 
  {x["KLRD1",] > 0 |
  x["NCAM1",] > 0})
