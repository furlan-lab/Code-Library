#run from here to...

install.packages('reticulate')
reticulate::py_config()
reticulate::py_available("scanpy")
install.packages('IRkernel')
install.packages(c("devtools", "rlist","heatmaply", "hdf5r", "tidyverse", 
                   "openxlsx", "prettydoc", "kableExtra", "piano", "pals", 
                   "ClusterR", "mixtools","Cairo", "refGenome", "ggplot2"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
packages<-c('BiocGenerics','DelayedArray','DelayedMatrixStats','limma',
            'S4Vectors','SingleCellExperiment','SummarizedExperiment',
            'batchelor','apeglm','EnhancedVolcano','ashr','fgsea',
            'Seurat','TxDb.Hsapiens.UCSC.hg38.knownGene','multtest',
            'edgeR','simpleaffy','sva','made4','RBGL','reshape','WGCNA','GEOquery',
            'oligo','GSEABase','Rsubread','refGenome','statmod','genbankr','Gviz',
            'ggbio','vcfR','AnnotationHub','EnrichedHeatmap','GenomicRanges','motifmatchr',
            'chromVAR','Rsamtools','Biostrings','ComplexHeatmap','ShortRead','BSgenome',
            'EnsDb.Hsapiens.v86','BSgenome.Hsapiens.UCSC.hg38','JASPAR2018',
            'patchwork','R453Plus1Toolbox','rgr','GenomicFeatures','BiocStyle','HPAanalyze')
BiocManager::install(packages)
packageurl <- "https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.9.900.3.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source") 
BiocManager::install(c('DESeq2','TFBSTools'))
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3', ref = "develop")
devtools::install_github('scfurl/m3addon')
devtools::install_github('scfurl/ArchR')
ArchR::installExtraPackages()
BiocManager::install("rhdf5")
devtools::install_github("timoast/signac", ref = "develop")
devtools::install_github('scfurl/seqGlue')
install.packages("Vennerable", repos="https://R-Forge.R-project.org", type="source")
devtools::install_github('scfurl/probedeeper')
devtools::install_github('ncborcherding/scRepertoire')
packages<-c('ChIPpeakAnno')
BiocManager::install(packages)
packages<-c('org.Mm.eg.db')
BiocManager::install(packages) 
packages<-c('factoextra')
install.packages(packages)
packages<-c('ChIPQC')
BiocManager::install(packages) 
packages<-c('rGREAT')
BiocManager::install(packages) 
packages<-c('soGGi')
BiocManager::install(packages) 
packages<-c('DiffBind')
BiocManager::install(packages) 
BiocManager::install("piano")
# ...to here!!



###If some do not; install individual script download below // also within the 'Tools' tab you can manually install (from CRAN)

install.packages('Seurat')
library(Seurat)

# Install the remotes package // dev version
install.packages('remotes')
remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)

##Seurat current release
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Seurat")
library(Seurat)
##Signac dev version
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("timoast/signac", ref = "develop", force = T)
library(Signac)

##devtools install; basic needed
devtools::check()
devtools::install_github("r-lib/devtools")

##scCustomize install
devtools::install_github(repo = "samuel-marsh/scCustomize")
remotes::install_github(repo = "samuel-marsh/scCustomize")
library(scCustomize)

##install.packages('ggrastr')//dev version
install.packages('devtools')
devtools::install_github('VPetukhov/ggrastr', build_vignettes = TRUE)
library("devtools") 

##required to run "h_cols <-rev(brewer.pal(name = "RdYlBu", n = 7))"
library(RColorBrewer)

##BiocManager install
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BSgenome.Hsapiens.UCSC.hg38 install
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
devtools::install_github("crazyhottommy/scATACutils")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("JASPAR2020")                     
                     
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnsDb.Hsapiens.v86")

devtools::install_github(repo = "samuel-marsh/scCustomize", ref = "develop")

remotes::install_github(repo = "samuel-marsh/scCustomize", ref = "develop") 
##restart and run...
install.packages("systemfonts", type = "source")
library(systemfonts)

install.packages('ggrastr')
library(ggrastr)
install.packages('devtools')
devtools::install_github('VPetukhov/ggrastr', build_vignettes = TRUE)

install.packages('nloptr')
library(nloptr)
devtools::install_github("r-lib/textshaping")
library(textshaping)

library(ArchR)
devtools::install_github("slowkow/ggrepel")
library(ggrepel)

#pbmcapply install for GSEA of cluster markers chunk
install.packages('pbmcapply')

install.packages('qlcMatrix')

devtools::install_github("hypercompetent/colorway", force = T)
library(colorway)

#FigR install (TF Motif analysis)
install.packages("remotes")
remotes::install_github("buenrostrolab/FigR")


#fun colors
install.packages("viridis")
library(viridis)
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()
devtools::install_github("slowkow/ggrepel")
library(ggrepel)


#install for heatmaps
library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
#CRAN
install.packages("circlize")
#Github
devtools::install_github("jokergoo/circlize")
library(circlize)
BiocManager::install("EnrichedHeatmap")
library(EnrichedHeatmap)




##IDK not req
#unlink("~/.R/Makevars")
#unlink("~/.Renviron")
#PATH="/usr/local/clang7/bin:${PATH}"
##System home set up; directory issue
#setwd("~")
#Sys.getenv()
#Sys.setenv(R_LIBS_SITE = "/home/cmartin3/R/x86_64-pc-linux-gnu-library/4.1")