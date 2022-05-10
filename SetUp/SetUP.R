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
