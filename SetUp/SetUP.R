#run from here to...

install.packages('reticulate')
reticulate::py_config()
reticulate::py_available("scanpy")
install.packages('IRkernel')
install.packages('remotes')
install.packages(c("devtools", "rlist","heatmaply", "hdf5r", "tidyverse", 
                   "openxlsx", "prettydoc", "kableExtra", "piano", "pals", 
                   "ClusterR", "mixtools","Cairo", "refGenome", "ggplot2", "RColorBrewer", "systemfonts", 'nloptr', 'pbmcapply', 'qlcMatrix'
                  'parallel', 'stringr', 'viridis', 'tibble', 'dplyr', 'future', 'ggplot2', 'Matrix', 'xfun',"ggpubr", "hdf5r", "rliger"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
packages<-c('BiocGenerics','DelayedArray','DelayedMatrixStats','limma',
            'S4Vectors','SingleCellExperiment','SummarizedExperiment',
            'batchelor','apeglm','EnhancedVolcano','ashr','fgsea','lme4',
            ,'TxDb.Hsapiens.UCSC.hg38.knownGene','multtest','HDF5Array',
            'terra', 'ggrastr','edgeR','simpleaffy','sva','made4','RBGL','reshape','WGCNA','GEOquery',
            'oligo','GSEABase','Rsubread','refGenome','statmod','genbankr','Gviz',
            'ggbio','vcfR','AnnotationHub','EnrichedHeatmap','GenomicRanges','motifmatchr',
            'chromVAR','Rsamtools','Biostrings','ComplexHeatmap','ShortRead','BSgenome',
            'EnsDb.Hsapiens.v86','BSgenome.Hsapiens.UCSC.hg38','JASPAR2018',
            'patchwork','R453Plus1Toolbox','rgr','GenomicFeatures','BiocStyle','HPAanalyze','DESeq2','TFBSTools',
           'rhdf5','ChIPpeakAnno','org.Mm.eg.db' ,'factoextra', 'msigdbr', 'JASPAR2020', "EnrichedHeatmap","dittoSeq", 
            "DropletUtils", "Nebulosa")
BiocManager::install(packages)
packageurl <- "https://cran.r-project.org/src/contrib/Archive/RcppArmadillo/RcppArmadillo_0.9.900.3.0.tar.gz"
install.packages(packageurl, repos=NULL, type="source") 

devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3', ref = "develop")
devtools::install_github('scfurl/m3addon')
devtools::install_github('GreenleafLab/ArchR')
ArchR::installExtraPackages()
devtools::install_github("timoast/signac", ref = "develop")
devtools::install_github('scfurl/seqGlue')
devtools::install_github('scfurl/probedeeper')
devtools::install_github('ncborcherding/scRepertoire')

devtools::install_github(repo = "samuel-marsh/scCustomize")
devtools::install_github('VPetukhov/ggrastr', build_vignettes = TRUE)
devtools::install_github("crazyhottommy/scATACutils")
devtools::install_github('VPetukhov/ggrastr', build_vignettes = TRUE)
devtools::install_github("r-lib/textshaping")
devtools::install_github("slowkow/ggrepel")
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("jokergoo/circlize")
devtools::install_github("enblacar/SCpubr", ref = "v1.1.1-dev-stable")
# Install the remotes package // dev version
remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')

##Signac dev version
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("timoast/signac", ref = "develop", force = T)

#FigR install (TF Motif analysis)
remotes::install_github("buenrostrolab/FigR")
