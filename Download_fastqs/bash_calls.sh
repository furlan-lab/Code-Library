sFH2
module load R
ml SRA-Toolkit
R
install.packages("colorout")
library(colorout)
#BiocManager::install("SRAdb")
# library(SRAdb)
# sqlfile <- "/fh/fast/furlan_s/grp/lib/SRAmetadb.sqlite"
# sra_con <- dbConnect(SQLite(),sqlfile)
# sra_tables <- dbListTables(sra_con)
# dbListFields(sra_con,"study")
# conversion <- sraConvert( c('SRP048560'), sra_con = sra_con )

##CHP
tab<-"/fh/fast/furlan_s/user/owalt/ewings/SraRunTable.txt"
tabd<-read.csv(tab)
td<-tabd[,] 
td$Run
system('fastq-dump -X 5 -Z SRR11006178')
#works

for(srr in td$Run){
	command<-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='fastq-dump --gzip -O /fh/fast/furlan_s/user/owalt/ewings/sra/sra %s'", srr)
	system(command)
}
cd /fh/fast/furlan_s/user/owalt/ewings/sra/sra


##install chromap
sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='chromap -i -r /fh/fast/furlan_s/grp/refs/GRCh38/GRCh38.p13.genome.fa -o /fh/fast/furlan_s/user/owalt/ewings/chromap/GRCh38.p13_chromap_index'



#awk -F '\t' '{$1 FS "x" FS $2;}1' OFS='\t' GRCh38.p13.chromsizes_main
module load MACS2
module load BEDTools
module load R
module load Kent_tools
LC_COLLATE=C
R
library(colorout)
tab<-"/fh/fast/furlan_s/user/owalt/ewings/SraRunTable.txt"
tabd<-read.csv(tab)
#td<-tabd[tabd$chip_antibody %in% c('FLI1 (Santa Cruz\\, sc-356)', 'WCE'),]
td<-tabd[,]
fqs<-file.path("/fh/fast/furlan_s/user/owalt/ewings/sra-cache/sra/bed", paste0(td$Run, ".fastq.gz"))
file.exists(fqs)


index<-"/fh/fast/furlan_s/user/owalt/ewings/chromap/GRCh38.p13_chromap_index"
fa<-"/fh/fast/furlan_s/grp/refs/GRCh38/GRCh38.p13.genome.fa"
out<-file.path("/fh/fast/furlan_s/user/owalt/ewings/sra-cache/sra/bed", paste0(td$Run, ".bed"))
#align
for(i in 1:length(fqs)){
	command <-sprintf("sbatch -n 1 -c 8 -p campus-new -M gizmo --wrap='chromap -t 8 --preset chip -x %s -r %s -1 %s -o %s'", index, fa, fqs[i], out[i])
	system(command)
}

g1<-td #[grep("SKMNC", td$source_name),]
g1t<-file.path("/fh/fast/furlan_s/user/owalt/ewings/ChIP/FLI1/bed", paste0("SRR14761052_RDES.bed"))
command <- sprintf("sbatch -n 1 -c 8 -p campus-new -M gizmo --wrap='macs2 callpeak -t %s -g hs --extsize 200 -n SRR14761052_RDES --outdir peaks'", g1t)
system(command)


g1name<-g1[,]$Run
#peaks
for(i in 1:length(g1t)){
	command <-sprintf("sbatch -n 1 -c 8 -p campus-new -M gizmo --wrap='macs2 callpeak -t %s -c %s -g hs --extsize 200 -n %s --outdir peaks'", g1t[i], g1name[i])
	system(command)
}
g1<-td[grep("A673", td$source_name),]
g1c<-file.path("bed", paste0(g1[grep("WCE", g1$chip_antibody),]$Run, ".bed"))
g1t<-file.path("bed", paste0(g1[!grepl("WCE", g1$chip_antibody),]$Run, ".bed"))
g1name<-g1[!grepl("WCE", g1$chip_antibody),]$Run
#peaks
for(i in 1:length(g1t)){
	command <-sprintf("sbatch -n 1 -c 8 -p campus-new -M gizmo --wrap='macs2 callpeak -t %s -c %s -g hs --extsize 200 -n %s --outdir peaks'", g1t[i], g1c, g1name[i])
	system(command)
}











beds <- file.path("/fh/fast/furlan_s/user/owalt/ewings/ChIP/FLI1/bed", paste0("SRR14761052_RDES.bed"))
file.exists(beds)
#sort beds
for(i in 1:length(beds)){
	fileout <- gsub("/fh/fast/furlan_s/user/owalt/ewings/ChIP/FLI1/bed/", "/fh/fast/furlan_s/user/owalt/ewings/ChIP/FLI1/bed/sorted/", beds[i])
	command <-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='sort -k 1,1  -k2,2n %s > %s'", beds[i], fileout)
	system(command)
}



for(bed in beds){
	fileout <- gsub("/bed/", "/bed/sorted/", bed)
	command <-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='sort -k 1,1  -k2,2n %s > %s'", bed, fileout)
	system(command)
}

 
mkdir filtered
for FILE in *; do awk '{ if ($1 == "chr1" ||
		   $1 == "chr2" ||
		   $1 == "chr3" ||
		   $1 == "chr4" ||
		   $1 == "chr5" ||
		   $1 == "chr6" ||
		   $1 == "chr7" ||
		   $1 == "chr8" ||
		   $1 == "chr9" ||
		   $1 == "chr10" ||
		   $1 == "chr11" ||
		   $1 == "chr12" ||
		   $1 == "chr13" ||
		   $1 == "chr14" ||
		   $1 == "chr15" ||
		   $1 == "chr16" ||
		   $1 == "chr17" ||
		   $1 == "chr18" ||
		   $1 == "chr19" ||
		   $1 == "chr20" ||
		   $1 == "chr21" ||
		   $1 == "chr22" ||
		   $1 == "chrX" ||
		   $1 == "chrY" ||
		   $1 == "chrM")
		   { print } }' $FILE > filtered/$FILE; done


module load MACS2
module load BEDTools
module load R
module load Kent_tools
LC_COLLATE=C
R
library(colorout)
tab<-"/fh/fast/furlan_s/grp/data/ddata/ewings/SraRunTable.txt"
tabd<-read.csv(tab)
td<-tabd[tabd$chip_antibody %in% c('FLI1 (Santa Cruz\\, sc-356)', 'WCE'),]
fqs<-file.path("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560", paste0(td$Run, ".fastq.gz"))
file.exists(fqs)

beds<-list.files("/fh/fast/furlan_s/user/owalt/ewings/ChIP/FLI1/bed/sorted/filtered/", full.name=T)
beds<- beds[9]
beds<-beds[grep(".bed$", beds)]
chromsizes <- "/fh/fast/furlan_s/grp/refs/chromsizes/GRCh38.p13.chromsizes_main"

#bedgraphs
for(bed in beds){
	outbg <- gsub( "fh/fast/furlan_s/user/owalt/ewings/ChIP/FLI1/bed/sorted/filtered/", "fh/fast/furlan_s/user/owalt/ewings/ChIP/FLI1/bg/", bed)
	outbg <- gsub( ".bed$", ".bedgraph", outbg)
	command <-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='bedtools genomecov -i %s -bg -g %s > %s'", bed, chromsizes, outbg)
	system(command)
}

bgs<-list.files("/fh/fast/furlan_s/user/owalt/ewings/ChIP/FLI1/bg/", full.name=T)
bgs<-bgs[grep(".bedgraph$", bgs)]
bgs<-bgs[7]


# for(bg in bgs){
# 	outbg <- gsub("/bg/", "/bg/sorted/", bg)
# 	command <-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='sort -k1,1  -k2,2n %s > %s'", bg, outbg)
# 	system(command)
# }

# bgs<-list.files("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bed/bg/sorted", full.name=T)
# bgs<-bgs[grep(".bedgraph$", bgs)]

for(bg in bgs){
	outbw <- gsub("/bed/bg/", "/bw/", bg)
	outbw <- gsub(".bedgraph", ".bw", outbw)
	command <-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='bedGraphToBigWig %s %s %s'", bg, chromsizes, outbw)
	system(command)
}



#fix SRR1593966 (SKNMC WCE control)

ml BEDTools
ml Kent_tools
sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='sort -k 1,1  -k2,2n bed/SRR1593966.bed > bed/sorted/SRR1593966.bed'
awk '{ if ($1 == "chr1" ||
		   $1 == "chr2" ||
		   $1 == "chr3" ||
		   $1 == "chr4" ||
		   $1 == "chr5" ||
		   $1 == "chr6" ||
		   $1 == "chr7" ||
		   $1 == "chr8" ||
		   $1 == "chr9" ||
		   $1 == "chr10" ||
		   $1 == "chr11" ||
		   $1 == "chr12" ||
		   $1 == "chr13" ||
		   $1 == "chr14" ||
		   $1 == "chr15" ||
		   $1 == "chr16" ||
		   $1 == "chr17" ||
		   $1 == "chr18" ||
		   $1 == "chr19" ||
		   $1 == "chr20" ||
		   $1 == "chr21" ||
		   $1 == "chr22" ||
		   $1 == "chrX" ||
		   $1 == "chrY" ||
		   $1 == "chrM")
		   { print } }' bed/sorted/SRR1593966.bed > bed/sorted/filtered/SRR1593966.bed

export GENOMES=/fh/fast/furlan_s/grp/refs/chromsizes/GRCh38.p13.chromsizes_main
sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='bedtools genomecov -i bed/sorted/filtered/SRR1593966.bed -bg -g $GENOMES > bed/bg/SRR1593966.bed'
sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='bedGraphToBigWig bed/bg/SRR1593966.bed $GENOMES bw/SRR1593966.bw'


###add labels
ml R
R
library(colorout)
library(tidyr)
tab<-"/fh/fast/furlan_s/grp/data/ddata/ewings/SraRunTable.txt"
tabd<-read.csv(tab)
td<-tabd[tabd$chip_antibody %in% c('FLI1 (Santa Cruz\\, sc-356)', 'WCE'),]


td$label<-paste0(int,"_", td$chip_antibody %>% strsplit(" ") %>% sapply("[[", 1))
td$replace<-paste0(td$Run, "_", td$label)
fin<-list.files("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bw", full.name=T)
ext<-".bw"
fixed<-td$replace[match(gsub(ext, "", basename(fin)), td$Run)]
unfixed<-td$Run[match(gsub(ext, "", basename(fin)), td$Run)]

for(i in 1:length(unfixed)){
	fout<-gsub(unfixed[i], fixed[i], basename(fin[i]))
	#print(paste0("Moving:\n", fin[i], "\nTo:\n", file.path(dirname(fin[i]), fout)))
	file.rename(fin[i], file.path(dirname(fin[i]), fout))
}




tab<-"/fh/fast/furlan_s/user/owalt/ewings/SraRunTable.txt"
tabd<-data.frame(read.csv(tab))

gsm <- read.csv("/fh/fast/furlan_s/user/owalt/ewings/GSM.csv")

td<-tabd
td <-gsub(" ", "", td$source_name)

td$source_name
int<-gsub(" ", "", td$source_name)

df <- cbind(gsm, td$Sample.Name)
write.csv(df, file = "/fh/fast/furlan_s/user/owalt/ewings/temp.csv")
df <- read.csv("/fh/fast/furlan_s/user/owalt/ewings/temp.csv")

condition <- df$Exp

td$condition<- paste0(condition)
td$source_name

write.csv(td, file = "/fh/fast/furlan_s/user/owalt/ewings/test.csv")


td$label<-paste0(td$condition)

td$replace<-paste0(td$Run, "_", td$label)
fin<-list.files("/fh/fast/furlan_s/user/owalt/ewings/sra-cache/sra/bw", full.name=T)
ext<-".bw"

fixed<-td$replace[match(gsub(ext, "", basename(fin)), td$Run)]
unfixed<-td$Run[match(gsub(ext, "", basename(fin)), td$Run)]

for(i in 1:length(unfixed)){
	fout<-gsub(unfixed[i], fixed[i], basename(fin[i]))
	#print(paste0("Moving:\n", fin[i], "\nTo:\n", file.path(dirname(fin[i]), fout)))
	file.rename(fin[i], file.path(dirname(fin[i]), fout))
}











fin<-list.files("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/peaks", full.name=T)
fixed<-td$replace[match(basename(fin) %>% strsplit("_") %>% sapply("[[", 1), td$Run)]
unfixed<-td$Run[match(basename(fin) %>% strsplit("_") %>% sapply("[[", 1), td$Run)]

for(i in 1:length(unfixed)){
	fout<-gsub(unfixed[i], fixed[i], basename(fin[i]))
	#print(paste0("Moving:   ", fin[i], "   To:   ", file.path(dirname(fin[i]), fout)))
	file.rename(fin[i], file.path(dirname(fin[i]), fout))
}

##################################################################################################################################################################
############################################################Histones Narrow#################################################################################
##################################################################################################################################################################
##CHP

cd /fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560

#awk -F '\t' '{$1 FS "x" FS $2;}1' OFS='\t' GRCh38.p13.chromsizes_main
module load MACS2
module load BEDTools
module load R
module load Kent_tools
ml SRA-Toolkit
ml Homer
LC_COLLATE=C
R
library(colorout)
library(tidyr)
tab<-"/fh/fast/furlan_s/grp/data/ddata/ewings/SraRunTable.txt"
tabd<-read.csv(tab)
td<-tabd[grepl("H3K", tabd$chip_antibody),]
fqs<-file.path("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560", paste0(td$Run, ".fastq.gz"))
file.exists(fqs)
index<-"/fh/fast/furlan_s/grp/refs/GRCh38/chromap/GRCh38.p13_chromap_index"
fa<-"/fh/fast/furlan_s/grp/refs/GRCh38/GRCh38.p13.genome.fa"
dir.create("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bed")
dir.create("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bed/sorted")
dir.create("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bed/sorted/filtered")
dir.create("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bed/bg")
dir.create("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bw")
dir.create("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/peaks")
out<-file.path("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bed", paste0(td$Run, ".bed"))
#download
for(srr in td$Run){
	command<-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='fastq-dump --gzip -O /fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560 %s'", srr)
	system(command)
}



#align
for(i in 1:length(fqs)){
	command <-sprintf("sbatch -n 1 -c 8 -p campus-new -M gizmo --wrap='chromap -t 8 --preset chip -x %s -r %s -1 %s -o %s'", index, fa, fqs[i], out[i])
	system(command)
}


#make control tags
g1c<-"/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/FLI1_done/bed/SRR1593966.bed"
dir.create("tags")
for(i in 1:length(g1c)){
	run<-gsub(".bed", "", basename(g1c[i]))
	command<-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='makeTagDirectory tags/%s %s'", run, g1c[i])
	system(command)
}
#make treatment tags
g1<-td[grep("SKMNC", td$source_name),]
g1t<-file.path("bed", paste0(g1[!grepl("WCE", g1$chip_antibody),]$Run, ".bed"))
for(i in 1:length(g1t)){
	run<-gsub(".bed", "", basename(g1t[i]))
	frun<-gsub(".bed/", "", run)
	command<-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='makeTagDirectory tags/%s %s'", frun, g1t[i])
	system(command)
}

#make control tags G2
g1c<-"/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/FLI1_done/bed/SRR1593991.bed"
for(i in 1:length(g1c)){
	run<-gsub(".bed", "", basename(g1c[i]))
	command<-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='makeTagDirectory tags/%s %s'", run, g1c[i])
	system(command)
}
#make treatment tags G2
g1<-td[grep("A673", td$source_name),]
g1t<-file.path("bed", paste0(g1[!grepl("WCE", g1$chip_antibody),]$Run, ".bed"))
for(i in 1:length(g1t)){
	run<-gsub(".bed", "", basename(g1t[i]))
	frun<-gsub(".bed/", "", run)
	command<-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='makeTagDirectory tags/%s %s'", frun, g1t[i])
	system(command)
}


#peaks G1
g1<-td[grep("SKMNC", td$source_name),]
g1t<-file.path("tags", g1[!grepl("WCE", g1$chip_antibody),]$Run)
g1c<-"tags/SRR1593966"
for(i in 1:length(g1t)){
	command <-sprintf("sbatch -n 1 -c 8 -p campus-new -M gizmo --wrap='findPeaks %s -style histone -o auto -i %s'", g1t[i], g1c)
	system(command)
}
#peaks G2
g1<-td[grep("A673", td$source_name),]
g1t<-file.path("tags", g1[!grepl("WCE", g1$chip_antibody),]$Run)
g1c<-"tags/SRR1593991"
for(i in 1:length(g1t)){
	command <-sprintf("sbatch -n 1 -c 8 -p campus-new -M gizmo --wrap='findPeaks %s -style histone -o auto -i %s'", g1t[i], g1c)
	system(command)
}



beds<-list.files("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bed", full.name=T)
beds<-beds[grep(".bed$", beds)]
chromsizes <- "/fh/fast/furlan_s/grp/refs/chromsizes/GRCh38.p13.chromsizes_main"

#sort beds
for(bed in beds){
	fileout <- gsub("/bed/", "/bed/sorted/", bed)
	command <-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='sort -k 1,1  -k2,2n %s > %s'", bed, fileout)
	system(command)
}

cd bed/sorted
for FILE in *; do awk '{ if ($1 == "chr1" ||
		   $1 == "chr2" ||
		   $1 == "chr3" ||
		   $1 == "chr4" ||
		   $1 == "chr5" ||
		   $1 == "chr6" ||
		   $1 == "chr7" ||
		   $1 == "chr8" ||
		   $1 == "chr9" ||
		   $1 == "chr10" ||
		   $1 == "chr11" ||
		   $1 == "chr12" ||
		   $1 == "chr13" ||
		   $1 == "chr14" ||
		   $1 == "chr15" ||
		   $1 == "chr16" ||
		   $1 == "chr17" ||
		   $1 == "chr18" ||
		   $1 == "chr19" ||
		   $1 == "chr20" ||
		   $1 == "chr21" ||
		   $1 == "chr22" ||
		   $1 == "chrX" ||
		   $1 == "chrY" ||
		   $1 == "chrM")
		   { print } }' $FILE > filtered/$FILE; done


module load MACS2
module load BEDTools
module load R
module load Kent_tools
LC_COLLATE=C
R
library(colorout)
tab<-"/fh/fast/furlan_s/grp/data/ddata/ewings/SraRunTable.txt"
tabd<-read.csv(tab)
td<-tabd[grepl("H3K", tabd$chip_antibody),]
fqs<-file.path("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560", paste0(td$Run, ".fastq.gz"))
file.exists(fqs)

beds<-list.files("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bed/sorted/filtered", full.name=T)
beds<-beds[grep(".bed$", beds)]
chromsizes <- "/fh/fast/furlan_s/grp/refs/chromsizes/GRCh38.p13.chromsizes_main"

#bedgraphs
for(bed in beds){
	outbg <- gsub( "/bed/sorted/filtered/", "/bed/bg/", bed)
	outbg <- gsub( ".bed$", ".bedgraph", outbg)
	command <-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='bedtools genomecov -i %s -bg -g %s > %s'", bed, chromsizes, outbg)
	system(command)
}

bgs<-list.files("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bed/bg", full.name=T)
bgs<-bgs[grep(".bedgraph$", bgs)]


for(bg in bgs){
	outbw <- gsub("/bed/bg/", "/bw/", bg)
	outbw <- gsub(".bedgraph", ".bw", outbw)
	command <-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='bedGraphToBigWig %s %s %s'", bg, chromsizes, outbw)
	system(command)
}

###add labels
int<-gsub(" ", "", td$source_name)

td$label<-paste0(int,"_", td$chip_antibody %>% strsplit(" ") %>% sapply("[[", 1))
td$replace<-paste0(td$Run, "_", td$label)
fin<-list.files("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/bw", full.name=T)
ext<-".bw"
fixed<-td$replace[match(gsub(ext, "", basename(fin)), td$Run)]
unfixed<-td$Run[match(gsub(ext, "", basename(fin)), td$Run)]

for(i in 1:length(unfixed)){
	fout<-gsub(unfixed[i], fixed[i], basename(fin[i]))
	#print(paste0("Moving:\n", fin[i], "\nTo:\n", file.path(dirname(fin[i]), fout)))
	file.rename(fin[i], file.path(dirname(fin[i]), fout))
}

fin<-list.files("/fh/fast/furlan_s/grp/data/ddata/ewings/SRP048560/tags", full.name=T)
tf<-basename(fin) %in% td$Run
fin<-fin[tf]
unfixed<-file.path(fin, "regions.txt")
fixed<-file.path(fin, paste0(td$replace[match(basename(fin), td$Run)], "_regions.txt"))
dir.create("regions")
for(i in 1:length(unfixed)){
	#print(paste0("Moving:   ", unfixed[i], "   To:   ", fixed[i]))
	file.rename(unfixed[i], fixed[i])
	file.copy(fixed[i], file.path("regions", basename(fixed[i])))
}































##HiC
cd /fh/fast/furlan_s/grp/data/ddata/ewings/HiC/fastq
sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='wget https://www.encodeproject.org/files/ENCFF205PEP/@@download/ENCFF205PEP.fastq.gz'
sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='wget https://www.encodeproject.org/files/ENCFF940BBC/@@download/ENCFF940BBC.fastq.gz'
sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='wget https://www.encodeproject.org/files/ENCFF253NHH/@@download/ENCFF253NHH.fastq.gz'
sbatch -n 1 -c 1 -p campus-new -M gizmo --wrap='wget https://www.encodeproject.org/files/ENCFF811IQJ/@@download/ENCFF811IQJ.fastq.gz'


# cd /fh/fast/furlan_s/grp/data/ddata/ewings/HiC
# #chromap --preset hic -x index -r ref.fa -1 fastq/ENCFF205PEP.fastq.gz -2 fastq/ENCFF940BBC.fastq.gz -o aln.pairs
# ml Python/3.8.6-GCCcore-10.2.0
# cd
# python -m venv ~/.virtualenvs/HiC-Pro
# source ~/.virtualenvs/HiC-Pro/bin/activate
# pip install --upgrade pip
# pip install --upgrade numpy
# pip install bx
# pip install bx_python
# pip install jsonschema
# pip install scipy
# pip install pysam
# pip install iced
# pip install pandas
# cd ~/software/HiC-Pro

# ml R
# ml SAMtools
# ml Bowtie2
# make configure

cd /fh/fast/furlan_s/grp/data/ddata/ewings/HiC
mkdir pairs
mkdir pairs/SKNMC1
mkdir pairs/SKNMC2
sbatch -n 1 -c 16 -p campus-new -M gizmo --wrap='chromap --preset hic -x /fh/fast/furlan_s/grp/refs/GRCh38/chromap/GRCh38.p13_chromap_index -t 16 -r /fh/fast/furlan_s/grp/refs/GRCh38/GRCh38.p13.genome.fa -1 fastq/ENCFF205PEP.fastq.gz -2 fastq/ENCFF940BBC.fastq.gz -o pairs/SKNMC1/aln.pairs'           # Hi-C reads and pairs output
sbatch -n 1 -c 16 -p campus-new -M gizmo --wrap='chromap --preset hic -x /fh/fast/furlan_s/grp/refs/GRCh38/chromap/GRCh38.p13_chromap_index -t 16 -r /fh/fast/furlan_s/grp/refs/GRCh38/GRCh38.p13.genome.fa -1 fastq/ENCFF253NHH.fastq.gz -2 fastq/ENCFF811IQJ.fastq.gz -o pairs/SKNMC2/aln.pairs'           # Hi-C reads and pairs output



