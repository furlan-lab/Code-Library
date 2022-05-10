wd=/fh/scratch/delete90/furlan_s/owalt/bulk_HSCT
mkdir -p $wd
cd $wd
mkdir fastq

ml R/4.1.0-foss-2020b
ml SRA-Toolkit
module load cutadapt
module load TrimGalore
module load FastQC
ml Python
R

library(data.table)
DIR<-"/fh/scratch/delete90/furlan_s/owalt/bulk_HSCT"

bch<-fread(file.path(DIR, "meta.txt"))
outdir<-file.path(DIR, "fastq")
setwd(outdir)
cmd<-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --mem-per-cpu=4000MB --wrap='fastq-dump --skip-technical --split-3 --origfmt --gzip %s'", bch$Run)
#system(cmd[2])
ls
fq<-c(paste0(bch$Run, "_1.fastq.gz"), paste0(bch$Run, "_2.fastq.gz"))
file.exists(fq)

fq="/fh/scratch/delete90/furlan_s/owalt/bulk_HSCT/fastq"
cd $fq
ml Python
piperna GENOMESFILE
piperna MAKERUNSHEET -gk GRCh38 -fq $fq -r1c _1 -r2c _2 -t pe -sc _ -b file
piperna ALIGN -r runsheet.csv 
cat piperna.log
piperna SUMMARIZE -r runsheet.csv -th 24 -m 16 -l pipernaS.log

#worked!!!

#PTCy bulk data
#---------------------------------------------------------------------------------------#

wd=/fh/scratch/delete90/furlan_s/owalt/PTCy
cd $wd
mkdir fastq

ml R/4.1.0-foss-2020b
ml SRA-Toolkit

R

library(data.table)
DIR<-"/fh/scratch/delete90/furlan_s/owalt/PTCy"

bch<-fread(file.path(DIR, "meta.txt"))
outdir<-file.path(DIR, "fastq")
setwd(outdir)
cmd<-sprintf("sbatch -n 1 -c 1 -p campus-new -M gizmo --mem-per-cpu=4000MB --wrap='fastq-dump --skip-technical --split-3 --origfmt --gzip %s'", bch$Run)
#system(cmd[2])
for(cm in cmd){
	system(cm)
}

fq<-c(paste0(bch$Run, "_1.fastq.gz"), paste0(bch$Run, "_2.fastq.gz"))
file.exists(fq)

fq="/fh/scratch/delete90/furlan_s/owalt/PTCy/fastq"
cd $fq
ls
ml Python

python3 -m pip install --user pipx
python3 -m pipx ensurepath
pipx install --include-deps piperna

pip install piperna



piperna GENOMESFILE
piperna MAKERUNSHEET -gk GRCh38 -fq $fq -r1c _1 -r2c _2 -t pe -sc _ -b file
piperna ALIGN -r runsheet.csv 
cat piperna.log
piperna SUMMARIZE -r runsheet.csv -th 24 -m 16 -l pipernaS.log



