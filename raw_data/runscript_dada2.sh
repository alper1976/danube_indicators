### Use freebio to run dada2

## Illumina run 160408_M00485_0264_000000000-AME7C

# use fastx_toolkit for jds2 data

# transfer files to scratch

mkdir /work/users/alexaei/input
mkdir /work/users/alexaei/output

rsync -avc --progress 01_raw_data/00_danube/01_jds3/160408_M00485_0264_000000000-AME7C/*/*.fastq.gz /work/users/alexaei/input/.

mkdir /work/users/alexaei/databases
mkdir /work/users/alexaei/databases/silva
cd /work/users/alexaei/databases/silva
wget https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz

##############################################
### run cutadapt to remove primer sequences
##############################################
# START TMUX
tmux

#start processes inside the started tmux session

DATADIR=/work/users/alexaei #Root of all data
DIRS=$DATADIR/input              #All directories in which to look for data files
mkdir -p $DATADIR/output/AdaptersRemoved       #Creating directory for new fastq files

for f in $DIRS/*.fastq.gz
    do
        FILE=${f#$DIRS}         #Removing directory form file name
        if [[ $FILE =~ R1 ]] #If file name contains R1
        then
            cutadapt -g NNNNCCTACGGGNGGCWGCAG -o $DATADIR/output/AdaptersRemoved/$FILE $f    #Removing second half of primer 574*f
        elif [[ $FILE =~ R2 ]]
        then
            cutadapt -g GACTACHVGGGTATCTAATCC -o $DATADIR/output/AdaptersRemoved/$FILE $f       #Removing second half of primer 1132r
        fi
    done

# detach from tmux session by typing Ctrl + b and then d

# if you want multiple sessions running in parallel you should name the session by Ctrl + b and $

# one can list all running sessions with tmux list-session

# move files to abel for further processing with dada2
rsync -avc --progress /work/users/alexaei/output/* /usit/abel/u1/alexaei/nobackup/01_data/00_danube/160408_M00485_0264_000000000-AME7C/.
rm -r /work/users/alexaei/input/*
rm -r /work/users/alexaei/output/*

# move files to abel work folder for further processing with dada2
cd /usit/abel/u1/alexaei/nobackup/01_data/00_danube/160408_M00485_0264_000000000-AME7C
mkdir /work/users/alexaei/output/160408_M00485_0264_000000000-AME7C
rsync -avc --progress /usit/abel/u1/alexaei/nobackup/01_data/00_danube/160408_M00485_0264_000000000-AME7C/AdaptersRemoved /work/users/alexaei/output/160408_M00485_0264_000000000-AME7
## run dada2 interactive in R
#qlogin --account=uio --ntasks=1
#
#module load R
#
#R

# or submit to slurm queue

sbatch /usit/abel/u1/alexaei/scripts/projects/00_danube/160408_M00485_0264_000000000-AME7C/run_dada2_abel.slurm

mv /work/users/alexaei/output/*.pdf /usit/abel/u1/alexaei/nobackup/01_data/00_danube/160408_M00485_0264_000000000-AME7C/.
mv /work/users/alexaei/output/*.rds /usit/abel/u1/alexaei/nobackup/01_data/00_danube/160408_M00485_0264_000000000-AME7C/.

###################################################################################################################################
# second batch
###################################################################################################################################

## Illumina run 161010_M00485_0309_000000000-ATHUN


rsync -avc --progress 01_raw_data/00_danube/01_jds3/161010_M00485_0309_000000000-ATHUN/*/*.fastq.gz /work/users/alexaei/input/.

##############################################
### run cutadapt to remove primer sequences
##############################################
# START TMUX
tmux

#start processes insiede the started tmux session

DATADIR=/work/users/alexaei #Root of all data
DIRS=$DATADIR/input              #All directories in which to look for data files
mkdir -p $DATADIR/output/AdaptersRemoved       #Creating directory for new fastq files

for f in $DIRS/*.fastq.gz
    do
        FILE=${f#$DIRS}         #Removing directory form file name
        if [[ $FILE =~ R1 ]] #If file name contains R1
        then
            cutadapt -g NNNNCCTACGGGNGGCWGCAG -o $DATADIR/output/AdaptersRemoved/$FILE $f    #Removing second half of primer 574*f
        elif [[ $FILE =~ R2 ]]
        then
            cutadapt -g GACTACHVGGGTATCTAATCC -o $DATADIR/output/AdaptersRemoved/$FILE $f       #Removing second half of primer 1132r
        fi
    done

mkdir /usit/abel/u1/alexaei/nobackup/01_data/00_danube/161010_M00485_0309_000000000-ATHUN
rsync -avc --progress /work/users/alexaei/output/* /usit/abel/u1/alexaei/nobackup/01_data/00_danube/161010_M00485_0309_000000000-ATHUN/.
rm -r /work/users/alexaei/input/*
rm -r /work/users/alexaei/output/*

# go to abel
mkdir /work/users/alexaei/output/161010_M00485_0309_000000000-ATHUN
rsync -avc --progress /usit/abel/u1/alexaei/nobackup/01_data/00_danube/161010_M00485_0309_000000000-ATHUN/AdaptersRemoved /work/users/alexaei/output/161010_M00485_0309_000000000-ATHUN/.

# fix files

sbatch /usit/abel/u1/alexaei/scripts/projects/00_danube/161010_M00485_0309_000000000-ATHUN/run_dada2_abel.slurm

mv /work/users/alexaei/output/*.pdf /usit/abel/u1/alexaei/nobackup/01_data/00_danube/161010_M00485_0309_000000000-ATHUN/.
mv /work/users/alexaei/output/*.rds /usit/abel/u1/alexaei/nobackup/01_data/00_danube/161010_M00485_0309_000000000-ATHUN/.

###################################################################################################################################
# jds2 131002_M00485_0074_000000000-A5GKR
###################################################################################################################################

## Illumina run jds2

mkdir /usit/abel/u1/alexaei/nobackup/01_data/00_danube/jds2
cd /usit/abel/u1/alexaei/nobackup/01_data/00_danube/jds2

### download files from SRA

wget 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&term=PRJNA256993' -O - | grep -v '^Run' |cut -d',' -f1 | xargs fastq-dump --gzip --split-3

### process with cutadapt - freebio

mkdir /work/users/alexaei/output/jds2
rsync -avc --progress /usit/abel/u1/alexaei/nobackup/01_data/00_danube/jds2 /work/users/alexaei/output/jds2/.

# START TMUX
tmux

#start processes insiede the started tmux session

DATADIR=/work/users/alexaei/output #Root of all data
DIRS=$DATADIR/jds2              #All directories in which to look for data files
mkdir -p $DATADIR/jds2/AdaptersRemoved     #Creating directory for new fastq files

for f in $DIRS/*.fastq.gz
  do
    FILE=${f#$DIRS}         #Removing directory form file name
    if [[ $FILE =~ _1. ]] #If file name contains R1
    then
        cutadapt -g NNNNCCTACGGGNGGCWGCAG -o $DATADIR/jds2/AdaptersRemoved/$FILE $f    #Removing second half of primer 574*f
    elif [[ $FILE =~ _2. ]]
    then
        cutadapt -g GACTACHVGGGTATCTAATCC -o $DATADIR/jds2/AdaptersRemoved/$FILE $f       #Removing second half of primer 1132r
    fi
  done

rsync -avc --progress $DATADIR/jds2/AdaptersRemoved /usit/abel/u1/alexaei/nobackup/01_data/00_danube/jds2/.

# move files to abel work folder for further processing with dada2
mkdir /work/users/alexaei/output/jds2
rsync -avc --progress /usit/abel/u1/alexaei/nobackup/01_data/00_danube/jds2/AdaptersRemoved /work/users/alexaei/output/jds2/.

module load python3

cd /work/users/alexaei/output/jds2/AdaptersRemoved #Root of all data

for f in *_1.fastq.gz
    do
        RFILE=${f%_1.fastq.gz}_2.fastq.gz
        gunzip $f
        gunzip $RFILE
        python3 /usit/abel/u1/alexaei/scripts/paired_end_pairer.py -l ${f%.gz} -r ${RFILE%.gz}
        gzip ${f%.gz}
        gzip ${RFILE%.gz}
    done

mkdir cd /work/users/alexaei/output/jds2/pairs
mv /work/users/alexaei/output/jds2/AdaptersRemoved/*paired* /work/users/alexaei/output/jds2/pairs/.

find /work/users/alexaei/output/jds2/pairs -maxdepth 1 -type f -size 0 -delete

cd /work/users/alexaei/output/jds2/pairs
for i in *paired*; do gzip $i; done

# two files

cd /usit/abel/u1/alexaei/scripts/projects/00_danube/jds2

sbatch run_dada2_abel.slurm




