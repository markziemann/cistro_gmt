#!/bin/bash
#set -x

###############################################################
#
#  Curate TFBS gene sets from cistrome ENCODE data
#  by Mark Ziemann and Antony Kaspi
#  2017-2018
#  This software is distributed with and may be used according
#  to 2-clause FreeBSD License
#
###############################################################

# Description: Will download TFBS peak sets and will
# annotate TSS and enhancers bound by transcription factors
# and create a gene matrix that can be used in GSEA and other
# pathway analysis tools

###############################################################
echo "Testing dependancies"
###############################################################

PARALLEL=`which parallel | wc -l`
if [ $PARALLEL -eq "0" ] ; then
  echo 'Error: gnu parallel was not found.
It can be downloaded from the following URL:
https://www.gnu.org/software/parallel/
Also available from the Ubuntu software centre:
sudo apt-get install parallel'
  exit
fi

LIFTOVER=`which liftOver | wc -l`
if [ $LIFTOVER -eq "0" ] ; then
  echo 'Error: liftOver was not found.
It can be downloaded from the following URL:
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver'
  exit
fi

BEDTOOLS=`which bedtools | wc -l`
if [ $BEDTOOLS -eq "0" ] ; then
  echo 'Error: bedtools was not found.
It can be downloaded from the following URL:
http://bedtools.readthedocs.io/en/latest/
Older version is available from the Ubuntu software centre:
sudo apt-get install bedtools'
  exit
fi

###############################################################
# Lets declare a bunch of vars for later
###############################################################
GTFURL=ftp://ftp.ensembl.org/pub/release-86/gtf/mus_musculus/Mus_musculus.GRCm38.86.gtf.gz
GTF=`basename $GTFURL`
TSSBED=`echo $GTF | sed 's/.gtf.gz/.tss.bed/'`
SIZE_RANGE='2500'
GNAMES=$(basename $GTF .gtf.gz).gnames.txt
#this encode stuff will have to be modified
METADATA="mouse_factor_full_QC.txt"
METADATA_SUMMARY=metadata_summary.tsv
URLLIST=files.txt

NUMRANGE='500'
UPPER_LIMIT='5000'
NRCPU=`nproc`

###############################################################
echo "Dependancies OK. Downloading required annotation files"
###############################################################
wget -N $GTFURL

###############################################################
echo "extracting gene names from GTF"
###############################################################
zcat $GTF | grep -w gene | cut -f9 | cut -d '"' -f2,6 \
| tr '"' '\t' | sort -uk 2b,2 > $GNAMES

###############################################################
echo "mm10 Enhancers obtained from Enhancerdb"
###############################################################
#Get coordinates of plus strand genes
ENH=enhancers.bed

###############################################################
echo "intersect peaks with TSS"
###############################################################

map_pk(){
BEDIN=$1
CODENUM=$(echo $(basename $BEDIN) | cut -d '_' -f1)
METADATA=$2
TSS=$3
DIST=$4
MAX_GENES=$5

OUT=$(echo $BEDIN | sed 's/narrowPeak.bed/narrowPeak.gmt/')
NAME=$(grep -w ^$CODENUM $METADATA | cut -f4)_$(grep -w ^$CODENUM $METADATA | cut -f5)_$(grep -w ^$CODENUM $METADATA | cut -f6)_$CODENUM
NAME=$(echo $NAME | sed 's#/#_#g')

bedtools intersect -wb \
 -a <(sed 's/^chr//' $BEDIN ) \
 -b <(awk -v d=$DIST '{OFS="\t"} {print $1,$2-d,$3+d,$4}' $TSS \
| awk '{OFS="\t"} { if ($2<1) print $1,0,$3+1000,$4 ; else print $0}')  \
| sort -k8gr | awk '!arr[$14]++ {print $14}' \
| cut -d '_' -f2- | awk '!arr[$1]++ {print $1}' | head -$MAX_GENES \
| tr '\n' '\t' | sed "s#^#${NAME}\tCistromeDB\t#" | sed 's/$/\n/'
}
export -f map_pk
parallel -j 16 map_pk ::: TF_mouse/*narrowPeak.bed ::: mouse_factor_full_QC.txt ::: $ENH ::: 1000 ::: 1000 \
| sed 's/ //g' > MouseTfPeaks_enh.gmt
