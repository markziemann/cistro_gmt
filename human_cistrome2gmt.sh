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
GTFURL=ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
GTF=`basename $GTFURL`
TSSBED=`echo $GTF | sed 's/.gtf.gz/.tss.bed/'`
SIZE_RANGE='2500,1000,500'
ENHANCERURL=https://genecards.weizmann.ac.il/geneloc/gh.gff
ENHANCERZIP=`basename $ENHANCERURL`
GNAMES=$(basename $GTF .gtf.gz).gnames.txt
ENH_BED19=enhancers.hg19.bed
ENH_BED38=enhancers.hg38.bed
#this encode stuff will have to be modified
METADATA_URL="https://www.encodeproject.org/metadata/type=Experiment&assay_title=ChIP-seq&target.investigated_as=transcription+factor&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&files.file_type=bed+narrowPeak/metadata.tsv"
METADATA="TF_human_data_information.txt"
METADATA_SUMMARY=metadata_summary.tsv
URLLIST=files.txt

NUMRANGE='500,1000'
UPPER_LIMIT='5000'
TSSGMT=human_tss_TFBS.gmt
ENHGMT=human_enhancer_TFBS.gmt
NRCPU=`nproc`

###############################################################
echo "Dependancies OK. Downloading required annotation files"
###############################################################
wget -N $GTFURL
wget -N $ENHANCERURL

###############################################################
echo "extracting gene names from GTF"
###############################################################
zcat $GTF | grep -w gene | cut -f9 | cut -d '"' -f2,6 \
| tr '"' '\t' | sort -uk 2b,2 > $GNAMES

###############################################################
echo "extracting enhancer coordinates"
###############################################################

myfunc(){
LINE="$*"
COORD=$(echo $LINE | cut -d ' ' -f1)
echo $LINE | cut -d ' ' -f2- | tr ' ' '\n' | sed "s/^/${COORD}\t/;s/:/\t/;s/-/\t/"
}
export -f myfunc

sed 1d gh.gff  | tr '=' '\t' \
| awk -F"\t" '{print $1,$4,$5,$NF}' \
| sed 's/ /:/;s/ /-/;s/,//g' \
| parallel -X myfunc \
| sed 's/chr//' | bedtools sort > $ENH_BED38

grep -v ENSG $ENH_BED38  | sort -k 4b,4 \
| join -1 4 -2 2 - $GNAMES | tr ' ' '\t' | cut -f2- > $ENH_BED38.tmp

grep ENSG $ENH_BED38 >> $ENH_BED38.tmp

bedtools sort -i $ENH_BED38.tmp > $ENH_BED38 && rm $ENH_BED38.tmp

###############################################################
echo "extracting TSS positions from ensembl GTF"
###############################################################
#Get coordinates of plus strand genes
zgrep -w 'exon_number "1"' $GTF | awk '$7=="+"' \
| cut -f1,4,5,7,9 | cut -d '"' -f-2,12 \
| sed 's!gene_id "!!' | tr '"' '_'  \
| awk '{OFS="\t"} {print $1,$2,$2+1,$5,$4}' \
| awk '{OFS="\t"} { if ($2<1) print $1,"1",$3,$5,$4 ; else print $0 }' > $TSSBED

#Get coordinates of minus strand genes
zgrep -w 'exon_number "1"' $GTF | awk '$7=="-"' \
| cut -f1,4,5,7,9 | cut -d '"' -f-2,12 \
| sed 's!gene_id "!!' | tr '"' '_'  \
| awk '{OFS="\t"} {print $1,$3-1,$3,$5,$4}' \
| awk '{OFS="\t"} { if ($2<1) print $1,"1",$3,$5,$4 ; else print $0 }' >> $TSSBED

###############################################################
echo "intersect peaks with TSS"
###############################################################

map_pk(){
BEDIN=$1
METADATA=$2
TSS=$3
DIST=$4
MAX_GENES=$5

OUT=$(echo $BEDIN | sed 's/narrowPeak.bed/narrowPeak.gmt/')
NAME=$(grep -w $(basename $BEDIN) $METADATA | cut -f7)_$(grep -w $(basename $BEDIN) $METADATA | cut -f4)_$(grep -w $(basename $BEDIN) $METADATA | cut -f2)
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
parallel map_pk ::: TF_human/*narrowPeak.bed ::: TF_human_data_information.txt ::: Homo_sapiens.GRCh38.86.tss.bed ::: 1000 ::: 1000 \
| sed 's/ //g' > HumanTfPeaks.gmt
