#!/bin/bash

conv(){
line=$1
  NAME_DESC=`echo $line | cut -d ' ' -f-2`

  GENES=`echo $line | cut -d ' ' -f3- \
  | tr ' ' '\n' | sort -uk 1b,1 \
  | join -1 1 -2 1 - \
  <(cut -f3,5 mouse2hum_biomart_ens87.txt \
  | sed 1d | awk '$1!="" && $2!=""' \
  | sort -uk 1b,1) | cut -d ' ' -f2 \
  | sort -u | tr '\n' '\t' \
  | sed 's/\t$/\n/'`

  echo $NAME_DESC $GENES | tr ' ' '\t'
}
export -f conv
for GMT in `ls *gmt | grep -i mouse *gmt ` ; do
  NAME=`echo $GMT | sed 's/.gmt/_human.gmt/'`
  parallel -k conv < $GMT > $NAME
done

exit
