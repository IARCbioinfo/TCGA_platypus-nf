#!/bin/bash

SM=`grep "^#CHROM" $1 | awk '{print $10}'`	#get the sample name from the header of input vcf

uniqID=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 16 | head -n 1)
res=$SM"_"$uniqID"_reformat.tsv"						#get the name of the output

grep -v "^#" $1 | awk '{ 							#duplicate the column start to have start and end (both equal at this time)
        OFS="	"; FS="	"
    } {
        $2=$2 OFS $2;
    } {
        print $0
    }
' | awk '{ 									#compute the true end for deletions
        OFS="	"; FS="	"
    }
    length($5) > length($6) {
       $3 = $2 + length($5) - 1
    } { 									# move the column 4 at the end (ID), to have chr-start-end-ref-alt-[...]
       a=$4; for (i=4;i<NF; i++) $i=$(i+1); $NF=a; print
    }
' > $res

sed -i "s/$/\t$SM/" $res							#add the sample name as a last column

TSS=`echo $SM | tr "-" "     " | awk '{print $2}'`

sourceSite=`awk -F "	" -v tss="$TSS" '$1 == tss {print $3}' $2`		#get the tissue source site from sample barcode

studyAbb=`awk -F "	" -v ss="$sourceSite" '$2 == ss {print $1}' $3`		#get the cancer tissue abbreviation

sed -i "s/$/\t$studyAbb/" $res							#add the cancer tissue abbreviation as a last column
