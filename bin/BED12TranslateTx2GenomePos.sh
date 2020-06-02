#!/usr/bin/env bash
#BED12TranslateTx2GenomePos.sh
#This is a translator of transcript locations to gemomic locations, based on the BED12 annotation
#The column 4 (names) in the annotation file must be unique, and match the names in the transcript file
#Awk commands, running on linux/mac terminal
#Both of the two input files are 0-based
#Version: Yu Sun, 2018-06-05 ~ 2018-06-08

array=( "$@" ) #read all command line variable to "array"
if [ ! -n "$6" ]
then
  echo "    This is a translator of transcript locations to gemomic locations, based on the BED12 annotation"
  echo "    The column 4 (names) in the annotation file must be unique, and match the names in the transcript file"
  echo "    Both of the two input files are 0-based, and please try to avoid white characters in names"
  echo "    Usage: `basename $0` -i [Genes.bed12] -t [Transcript.bed|Only require the first three columns :TxName start end] -o [Output.bed12]"
  echo "    For each record in the Transcript.bed, it will generate a corresponding bed12 record in the Output.bed12"
  echo "    Example of back-calculating test: gene.cds.bed12 should be same as gene.cds.bed12.validate"
  echo "       BED12Extractor.sh -a cds -i gene.bed12 -o gene.cds.bed12"
  echo "       BED12TranslateGenomePos2Tx.sh -i gene.bed12 -t gene.cds.bed12 -o gene.cds.tx"
  echo "       BED12TranslateTx2GenomePos.sh -i gene.bed12 -t gene.cds.tx -o gene.cds.bed12.validate"
else
  echo "Starting BED12TranslateTx2GenomePos.sh pipeline"
  echo "1. Getting parameters"

for arg in "$@"
do
 if [[ $arg == "-i" ]]
  then
    Anno=${array[$counter+1]}
    echo '   Getting annotation: '$Anno
 elif [[ $arg == "-t" ]]    
  then
    Data=${array[$counter+1]}
    echo '   Getting Tx bed: '$Data
 elif [[ $arg == "-o" ]]
  then
    Output=${array[$counter+1]}
    echo '   Getting output name: '$Output
 fi
  let counter=$counter+1
done

  echo "2. Translating positions"
rm -rf $Output
RecordNumber=`wc -l $Data|awk '{print $1}'`
for i in $(eval echo {1..$RecordNumber})
do
  #echo $i
  Tx=`sed -n ${i}p $Data`
  TxName=`echo $Tx | awk '{print $1}'`
  TxStart=`echo $Tx | awk '{print $2}'`
  TxEnd=`echo $Tx | awk '{print $3}'`
  CurrAnno=`awk -v name=$TxName '{if ($4==name) print $0}' $Anno`
  #echo $CurrAnno
  #echo $TxStart,$TxEnd

  echo $CurrAnno | awk -v TxStart=$TxStart -v TxEnd=$TxEnd '{OFS="\t";split($11,blockSizes,","); split($12,blockStarts,","); blockCount=$10;
  error=0;A=""; B="";start=0;end=0;N=0;TxLength=TxEnd-TxStart;
  if (TxStart >= TxEnd || TxStart<0 || TxEnd<0 ) {error=1}
  for (i=1;i<=$10;i++) {FullLength=FullLength+blockSizes[i]}
  if ($6=="+") {
    blocksum[0]=0;everstart=0;everend=0;
    for (i=1;i<=$10;i++){ blocksum[i]=blocksum[i-1]+blockSizes[i];
      if (TxStart<blocksum[i] && everstart==0){startexopos=i;everstart=1;}
      if (TxEnd<=blocksum[i] && everend==0){endexopos=i;everend=1;}
    }
    if (startexopos==endexopos){
      N=1;A=TxEnd-TxStart",";B=0",";start=blockStarts[startexopos]+TxStart-blocksum[startexopos-1]+$2;end=blockStarts[endexopos]+TxEnd-blocksum[endexopos-1]+$2;}
    if (startexopos<endexopos){N=endexopos-startexopos+1;
      for (j=1;j<=N;j++) {
        if (j==1){tempA=blocksum[startexopos]-TxStart;A=tempA;B=B"0";
        }else if (j==N) {tempA=TxEnd-blocksum[endexopos-1];A=A","tempA;
          B=B","blockStarts[endexopos]-blockStarts[startexopos]-(TxStart-blocksum[startexopos-1]);
        }else {A=A","blockSizes[j-1+startexopos];
          B=B","blockStarts[j-1+startexopos]-blockStarts[startexopos]-(TxStart-blocksum[startexopos-1]);
        }  
      }
      A=A",";B=B",";start=blockStarts[startexopos]+TxStart-blocksum[startexopos-1]+$2;end=blockStarts[endexopos]+TxEnd-blocksum[endexopos-1]+$2;
    }
    if (everend==0){error=1}
  }
  if ($6=="-" && blockCount>=1) {
    blocksum[0]=0;everstart=0;everend=0;
    for (i=1;i<=$10;i++){ blocksum[i]=blocksum[i-1]+blockSizes[$10-i+1];
      if (TxStart<blocksum[i] && everstart==0){startexopos=i;everstart=1;}
      if (TxEnd<=blocksum[i] && everend==0){endexopos=i;everend=1;}
    }
    if (startexopos==endexopos){
      N=1;A=TxEnd-TxStart",";B=0",";start=blockStarts[$10-startexopos+1]+blocksum[startexopos]-TxEnd+$2;end=blockStarts[$10-endexopos+1]+blocksum[endexopos]-TxStart+$2;
    }
    if (startexopos<endexopos){N=endexopos-startexopos+1;
      for (j=1;j<=N;j++) {
        if (j==1){tempA=blocksum[startexopos]-TxStart;A=tempA",";
          B=blockStarts[$10-startexopos+1]-blockStarts[$10-endexopos+1]-(blocksum[endexopos]-TxEnd)",";
        }else if (j==N) {tempA=TxEnd-blocksum[endexopos-1];A=tempA","A;B="0,"B;
        }else {A=blockSizes[$10-(j-1+startexopos)+1]","A;
          B=blockStarts[$10-(j-1+startexopos)+1]-blockStarts[$10-endexopos+1]-(blocksum[endexopos]-TxEnd)","B;
        }  
      }
      start=blockStarts[$10-endexopos+1]+(blocksum[endexopos]-TxEnd)+$2;end=blockStarts[$10-startexopos+1]+(blocksum[startexopos]-TxStart)+$2;
    }
    if (everend==0){error=1}
  }
  if (error==0) {
    print $1,start,end,$4,$5,$6,start,end,$9,N,A,B
  }else print "Out range error: "$4,TxStart,TxEnd | "cat 1>&2"  }' >> $Output

done

if [ -s $Output ]
  then
    echo "3. Done"
else
   echo "Error occurred!"
fi
  
fi
