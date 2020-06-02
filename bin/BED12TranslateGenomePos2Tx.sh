#!/usr/bin/env bash
#BED12TranslateGenomePos2Tx.sh
#This is a translator of gemomic locations to transcript positions, based on the BED12 annotation
#The column 4 (names) in the annotation file must be unique, and match the names in the genomic location file
#Awk commands, running on linux/mac terminal
#Both of the two input files are 0-based
#Version: Yu Sun, 2018-06-10
#Version: Yu Sun, 2018-06-12, fix intron location bug (if the genomic pos is in intron), output error

array=( "$@" ) #read all command line variable to "array"
if [ ! -n "$6" ]
then
  echo "    This is a translator of gemomic locations to transcript positions, based on the BED12 annotation"
  echo "    The column 4 (names) in the annotation file must be unique, and match the names in the genomic location file"
  echo "    Input and output files are 0-based, and please try to avoid white characters in names"
  echo "    Records in GenomePos.bed4 can duplicate: multiple records with same column 4(names) value, and you can directly use the bed12 file as GenomePos.bed4 (such as extracted cds bed12 annotation)"
  echo "    Usage: `basename $0` -i [Genes.bed12] -t [GenomePos.bed4|Only require the first four columns : chr start end TxName] -o [Tx.bed3]"
  echo "    For each record in the GenomePos.bed4, it will generate a corresponding bed12 record in the Tx.bed3"
  echo "    Example of back-calculating test: gene.cds.bed12 should be same as gene.cds.bed12.validate"
  echo "       BED12Extractor.sh -a cds -i gene.bed12 -o gene.cds.bed12"
  echo "       BED12TranslateGenomePos2Tx.sh -i gene.bed12 -t gene.cds.bed12 -o gene.cds.tx"
  echo "       BED12TranslateTx2GenomePos.sh -i gene.bed12 -t gene.cds.tx -o gene.cds.bed12.validate"
else
  echo "Starting BED12TranslateGenomePos2Tx.sh pipeline"
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
    echo '   Getting GenomePos bed: '$Data
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
  chr=`echo $Tx | awk '{print $1}'`
  Name=`echo $Tx | awk '{print $4}'`
  GenoStart=`echo $Tx | awk '{print $2}'`
  GenoEnd=`echo $Tx | awk '{print $3}'`
  CurrAnno=`awk -v name=$Name -v chr=$chr '{if ($4==name && $1==chr) print $0}' $Anno`
  #echo $Tx
  #echo $CurrAnno
  #echo $GenoStart,$GenoEnd

  echo $CurrAnno | awk -v GenoStart=$GenoStart -v GenoEnd=$GenoEnd '{OFS="\t";split($11,blockSizes,","); split($12,blockStarts,","); blockCount=$10;
  error=0;start=0;end=0;everstart=0;everend=0;
  if (GenoStart >= GenoEnd || GenoStart<$2 || GenoEnd>$3 ) {error=1}
  for (i=1;i<=$10;i++) {FullLength=FullLength+blockSizes[i]}
  if ($6=="+") {
    blocksum[0]=0;
    for (i=1;i<=$10;i++) {
      blocksum[i]=blocksum[i-1]+blockSizes[i];
      if (GenoStart>=$2+blockStarts[i] && GenoStart <=$2+blockStarts[i]+blockSizes[i]) {
        start=blocksum[i-1]+GenoStart-blockStarts[i]-$2;
        everstart=1;
      }
      if (GenoEnd>=$2+blockStarts[i] && GenoEnd <=$2+blockStarts[i]+blockSizes[i]) {
        end=blocksum[i-1]+GenoEnd-blockStarts[i]-$2;
        everend=1;
      }
    }
  }
  if ($6=="-") {
    blocksum[0]=0;
    for (i=1;i<=$10;i++) {
      blocksum[i]=blocksum[i-1]+blockSizes[$10-i+1];
      if (GenoStart>=$2+blockStarts[$10-i+1] && GenoStart <=$2+blockStarts[$10-i+1]+blockSizes[$10-i+1]) {
        end=blocksum[i-1]+blockStarts[$10-i+1]+blockSizes[$10-i+1]-(GenoStart-$2);
        everend=1;
      }
      if (GenoEnd>=$2+blockStarts[$10-i+1] && GenoEnd <=$2+blockStarts[$10-i+1]+blockSizes[$10-i+1]) {
        start=blocksum[i-1]+blockStarts[$10-i+1]+blockSizes[$10-i+1]-(GenoEnd-$2);
        everstart=1;
      }
    }
  }
  if (everstart+everend<2) {error=1}
  if (error==0) {
    print $4,start,end
  }else print "Out range error: "$4,GenoStart,GenoEnd | "cat 1>&2"  }' >> $Output

done

if [ -s $Output ]
  then
    echo "3. Done"
else
   echo "Error occurred!"
fi
  
fi
