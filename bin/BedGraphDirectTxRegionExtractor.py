#!/usr/bin/python
#BedGraphDirectTxRegionExtractor.py
#This is an updated version of BedGraphRegionExtractor.py
#This script takes bed12 and bedGraph as input, generate a spliced single base resolution of counts
#The output uses genomic coordinates, but the intron regions will be spliced out
#The fourth column can be directly plotted by R
#Version: 2018-09-10, Yu Sun
                                                                          
import sys

def Counter():
    fi=open(sys.argv[1],'r')
    fibed=open(sys.argv[2],'r')
    fo=open(sys.argv[3],'w')
    
    bed12=fi.readlines()  #read all lines into array
    Chr=bed12[0].split()[0]
    Start=int(bed12[0].split()[1])
    End=int(bed12[0].split()[2])
    Name=bed12[0].split()[3]
    Strand=bed12[0].split()[5]
    ExonNum=int(bed12[0].split()[9])
    BlocksLength=bed12[0].split()[10]
    StartBlocks=bed12[0].split()[11]

    BlocksLengthList=BlocksLength.split(",")
    StartBlocksList=StartBlocks.split(",")
    GenomeStarts=["z"]*ExonNum
    GenomeEnds=["z"]*ExonNum

    #Get bed6 coordinates
    if ExonNum == 1 :
        GenomeStarts[0]=int(Start)+int(StartBlocksList[0])
        GenomeEnds[0]=int(Start)+int(StartBlocksList[0])+int(BlocksLengthList[0])
    else:
        for i in range(0,ExonNum):
            GenomeStarts[i]=int(Start)+int(StartBlocksList[i])
            GenomeEnds[i]=int(Start)+int(StartBlocksList[i])+int(BlocksLengthList[i])

    #Manipulate bedGraph
    bedGraphs=fibed.readlines()
    
    #Compare single-base coordinates and bedGraph
    for i in range(0,ExonNum):
        for j in range(GenomeStarts[i],GenomeEnds[i]):
            for line in bedGraphs:
#                print line
                if (int(line.split()[1]) <= j) & (j+1 <= int(line.split()[2])):    #The parentheses is required!!!
#                    print int(line.split()[1]) <= j & j+1 <= int(line.split()[2])
#                    print j
                    if line.split()[0] == Chr:
#                        print Chr+"\t"+str(j)+"\t"+str(j+1)+"\t"+line.split()[3]
                        fo.write(Chr+"\t"+str(j)+"\t"+str(j+1)+"\t"+line.split()[3]+"\n")

#    print GenomeStarts
#    print GenomeEnds

    fi.close()
    fibed.close()
    fo.close()

if len(sys.argv) != 4: #if the length of argv is not equal to 3, then print warning message
    print "This is an updated version of BedGraphRegionExtractor.py, splicing intron directly"
    print "This script takes single-line bed12 and a bedGraph as input, generate a spliced single base resolution of counts"
    print "The output uses genomic coordinates, but the intron regions will be spliced out"
    print "The fourth column can be directly plotted by R"
    print "Usage: [BedGraphDirectTxRegionExtractor.py] [Region.bed12] [Input.bedGraph] [Output]"
else:
    Counter()
