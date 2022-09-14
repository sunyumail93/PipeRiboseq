#!/usr/bin/python
#BedGraphDirectTxRegionSenseCoverage.py
#This script takes BedGraphDirectTxRegionExtractor.py calculated coverage and single gene bed12 annotation as input, decides strand and gets sense strand coverage
#If the gene is sense strand, nothing changed, if antise strand, then all lines will be reversed, and the minus sign will be removed
#Version: Yu Sun, 2018-09-11
#Version: Yu Sun, 2022-09-14, update for Python3
                                                                          
import sys

def Counter():
    fgene=open(sys.argv[1],'r')
    fsense=open(sys.argv[2],'r')
    fanti=open(sys.argv[3],'r')
    fo=open(sys.argv[4],'w')

    bed12=fgene.readlines()
    Strand=bed12[0].split()[5]

    if Strand == "+":
        SenseCoverage=fsense.readlines()
        for line in SenseCoverage:
            fo.write(line)
    else:
        AntiseneCoverage=fanti.readlines()
        for j in range(0,len(AntiseneCoverage)):
            CurrLine=AntiseneCoverage[len(AntiseneCoverage)-1-j].strip()
            CurrLineSplit=CurrLine.split()
            NewValue=CurrLineSplit[3].replace("-","")
            fo.write(CurrLineSplit[0]+"\t"+CurrLineSplit[1]+"\t"+CurrLineSplit[2]+"\t"+NewValue+"\n")

if len(sys.argv) != 5: #if the length of argv is not equal to 5, then print warning message
    print("This script takes BedGraphDirectTxRegionExtractor.py calculated coverage and single gene bed12 annotation as input, decides strand and gets sense strand coverage")
    print("If the gene is sense strand, nothing changed, if antise strand, then all lines will be reversed, and the minus sign will be removed")
    print("Usage: [BedGraphDirectTxRegionSenseCoverage.py] [Singel gene: Gene.bed12] [Gene.plus.bedGraph.bb] [Gene.minus.bedGraph.bb] [Output: Gene.sense.bedGraph]")
else:
    Counter()
