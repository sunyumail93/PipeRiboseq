#!/usr/bin/python
#BED13Finalizer.py
#This script takes original bed13 as input, finalize it by reversing col13 to be the real sequence when the read mapped to reverse strand
#Input: get from bed12 (convert bam to bed12) and paste the col10 of sam/bam file
#Output: a correct bed13
#Version: Yu Sun, 2018-12-29
#Version: Yu Sun, 2022-09-14, update for Python3
                                                                          
import sys

def Complementer(String):
    StringList=list(String)
    for i in range(len(StringList)):
    #Substitution based on 'Nucleic acid notation'
    #For canonical bases:
        if StringList[i]=="A": StringList[i]="T"
        elif StringList[i]=="a": StringList[i]="t"
        elif StringList[i]=="T": StringList[i]="A"
        elif StringList[i]=="t": StringList[i]="a"
        elif StringList[i]=="G": StringList[i]="C"
        elif StringList[i]=="g": StringList[i]="c"
        elif StringList[i]=="C": StringList[i]="G"
        elif StringList[i]=="c": StringList[i]="g"
    #For degenerated bases:
        elif StringList[i]=="W": StringList[i]="S"
        elif StringList[i]=="w": StringList[i]="s"
        elif StringList[i]=="S": StringList[i]="W"
        elif StringList[i]=="s": StringList[i]="w"
        elif StringList[i]=="M": StringList[i]="K"
        elif StringList[i]=="m": StringList[i]="k"
        elif StringList[i]=="K": StringList[i]="M"
        elif StringList[i]=="k": StringList[i]="m"
        elif StringList[i]=="R": StringList[i]="Y"
        elif StringList[i]=="r": StringList[i]="y"
        elif StringList[i]=="Y": StringList[i]="R"
        elif StringList[i]=="y": StringList[i]="r"
        elif StringList[i]=="B": StringList[i]="V"
        elif StringList[i]=="b": StringList[i]="v"
        elif StringList[i]=="V": StringList[i]="B"
        elif StringList[i]=="v": StringList[i]="b"
        elif StringList[i]=="D": StringList[i]="H"
        elif StringList[i]=="d": StringList[i]="h"
        elif StringList[i]=="H": StringList[i]="D"
        elif StringList[i]=="h": StringList[i]="d"
        elif StringList[i]=="N": StringList[i]="N"
        elif StringList[i]=="n": StringList[i]="n"
        else: StringList[i]="N"
    return "".join(map(str,StringList))

def Counter():
    fi=open(sys.argv[1],'r')
    fo=open(sys.argv[2],'w')
    
    for line in fi:
        SplitLine=line.strip().split()
        Strand=SplitLine[5]
        Sequence=SplitLine[12]
        Reversed=Sequence[::-1]
        if Strand == "-":
            Comp=Complementer(Reversed)
        else:
            Comp=Sequence
#        print(Strand+"\t"+Comp)
        fo.write(SplitLine[0]+"\t"+SplitLine[1]+"\t"+SplitLine[2]+"\t"+SplitLine[3]+"\t"+SplitLine[4]+"\t"+SplitLine[5]+"\t"+SplitLine[6]+"\t"+SplitLine[7]+"\t"+SplitLine[8]+"\t"+SplitLine[9]+"\t"+SplitLine[10]+"\t"+SplitLine[11]+"\t"+Comp+"\n")

    fi.close()
    fo.close()

if len(sys.argv) != 3:
    print("This script takes original bed13 as input, finalize it by reversing col13 to be the real sequence when the read mapped to reverse strand")
    print("Input: get from bed12 (convert bam to bed12) and paste the col10 of sam/bam file")
    print("Usage: [BED13Finalizer.py] [InputBED13-Like File] [Output BED13]")
else:
    Counter()
