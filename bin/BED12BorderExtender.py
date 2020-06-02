#!/usr/bin/python
#BED12BorderExtender.py
#This script takes Annotation bed12 as input, extend col2 and col3 borders, considering the splicing
#This is similar to bedtools slop, but considering the splicing. So only the first and last exon length may change.
#Version: Yu H. Sun, 2020-06-02
                                                                          
import sys

def Counter():
    fi=open(sys.argv[1],'r')
    Left=int(sys.argv[2])
    Right=int(sys.argv[3])
    fo=open(sys.argv[4],'w')
    
    for line in fi:
        Coor=line.strip()
        CoorList=Coor.split()
        Start=int(CoorList[1])
        End=int(CoorList[2])
        Strand=CoorList[5]
        ORFStart=int(CoorList[6])
        ORFEnd=int(CoorList[7])
        ExonNum=int(CoorList[9])
        ExonLenBlock=CoorList[10]
        ExonStartBlock=CoorList[11]

        if ExonLenBlock[-1:] == ",":
            ExonLenBlock_clean=ExonLenBlock[:-1]
        else:
            ExonLenBlock_clean=ExonLenBlock

        if ExonStartBlock[-1:] == ",":
            ExonStartBlock_clean=ExonStartBlock[:-1]
        else:
            ExonStartBlock_clean=ExonStartBlock    

        if ExonNum == 1:
            NewExonBlock=int(ExonLenBlock_clean)+Left+Right
            if Strand == "+":
                fo.write(CoorList[0]+"\t"+str(Start-Left)+"\t"+str(End+Right)+"\t"+CoorList[3]+"\t"+CoorList[4]+"\t"+Strand+"\t"+CoorList[6]+"\t"+CoorList[7]+"\t"+CoorList[8]+"\t"+str(ExonNum)+"\t"+str(NewExonBlock)+"\t"+ExonStartBlock+"\n")
            elif Strand == "-":
                fo.write(CoorList[0]+"\t"+str(Start-Right)+"\t"+str(End+Left)+"\t"+CoorList[3]+"\t"+CoorList[4]+"\t"+Strand+"\t"+CoorList[6]+"\t"+CoorList[7]+"\t"+CoorList[8]+"\t"+str(ExonNum)+"\t"+str(NewExonBlock)+"\t"+ExonStartBlock+"\n")
        else:
            ExonLenBlock_clean_list=ExonLenBlock_clean.split(",")
            ExonStartBlock_clean_list=ExonStartBlock_clean.split(",")
            if Strand == "+":
                #Mutate start and end exons
                ExonLenBlock_clean_list[0]=int(ExonLenBlock_clean_list[0])+Left
                ExonLenBlock_clean_list[-1]=int(ExonLenBlock_clean_list[-1])+Right
                for i in range(len(ExonStartBlock_clean_list)):
                    if i!=0:
                        ExonStartBlock_clean_list[i]=int(ExonStartBlock_clean_list[i])+Left
                fo.write(CoorList[0]+"\t"+str(Start-Left)+"\t"+str(End+Right)+"\t"+CoorList[3]+"\t"+CoorList[4]+"\t"+Strand+"\t"+CoorList[6]+"\t"+CoorList[7]+"\t"+CoorList[8]+"\t"+str(ExonNum)+"\t"+",".join(map(str,ExonLenBlock_clean_list))+"\t"+",".join(map(str,ExonStartBlock_clean_list))+"\n")                
            else:
                ExonLenBlock_clean_list[0]=int(ExonLenBlock_clean_list[0])+Right
                ExonLenBlock_clean_list[-1]=int(ExonLenBlock_clean_list[-1])+Left
                for i in range(len(ExonStartBlock_clean_list)):
                    if i!=0:
                        ExonStartBlock_clean_list[i]=int(ExonStartBlock_clean_list[i])+Right
                fo.write(CoorList[0]+"\t"+str(Start-Right)+"\t"+str(End+Left)+"\t"+CoorList[3]+"\t"+CoorList[4]+"\t"+Strand+"\t"+CoorList[6]+"\t"+CoorList[7]+"\t"+CoorList[8]+"\t"+str(ExonNum)+"\t"+",".join(map(str,ExonLenBlock_clean_list))+"\t"+",".join(map(str,ExonStartBlock_clean_list))+"\n")

if len(sys.argv) != 5: #if the length of argv is not equal to 5, then print warning message
    print "This script takes Annotation bed12 as input, extend col2 and col3 borders, considering the splicing"
    print "This is similar to bedtools slop, but considering the splicing"
    print "Usage: [BED12BorderExtender.py] [Gene.bed12] [LeftExtension] [RightExtension] [Output.ext.bed12]"
else:
    Counter()
