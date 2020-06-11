#!/usr/bin/python
#Shift5endOffset_BED13RPF.py
#This script takes RPF bed14 file as input (col14 is Y/N 5end mismatch tag), modify the reads based on the given p site offset.
#Normally the offset is between [10,15], and won't be longer than the minimal read length (usually 18 nt).
#Column 2,3,7,8,10,11,12 may be modified, and the output read coordinate length is Original length - P offset (You can know the real length from the col13).
#One possible drawback is that we cannot extend the right end since we don't know the annotation (whether is may cross splicing site or not).
#But for RPF analysis, we are more interested in the 5end, so this method is still very useful.
#Suggested output name: BED13.XntShifted.RPF
#Version: Yu Sun, 2018-09-27 
#Version: Yu Sun, 2018-12-29

import sys

def ReadPOffset(Default):
    DefaultOffset = Default
    fi = open(sys.argv[3], "r")
    data = fi.readlines()[0]
    print("Reading in the following P offset dict: \n" + data.strip())
    offdict = {}
    MMoffdict = {}
    exec (data)
    MMoffdict = offdict['m0']
    for i in range(26, 33):
        if i in offdict:
            if i in MMoffdict:
                print "Length: " + str(i) + "nt, 5end match Offset (From file)  : " + str(offdict[i]) + ", 5end mismatch Offset (From file)  : " + str(MMoffdict[i])
            else:
                MMoffdict[i] = DefaultOffset + 1
                print "Length: " + str(i) + "nt, 5end match Offset (From file)  : " + str(offdict[i]) + ", 5end mismatch Offset (Use default): " + str(MMoffdict[i])
        else:
            offdict[i] = DefaultOffset
            if i in MMoffdict:
                print "Length: " + str(i) + "nt, 5end match Offset (Use default): " + str(offdict[i]) + ", 5end mismatch Offset (From file)  : " + str(MMoffdict[i])
            else:
                MMoffdict[i] = DefaultOffset + 1
                print "Length: " + str(i) + "nt, 5end match Offset (Use default): " + str(offdict[i]) + ", 5end mismatch Offset (Use default): " + str(MMoffdict[i])
    return offdict, MMoffdict

def Converter(offdict, MMoffdict):
    fi=open(sys.argv[1],'r')
    Offset=int(sys.argv[2])      #All numbers are strings when read in
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
        CurrSeqLength=len(CoorList[12])
        Tag=CoorList[13]

        if Tag == "Y":
            #This record has 5end mismatch, use MM dict:
            Currdict=MMoffdict.copy()
        else:
            Currdict=offdict.copy()

        #Get the offset for the current record:
        Offset=Currdict[CurrSeqLength]

        if ExonLenBlock[-1:] == ",":
            ExonLenBlock_clean=ExonLenBlock[:-1]
        else:
            ExonLenBlock_clean=ExonLenBlock

        if ExonStartBlock[-1:] == ",":
            ExonStartBlock_clean=ExonStartBlock[:-1]
        else:
            ExonStartBlock_clean=ExonStartBlock

        if ExonNum == 1:
            if Strand == "+":
                NewExonBlock=int(ExonLenBlock_clean)-Offset
                fo.write(CoorList[0]+"\t"+str(Start+Offset)+"\t"+str(End)+"\t"+CoorList[3]+"\t"+CoorList[4]+"\t"+Strand+"\t"+str(Start+Offset)+"\t"+str(End)+"\t"+CoorList[8]+"\t"+str(ExonNum)+"\t"+str(NewExonBlock)+"\t"+ExonStartBlock+"\t"+CoorList[12]+"\n")
            elif Strand == "-":
                NewExonBlock=int(ExonLenBlock_clean)-Offset
                fo.write(CoorList[0]+"\t"+str(Start)+"\t"+str(End-Offset)+"\t"+CoorList[3]+"\t"+CoorList[4]+"\t"+Strand+"\t"+str(Start)+"\t"+str(End-Offset)+"\t"+CoorList[8]+"\t"+str(ExonNum)+"\t"+str(NewExonBlock)+"\t"+ExonStartBlock+"\t"+CoorList[12]+"\n")
        else:
            if Strand == "+":
                ExonLenBlock_clean_list=ExonLenBlock_clean.split(",")
                ExonStartBlock_clean_list=ExonStartBlock_clean.split(",")
                SumBlocks=0
                for i in range(0,ExonNum):
                    SumBlocks=SumBlocks+int(ExonLenBlock_clean_list[i])
                    if SumBlocks > Offset:
                        ExonID=i  #Real Exon that the 5end falls in: i+1
                        ExonOverlap=int(ExonLenBlock_clean_list[i])-(SumBlocks-Offset)
                        break
                NewStart=Start+int(ExonStartBlock_clean_list[ExonID])+ExonOverlap
                NewExonNum=ExonNum-ExonID
                NewLenBlock=[0]*NewExonNum
                NewStartBlock=[0]*NewExonNum
                for j in range(0,NewExonNum):
                    if j==0:
                        NewLenBlock[j]=SumBlocks-Offset
                        NewStartBlock[j]=0
                    else:
                        NewLenBlock[j]=int(ExonLenBlock_clean_list[ExonID+j])
                        NewStartBlock[j]=int(ExonStartBlock_clean_list[ExonID+j])-int(ExonStartBlock_clean_list[ExonID])-ExonOverlap
                fo.write(CoorList[0]+"\t"+str(NewStart)+"\t"+str(End)+"\t"+CoorList[3]+"\t"+CoorList[4]+"\t"+Strand+"\t"+str(NewStart)+"\t"+str(End)+"\t"+CoorList[8]+"\t"+str(NewExonNum)+"\t"+",".join(map(str,NewLenBlock))+"\t"+",".join(map(str,NewStartBlock))+"\t"+CoorList[12]+"\n")
            else:
                ExonLenBlock_clean_list=ExonLenBlock_clean.split(",")
                ExonStartBlock_clean_list=ExonStartBlock_clean.split(",")
                SumBlocks=0
                for i in range(0,ExonNum):
                    SumBlocks=SumBlocks+int(ExonLenBlock_clean_list[ExonNum-i-1])
                    if SumBlocks > Offset:
                        ExonID=i  #Real Exon that the 5end falls in (count from right to left): i+1
                        ExonOverlap=SumBlocks-Offset
                        break
                NewEnd=Start+int(ExonStartBlock_clean_list[ExonNum-ExonID-1])+ExonOverlap
                NewExonNum=ExonNum-i
                NewLenBlock=[0]*NewExonNum
                NewStartBlock=[0]*NewExonNum
                for j in range(0,NewExonNum):
                    if j==NewExonNum-1:
                        NewLenBlock[j]=SumBlocks-Offset
                        NewStartBlock[j]=int(ExonStartBlock_clean_list[ExonNum-ExonID-1])
                    else:
                        NewLenBlock[j]=int(ExonLenBlock_clean_list[j])
                        NewStartBlock[j]=int(ExonStartBlock_clean_list[j])
                fo.write(CoorList[0]+"\t"+str(Start)+"\t"+str(NewEnd)+"\t"+CoorList[3]+"\t"+CoorList[4]+"\t"+Strand+"\t"+str(Start)+"\t"+str(NewEnd)+"\t"+CoorList[8]+"\t"+str(NewExonNum)+"\t"+",".join(map(str,NewLenBlock))+"\t"+",".join(map(str,NewStartBlock))+"\t"+CoorList[12]+"\n")

if len(sys.argv) != 5: #if the length of argv is not equal to 5, then print warning message
    print "This script takes RPF bed14 file as input (col14 is Y/N 5end mismatch tag), modify the reads based on the RoTISH P offset, and using a given p site offset as default (12nt usually)."
    print "Normally the offset is between [10,15], and won't be longer than the minimal read length (usually 26 nt for size selected bed13.RPF files)."
    print "Column 2,3,7,8,10,11,12 may be modified, and the output read coordinate length is Original length - P offset (You can know the real length from the col13)."
    print "One possible drawback is that we cannot extend the right end since we don't know the annotation (whether is may cross splicing site or not)."
    print "But for RPF analysis, we are more interested in the 5ends, so this method is still beneficial."
    print "Providing an empty offset file with only this line: offdict = {'m0':{}}, then this program will use Default P_Offset"
    print "Usage: [Shift5endOffset_BED14RPF.py] [BED14.RPF] [Default P_Offset] [Offset dict|RiboTISH.para.py] [Output|BED13.PShifted.RPF]"
    print "     Suggested output name: BED13.PShifted.RPF"
else:
    DefaultOffset=int(sys.argv[2])
    offdict, mmoffdict = ReadPOffset(DefaultOffset)
    print "Offsets used for 5end correction:"
    print offdict
    print "Offsets for 5end mismatches reads:"
    print mmoffdict
    Converter(offdict, mmoffdict)
