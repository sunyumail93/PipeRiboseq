#!/usr/bin/python
#TagBED14Mismatch5end.py
#This script adds a tag to indicate 5end mismatch from bed14 input (col13 is original sequenced reads, and col14 is genome sequences)
#The output has a modified col14 with Y (for 5-end mismatch) and N for non-mismatch
#Version: Yu Sun, 2018-12-29
#Version: Yu Sun, 2022-09-14, update for Python3
                                                                          
import sys

def Counter():
    fi=open(sys.argv[1],'r')
    fo=open(sys.argv[2],'w')
    
    for line in fi:
        Data=line.strip().split()
        NT_Seq=Data[12][0]
        NT_Gen=Data[13][0]
        if NT_Seq == NT_Gen:
            Tag="N"
        else:
            Tag="Y"
        fo.write(Data[0]+"\t"+Data[1]+"\t"+Data[2]+"\t"+Data[3]+"\t"+Data[4]+"\t"+Data[5]+"\t"+Data[6]+"\t"+Data[7]+"\t"+Data[8]+"\t"+Data[9]+"\t"+Data[10]+"\t"+Data[11]+"\t"+Data[12]+"\t"+Tag+"\n")
    fi.close()
    fo.close()

if len(sys.argv) != 3:
    print("This script adds a tag to indicate 5-end mismatch from bed14 input (col13 is original sequenced reads, and col14 is genome sequences)")
    print("The output has a modified col14 with Y (for 5-end mismatch) and N for non-mismatch")
    print("Usage: [TagBED14Mismatch5end.py] [Input.bed14] [Output.bed14]")
else:
    Counter()
