#!/usr/bin/python
#BedGraphRegionFiller.py
#This script takes benGraph as input, conver it into a single-base bedGraph for figure plotting
#There is another script written using bash which does exact the same thing: BedGraphRegionFiller.sh, but this python script is much much faster.
#Version: Yu Sun, 2018-10-27
#Version: Yu Sun, 2022-09-14, update for Python3
                                                                          
import sys

def Counter():
    fi=open(sys.argv[1],'r')
    fo=open(sys.argv[2],'w')
    for Line in fi:
        Data=Line.strip()
        line=Data.split()
        for i in range(int(line[2])-int(line[1])):
            fo.write(line[0]+"\t"+str(int(line[1])+i)+"\t"+str(int(line[1])+i+1)+"\t"+line[3]+"\n")
    fi.close()
    fo.close()
if len(sys.argv) != 3:
    print("This script takes benGraph as input, conver it into a single-base bedGraph for figure plotting")
    print("There is another script written using bash which does exact the same thing: BedGraphRegionFiller.sh, but this python script is much much faster.")
    print("Usage: [Template.py] [Input.bedGraph] [Output.bedGraph.sb]")
else:
    Counter()
