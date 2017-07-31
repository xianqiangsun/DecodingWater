############################################################################
##This file is used to calculate the exchange rate of the water molecules
##with respect to each center to check whether it is equilibrited
##This script reads the water coordinates according to each center
##and the exchange rate was analysised every 1 ns.
##
##
##Xianqiang Sun
##TheoChem&Bio
##KTH
##2012-10-25
##
###########################################################################

import numpy
import sys
sys.path.append('/home/x/xiansu/pfs/program/numpy/lib/python2.6/site-packages')
sys.path.append('/bubo/home/h8/xians/glob/program/python_pkg/lib64/python2.6/site-packages')
#from Numeric import *
from numpy import pi,cos,sin,arccos
from datetime import datetime

centreFile=open('optimizedCentre.txt','r')
waterFile=open('waterInforByCentre.txt','r')
H2OInforFile=open('WW_allCentre_H2O.txt','r')
##print centreFile

outputFile=open('exchangeRate.dat','w')

k=1.380648813*(10**(-23))   ## The unit of the boltzmann constant is J/K.
weiH2O=18.0154
mol=6.02214179*(10**23)
pw=0.0331725               ## The unit of this is No. of molecules in per A2.  
frameNo=6000
framePerNs=500          ##this value shows how many frames do you save in the simulation

def getCentreCoord(self):
    centreCoord=[]
    for centre in self:
        centre=centre.split()
        centreFloat=[]
        for coord in centre:
            coord =float(coord)
            centreFloat.append(coord)

        centreCoord.append(centreFloat)
##    print centreCoord
    return centreCoord

##getCentreCoord read the output file from the 'calculatCentre_oop.py'. Read the
##coorditates and save them as list with the float format. The output format is
##the format is:[[x1,y1,z2],[x2,y2,z2]...]

def getWaterInforWithCentre(self):
    waterInforWithCentre=[]
    for waterInfor in self:
        waterInfor=waterInfor.split()
        waterInforWithCentre.append(waterInfor)
    return waterInforWithCentre
##getWaterInforWithCentre is a function that can be used to read the output waterInfor
##file from 'OptimizeCentre_oop.py'. The water infor were saved as list,[[waterInfor]..]

def orangeWaterInforWithCentre(self,centre):
    waterInforAccCentre=[]
    for centreNo in range(len(centre)):
        eachWaterInforAccCentre=[]
        for waterInfor in self:
##            print waterInfor[-1]
            if waterInfor[-1]==str(centreNo):
                eachWaterInforAccCentre.append(waterInfor)
        waterInforAccCentre.append(eachWaterInforAccCentre)
    return waterInforAccCentre

##getWaterInforWithCentre is organise the water information to lists according the index of
##centres. the water information were saved as[[waterInforInCentre1],[waterInforInCentre1]..].
##WaterInfor in the final result were all saved as list.


def calculateExchangeRate(self,framePerNs):
    waterExchangeRate=[]
    for eachWaterNo in range(len(self)):
        eachWaterExchange=[0]*12
        eachWaterInforAccCenter=self[eachWaterNo]
##        oriFrameNo=eachWaterInforAccCenter[0][-2]
        oriSegname=eachWaterInforAccCenter[0][-3]
        oriResnumber=eachWaterInforAccCenter[0][2]
        oriEachWaterExchange=0
        oriTimeNs=0
        
        for eachFrameWater in eachWaterInforAccCenter:
            frameNo=int(eachFrameWater[-2])-1
            segname=eachFrameWater[-3]
            resnumber=eachFrameWater[2]
            timeNs=int(frameNo)/int(framePerNs)
            if segname==oriSegname and resnumber==oriResnumber and timeNs==oriTimeNs:
                continue
            elif segname!=oriSegname and timeNs==oriTimeNs:
                oriEachWaterExchange+=1
##                frameNo=int(eachFrameWater[-2])-1
                oriSegname=eachFrameWater[-3]
                #resnumber=eachFrameWater[2]
                continue
            elif resnumber!=oriResnumber and timeNs==oriTimeNs:
                oriEachWaterExchange+=1
##                frameNo=int(eachFrameWater[-2])-1
              #  segname=eachFrameWater[-3]
                oriResnumber=eachFrameWater[2]
                continue
            elif timeNs!=oriTimeNs:
                eachWaterExchange[oriTimeNs]=oriEachWaterExchange
                oriTimeNs=timeNs
                oriEachWaterExchange=0
                continue
        waterExchangeRate.append(eachWaterExchange)
    return waterExchangeRate
##This function reads the input of water information according to the center.
##Then, we will study the exchange rate of each water molecules in every nanosecond
##As I need to calculate the time according to how many frames did we saved in every nanosecond
##So, framePerNs is needed for the function.

def calculateAverageExchangeRate(self):
    allExchange=[]
    for eachWaterChange in self:
        eachAverage=float(sum(eachWaterChange))/len(eachWaterChange)
        eachWaterChange.append(eachAverage)
        allExchange.append(eachWaterChange)
    return allExchange
        
        
        
    

def writeExchangeToFile(self,outputFile):
    outputFile.write('centerNo 1 2 3 4 5 6 7 8 9 10 11 12  average (ns) \n')
    for waterNo in range(len(self)):
        eachExchange=self[waterNo]
        strExchange=''
        for i in eachExchange:
            i=str(i)
            strExchange=strExchange+'   '+i
##        eachExchange=str(eachExchange)[1:-1]
        strExchange= str(waterNo)+' '+strExchange+'\n'
        outputFile.write(strExchange)
    outputFile.close()
    
    
    

def main():
    centreCoord=getCentreCoord(centreFile)
    print 'there were in total',len(centreCoord),'coorditates'
    waters=getWaterInforWithCentre(waterFile)
    print 'there were in total',len(waters),'waters'
    waterInforCentre=orangeWaterInforWithCentre(waters,centreCoord)
    print len(waterInforCentre)
    waterExchangeRate=calculateExchangeRate(waterInforCentre,framePerNs)
    print waterExchangeRate
    waterExchangeRateAverage=calculateAverageExchangeRate(waterExchangeRate)
    writeExchangeToFile(waterExchangeRateAverage,outputFile)

    
    
if __name__=='__main__':
    main()
                
                
                
                
                
            
            
        
        
        
        
    

