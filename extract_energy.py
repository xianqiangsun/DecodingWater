import numpy
import sys
sys.path.append('/home/x/xiansu/pfs/program/numpy/lib/python2.6/site-packages')

from Numeric import *

centreFile=open('optimizedCentre.txt','r')
waterFile=open('waterInforByCentre.txt','r')
allH2OInforFile=open('WW_allCentre_H2O.txt','r')

interactionEnergy=open('interaction_Energy.dat','w')

frameNo=6000

k=1.380648813*(10**(-23))   ## The unit of the boltzmann constant is J/K.
weiH2O=18.0154
mol=6.02214179*(10**23)
pw=0.0331725               ## The unit of this is No. of molecules in per A2.

def getCentreCoord(self):
    centreCoord=[]
    for centre in self:
        centre=centre.split()
        centreFloat=[]
        for coord in centre:
            coord =float(coord)
            centreFloat.append(coord)
            
        centreCoord.append(centreFloat)
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

def getWaterEnergy(self,frameNo,threshold,whichPart,filePartName):
    energy=[]
    for waterInforEachCentre in self:
        energyEachCentre=[]
        waterNos=len(waterInforEachCentre)
        if waterNos/float(frameNo)<threshold:
            energy.append(energyEachCentre)
        elif waterNos/float(frameNo)>=threshold:
            folderName='Water'+str(whichPart)+'Energy/'
            for eachWater in waterInforEachCentre:
                segname=eachWater[-3]
                resid=eachWater[2]
                frameNo=int(eachWater[-2])
                fileName=folderName+segname+'_'+resid+'_'+filePartName+'.dat'
##                print fileName
                energyFile=open(fileName,'r')
                splitedEnergy=[]
                for line in energyFile:
                    line=line.split()
                    splitedEnergy.append(line)
                energyFile.close()
                eachEnergy=float(splitedEnergy[frameNo][-1])
##                print eachEnergy
                energyEachCentre.append(eachEnergy)
##                print energyEachCentre
            energy.append(energyEachCentre)
            print len(energy)
    return energy

##getWaterEnergy  is the function used to get the energy for each centre. the input includes
##self: the organized water information according to each centre. It can be the output of
##orangeWaterInforWithCentre;frameNo is the number of frames for the calculation. threshold
##is the threshould value of the smallest occuparation rate for calculation. This should be
##consistant with the 'tcl generation script'. whichPart means which part is used for the calculation
## in my script, It should be 'Water','Protein' or'POPC'. the filePartName should be 'waters',
##or 'protein',or 'POPC'. The whichPart and filePartName is determined by which part is used for the
##calculation. the output file is the interaction energy of each centre towards the part selected.

def printEnergy(WProteinEnergy,WWEnergy,WPOPCEnergy,waterInforCentre,outputFile):
    outputFile.write('CentreNO    Occupy    WproteinEnergy    WWEnergy    WPOPCEnergy    totalEnergy'+'\n')
    for CentreNo in range(len(waterInforCentre)):
        waterInforEachCentre=waterInforCentre[CentreNo]
        Occupy=len(waterInforEachCentre)
        print len(WProteinEnergy)
##        print WProteinEnergy
        eachWProteinEnergy=WProteinEnergy[CentreNo]
        eachWWEnergy=WWEnergy[CentreNo]
        eachWPOPCEnergy=WPOPCEnergy[CentreNo]
        if eachWProteinEnergy==[]:
            line=str(CentreNo)+'    '+str(Occupy)+'    '+'0.000'+'    '+'0.000'+'    '+'0.000'+'    '+'0.0000'+'\n'
            outputFile.write(line)
        elif eachWProteinEnergy!=[]:
            totalWProteinEnergy=sum(eachWProteinEnergy)/(float(Occupy))
            totalWWEnergy=sum(eachWWEnergy)/(float(Occupy))
            totalWPOPCEnergy=sum(eachWPOPCEnergy)/(float(Occupy))
            totalEnergy=totalWProteinEnergy+totalWWEnergy+totalWPOPCEnergy
            line=str(CentreNo)+'    '+str(Occupy)+'    '+str(totalWProteinEnergy)+'    '+str(totalWWEnergy)+'    '+str(totalWPOPCEnergy)+'    '+str(totalEnergy)+'\n'
            outputFile.write(line)
            print line
    outputFile.close()

centreCoord=getCentreCoord(centreFile)
print 'there were in total',len(centreCoord),'coorditates'
waters=getWaterInforWithCentre(waterFile)
print 'there were in total',len(waters),'waters'
waterInforCentre=orangeWaterInforWithCentre(waters,centreCoord)
print len(waterInforCentre)
##print waterInforCentre[0][0]

WaterProteinEnergy=getWaterEnergy(waterInforCentre,frameNo,0.8,'Protein','protein')
WaterWaterEnergy=getWaterEnergy(waterInforCentre,frameNo,0.8,'Water','waters')
WaterPOPCEnergy=getWaterEnergy(waterInforCentre,frameNo,0.8,'POPC','POPC')



print len(WaterProteinEnergy)

printEnergy(WaterProteinEnergy,WaterWaterEnergy,WaterPOPCEnergy,waterInforCentre,interactionEnergy)

                    
                    
                    
                
                
                
            
        
        
    
