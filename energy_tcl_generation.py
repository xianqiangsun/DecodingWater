import numpy
import sys
sys.path.append('/home/x/xiansu/pfs/program/numpy/lib/python2.6/site-packages')
sys.path.append('/bubo/home/h8/xians/glob/program/python_pkg/lib64/python2.6/site-packages')
from numpy import cos, sin, pi
#from Numeric import *

centreFile=open('optimizedCentre.txt','r')
waterFile=open('waterInforByCentre.txt','r')
allH2OInforFile=open('WW_allCentre_H2O.txt','r')


frameNo=6000

k=1.380648813*(10**(-23))   ## The unit of the boltzmann constant is J/K.
weiH2O=18.0154
mol=6.02214179*(10**23)
pw=0.0331725               ## The unit of this is No. of molecules in per A2.

WProt=open('water_protein.tcl','w')
WPOPC=open('water_POPC.tcl','w')
WWater=open('water_water.tcl','w')
WSys=open('water_system.tcl','w')

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

def getOccuWaterID(self):
    waterID=[]
    for waterInforEachCentre in self:
        waterIDEachCentre=[]
        for eachWater in waterInforEachCentre:
            eachWaterInfor=[]
            segname=eachWater[-3]
            resid=eachWater[2]
            eachWaterInfor.append(segname)
            eachWaterInfor.append(resid)
            if eachWaterInfor not in waterIDEachCentre:
                waterIDEachCentre.append(eachWaterInfor)
        waterID.append(waterIDEachCentre)
    return waterID
                
##getOccuWaterID(self) reads the output of orangeWaterInforWithCentre. According to each water
##information with centre infor, The segname and resid of waters according to each centre  is saved
##as list. the output look like[[[centre0watrer0segname, centre0water0resid],[...]...]]

def writeWaterProteinEnergyFile(self,outputFile,frameNo,waterInforCenter,threshold):
    outputFile.write('package require namdenergy')
    outputFile.write('set sel2 [atomselect top protein]'+'\n')
    outputFile.write('set namdroute "/home/x/xiansu/pfs/program/NAMD_CVS-2012-04-19_Linux-x86_64-multicore/namd2"'+'\n')
    for No in range(len(self)):
        waterNos=len(waterInforCenter[No])
        print waterNos/float(frameNo)
        if waterNos/float(frameNo)>=threshold:
            waterInforEachCentre=self[No]
            for eachWater in waterInforEachCentre:
##                print eachWater
                segname=eachWater[0]
                resid=eachWater[1]
                sel1='set sel1 [atomselect top "segname '+segname+' and resid '+resid+'"]'+'\n'
##                print sel1
                outputFile.write(sel1)
                commandLine='namdenergy -vdw -elec -nonb -sel $sel1 $sel2 -extsys namd.xsc -exe $namdroute -ofile "WaterProteinEnergy/'
                fileName=segname+'_'+resid+'_protein.dat"'+'\n'
                commandLine=commandLine+fileName
##                print commandLine
                outputFile.write(commandLine)
    outputFile.close()

def writeWaterPOPCEnergyFile(self,outputFile,frameNo,waterInforCenter,threshold):
    outputFile.write('package require namdenergy')
    outputFile.write('set sel2 [atomselect top "resname POPC"]'+'\n')
    outputFile.write('set namdroute "/home/x/xiansu/pfs/program/NAMD_CVS-2012-04-19_Linux-x86_64-multicore/namd2"'+'\n')
    for No in range(len(self)):
        waterNos=len(waterInforCenter[No])
        print waterNos/float(frameNo)
        if waterNos/float(frameNo)>=threshold:
            waterInforEachCentre=self[No]
            for eachWater in waterInforEachCentre:
##                print eachWater
                segname=eachWater[0]
                resid=eachWater[1]
                sel1='set sel1 [atomselect top "segname '+segname+' and resid '+resid+'"]'+'\n'
##                print sel1
                outputFile.write(sel1)
                commandLine='namdenergy -vdw -elec -nonb -sel $sel1 $sel2 -extsys namd.xsc -exe $namdroute -ofile "WaterPOPCEnergy/'
                fileName=segname+'_'+resid+'_POPC.dat"'+'\n'
                commandLine=commandLine+fileName
##                print commandLine
                outputFile.write(commandLine)
    outputFile.close()

def writeWaterWaterEnergyFile(self,outputFile,frameNo,waterInforCenter,threshold):
##    outputFile.write('set sel2 [atomselect top waters]'+'\n')
    outputFile.write('package require namdenergy')
    outputFile.write('set namdroute "/home/x/xiansu/pfs/program/NAMD_CVS-2012-04-19_Linux-x86_64-multicore/namd2"'+'\n')
    for No in range(len(self)):
        waterNos=len(waterInforCenter[No])
        print waterNos/float(frameNo)
        if waterNos/float(frameNo)>=threshold:
            waterInforEachCentre=self[No]
            for eachWater in waterInforEachCentre:
##                print eachWater
                segname=eachWater[0]
                resid=eachWater[1]
                sel1='set sel1 [atomselect top "segname '+segname+' and resid '+resid+'"]'+'\n'
                sel2='set sel2 [atomselect top "waters and not (segname '+segname+' and resid '+resid+')"]'+'\n'
##                print sel1
                outputFile.write(sel1)
                outputFile.write(sel2)
                commandLine='namdenergy -vdw -elec -nonb -sel $sel1 $sel2 -extsys namd.xsc -exe $namdroute -ofile "WaterWaterEnergy/'
                fileName=segname+'_'+resid+'_waters.dat"'+'\n'
                commandLine=commandLine+fileName
##                print commandLine
                outputFile.write(commandLine)
    outputFile.close()

##these three functions used to write out the tcl input script for energy calculation used by NAMD.
##The input incluse: self: the output of waterSegnameResID, which includes water segname and water
##resid according to each centre. the sencond input is the filename to write the tcl script. the
##the fourth one is the output of 'orangeWaterInforWithCentre' and the last one is the threshould value
##of which water centre to be considered.
def writeWaterSystemEnergyFile(self,outputFile,frameNo,waterInforCenter,threshold):
##    outputFile.write('set sel2 [atomselect top waters]'+'\n')
    outputFile.write('package require namdenergy')
    for No in range(len(self)):
        waterNos=len(waterInforCenter[No])
        print waterNos/float(frameNo)
        if waterNos/float(frameNo)>=threshold:
            waterInforEachCentre=self[No]
            for eachWater in waterInforEachCentre:
##                print eachWater
                segname=eachWater[0]
                resid=eachWater[1]
                sel1='set sel1 [atomselect top "segname '+segname+' and resid '+resid+'"]'+'\n'
                sel2='set sel2 [atomselect top "all and not (segname '+segname+' and resid '+resid+')"]'+'\n'
##                print sel1
                outputFile.write(sel1)
                outputFile.write(sel2)
                commandLine='namdenergy -vdw -elec -nonb -sel $sel1 $sel2 -extsys namd.xsc -ofile "WaterSystemEnergy/'
                fileName=segname+'_'+resid+'_system.dat"'+'\n'
                commandLine=commandLine+fileName
##                print commandLine
                outputFile.write(commandLine)
    outputFile.close()

##these three functions used to write out the tcl input script for energy calculation used by NAMD.
##The input incluse: self: the output of waterSegnameResID, which includes water segname and water
##resid according to each centre. the sencond input is the filename to write the tcl script. the
##the fourth one is the output of 'orangeWaterInforWithCentre' and the last one is the threshould value
##of which water centre to be considered.

            
 
centreCoord=getCentreCoord(centreFile)
print 'there were in total',len(centreCoord),'coorditates'
waters=getWaterInforWithCentre(waterFile)
print 'there were in total',len(waters),'waters'
waterInforCentre=orangeWaterInforWithCentre(waters,centreCoord)
print len(waterInforCentre)
print waterInforCentre[0][0]
waterSegnameResID=getOccuWaterID(waterInforCentre)
print waterSegnameResID[1]
writeWaterProteinEnergyFile(waterSegnameResID,WProt,frameNo,waterInforCentre,0.8)
writeWaterPOPCEnergyFile(waterSegnameResID,WPOPC,frameNo,waterInforCentre,0.8)
writeWaterWaterEnergyFile(waterSegnameResID,WWater,frameNo,waterInforCentre,0.8)
writeWaterSystemEnergyFile(waterSegnameResID,WSys,frameNo,waterInforCentre,0.8)


