############################################################################
##This file is used to calcululate the high density centre with more strict method
##. Then the H2O coordinates at each water Centre were obtained for the orientation
##calculation(eular calculation).
##
##
##Xianqiang Sun
##TheoChem&Bio
##KTH
##2012-05-16
###########################################################################

import numpy
from Numeric import *
from datetime import datetime
from MDAnalysis import Universe, Writer
from MDAnalysis.analysis.distances import distance_array
import MDAnalysis


DCD='restr_md_3.dcd'
PSF='ionized.psf'

centreFile=open('optimizedCentre.txt','r')
waterFile=open('waterInforByCentre.txt','r')
H2OFile=open('SW_allCentre_H2O.txt','w')

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

def removeLessOcuCentre(self,waterInforAccCen,frameNo):
    highOcuCentre=[]
    for i in range(len(self)):
        ocupyNo=len(waterInforAccCen[i])

        if (float(ocupyNo)/float(frameNo))>=0.95:
            highOcuCentre.append(self[i])
    return highOcuCentre

def removeLessOcuWaterInfor(self,waterInforAccCen,frameNo):
    highOcuWaterInfor=[]
    for i in range(len(self)):
        ocupyNo=len(waterInforAccCen[i])
        
        if (float(ocupyNo)/float(frameNo))>=0.95:
            highOcuWaterInfor.append(waterInforAccCen[i])
    return highOcuWaterInfor
##removeLessOcuCentre is used to obtain the high occupied centre with the ratio>=0.95. The input
##were self: all the centre coordidated; WaterInforAccCen: oranged water information according to
##the centre; frameNo is the total Number of framed you want to calculated.
##These highOcuCentre were saved as list of coordidates. The highOcuWaterInfor were saved according
##to the centre highOcuCentre with the format of [[[waterInforInCentre1],[waterInforInCentre1]]..].

def extractWaterCoor(self):
    waterCoor=[]
    for waterSet in self:
        waterCoorSet=[]
        for waterInfor in waterSet:
            eachWaterCoor=[]
            eachWaterCoorStr=waterInfor[3:-3]
            for i in eachWaterCoorStr:
                i=float(i)
                eachWaterCoor.append(i)
            waterCoorSet.append(eachWaterCoor)
        waterCoor.append(waterCoorSet)
    return waterCoor
##ExtractWaterCoor use the output of  function of removeLessOcuWaterInfor.The waterCoord were extracted
##according to each high occupied water coordidate centre. The format of the output is[[waterCoor,waterCoor]
##...]
            


def getWaterCoorWithH(self,centre,psf,dcd,outputFile):
    rho=Universe(psf,dcd)
    H2OCoordinate=[]
    no=0
    title='resname'+'    '+'atomid'+'    '+'resnumber'+'    X    Y     Z   '+'   '+'segname'+'  '+'frameNo'+'   '+'centreNo'+'\n'
    outputFile.write(title)
    for oxygenInforSet in self:
        H2OCoordinateSet=[]
        
        
        print 'There were',len(oxygenInforSet),'waters in the'
        for oxygenInfor in oxygenInforSet:
##            no1+=1
##            print no1
            frameNo=oxygenInfor[-2]
            frameNo=int(frameNo)-1
            segName=oxygenInfor[-3]
            resNumber=oxygenInfor[2]
            frame=rho.trajectory[frameNo]
            infor='segid '+segName+' and resid '+resNumber
            selected=rho.selectAtoms(infor)
            atomID=[]
            for atoms in selected.atoms:
                ID=str(atoms).split()[2][:-1]
                atomID.append(ID)
            selectedResId=selected.resids()
            selectedResNa=selected.resnames()
            coordsOH1H2=selected.coordinates()
            for i in range(3):
                atomInfor=str(selectedResNa[0])+'    '+str(atomID[i])+'    '+str(resNumber)+'    '+str(coordsOH1H2[i])[1:-1]+'   '+segName+'    '+str(frameNo)+'    '+str(no)+'\n'
                outputFile.write(atomInfor)
            H2OCoordinateSet.append(coordsOH1H2)

        no+=1

        H2OCoordinate.append(H2OCoordinateSet)
        print no,'is finished'
    outputFile.close()
    return H2OCoordinate

##getWaterCoorWithH is uses to extract the water coordinates at each centre. The oxygen atoms has extracted from the trajectroy
##Therefore, we just need to extract these H2O atom coordinates according to these information. The input of this function
##includes oxygen information which includeds the frameNo information for the calculation. You should remind that the frameNo
##changeing. Because the first frame is numbered as 1 in water infor. However, In the H2O coordinates, the frameNo is change
##to start with 0. Self: waterInfor which is saved as sets according to centre.Centre: The centres for the calculation. psf and
##dcd are the output files for the calculation. outputFile is the file for write in the calculaiton.

        

 
centreCoord=getCentreCoord(centreFile)
print 'there were in total',len(centreCoord),'coorditates'
waters=getWaterInforWithCentre(waterFile)
print 'there were in total',len(waters),'waters'
waterInforCentre=orangeWaterInforWithCentre(waters,centreCoord)
print len(waterInforCentre)
removedCentre=removeLessOcuCentre(centreCoord,waterInforCentre,965)

print len(removedCentre)
removedWaterInfor=removeLessOcuWaterInfor(centreCoord,waterInforCentre,965)
print len(removedWaterInfor)
for i in removedWaterInfor:
    print len(i)

removedWaterCoor=extractWaterCoor(removedWaterInfor)
print len(removedWaterCoor)
for i in removedWaterCoor:
    print len(i)

H2OCoordinates=getWaterCoorWithH(removedWaterInfor,removedCentre,PSF,DCD,H2OFile)

print len(H2OCoordinates)
for i in H2OCoordinates:
    print len(i)

    
            
        

            
        
    



        
    
