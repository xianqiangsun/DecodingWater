############################################################################
##This file is used to calculate the WW orientation entropy with respect to
##the each water pair. Firstly, each H2O coordinates were obtained for calculation
##Then, the angle between each H2O pair is obtained (FrameInfor should be take care
##)
##
##
##
##Xianqiang Sun
##TheoChem&Bio
##KTH
##2012-05-28
###########################################################################

#symmetry for the dictionary

import numpy
import sys
import random
sys.path.append('/home/x/xiansu/pfs/program/numpy/lib/python2.6/site-packages')

sys.path.append('/bubo/home/h8/xians/glob/program/python_pkg/lib64/python2.6/site-packages')
from numpy import cos,sin,arccos
from numpy import pi
#from Numeric import *
from datetime import datetime


centreFile=open('optimizedCentre.txt','r')
waterFile=open('waterInforByCentre.txt','r')
allH2OInforFile=open('WW_allCentre_H2O.txt','r')
##conservedH2OInforFile=open('SW_conservedCentre_H2O.txt','r')
angleFile=open('WW_angle.txt','w')
rdfFile=open('rdf_HO_HO.xvg','r')
SorOutput=open('WW_Orien_R5.dat','w')

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

def orangeH2OInfor(self):
    allH2OInfor=[]
    for H2OSet in self:
        H2OSetInfor=[]
        eachH2O=[]
        no=1
        for H2O in H2OSet:
            
            if no/3.0==1.0:
                eachH2O.append(H2O)
                
                H2OSetInfor.append(eachH2O)

                eachH2O=[]
                no=1
            else:
                eachH2O.append(H2O)
                no+=1
##        for i in H2OSetInfor:
##            print i
        print 'There were',len(H2OSetInfor), 'in H2OSetInfor!'
        allH2OInfor.append(H2OSetInfor)
    return allH2OInfor
##orangeH2OInfor can be used to organize a sets of waters which is saved as list H2O informations to sets,
##As there were three atoms in the H2O, therefore, each H2O can be saved as a list. Then, the H2O atom
##lists can be saved as lists. At last the list can be saved according to each centre.
##The output looks like [[[[AtomInfor]...WaterInfor]...CentreInfor]...All H2O]

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

def getH2OCoord(self):
    allH2OCoord=[]
    for H2OInforAccCentre in self:
##        print len(H2OInforAccCentre)
        H2OCoordAccCentre=[]

        for H2OInfor in H2OInforAccCentre:

            H2OCoord=[]

            for H2OAtom in H2OInfor:
                

                H2OAtomCoord=H2OAtom[3:-3]
                for number in range(len(H2OAtomCoord)):
                    H2OAtomCoord[number]=float(H2OAtomCoord[number])
                H2OCoord.append(H2OAtomCoord)
            H2OCoordAccCentre.append(H2OCoord)
        allH2OCoord.append(H2OCoordAccCentre)
    
                    
        print 'The coordinate length in each centre is :',len(H2OCoordAccCentre)
    return allH2OCoord
        
##getH2OCoord is used to extract all the H2O coordidates and save each H2O coordinates
##as list. Then H2O coordinates were saved according to each centre. At last, all
##the H2O coords were saved as a list.

def H2OInforNearbyHighCentre(self,allCentre,conservedH2OInfor,nearbyWaterInfor):
    H2OInforAndAround=[]
    for conservedCentreNo in range(len(self)):
        eachH2OAndAround=[]
        eachConservedCentre=self[conservedCentreNo]
        eachConservedCentreArray=numpy.array(eachConservedCentre)
        eachH2OAndAround.append(conservedH2OInfor[conservedCentreNo])
        for allCentreNo in range(len(allCentre)):
            aroundCentre=allCentre[allCentreNo]
            aroundCentreArray=numpy.array(aroundCentre)
            vectorCC=eachConservedCentreArray-aroundCentreArray
            vectorCCModulus=numpy.sqrt((vectorCC*vectorCC).sum())
            if 5.0>=vectorCCModulus>=1.5:
                
                eachH2OAndAround.append(nearbyWaterInfor[allCentreNo])
        if len(eachH2OAndAround)>1:
            H2OInforAndAround.append(eachH2OAndAround)
    print len(H2OInforAndAround)

    return H2OInforAndAround

##H2OInforNearbyHighCentre can be used to find centre around each high occupied centre.The input included
##several information to be calculated. self: the input of the high occupied centre. allcentre: the
##centres for all the optimized centres. conservedH2OInfor inculdes all the H2O infor according to the
##high occupied centre. The nearbyWaterInfor includes all the h2O infors consitant with allCentre
##the output of this functions includes each high occupied H2O infor
##and the H2O infor around each high occupied centre.[[[highOccupiedH2OInfor],[aroundH2OInfor2],[around
##H2OInfor2]],[...]...]. Moreover, the output centre is consistant with centre of findNearByCentre output.


def findNearByWater(self,allCentre,highWaterInfor,waterInfor):
    centreAndAround=[]
    WWCentreNo=[]

    for centreNo in range(len(self)):
        eachCentreAndAround=[]

        eachCentre=self[centreNo]
        eachCentre=numpy.array(eachCentre)
        eachCentreAndAround.append(highWaterInfor[centreNo])
        no=0
        for allCentreNo in range(len(allCentre)):
            aroundCentre=allCentre[allCentreNo]
            aroundCentre=numpy.array(aroundCentre)
            vectorCC=eachCentre-aroundCentre
            vectorCCModulus=numpy.sqrt((vectorCC*vectorCC).sum())
            
            if 5.0>=vectorCCModulus>=1.5:
                no+=1
                eachCentreAndAround.append(waterInfor[allCentreNo])
##        print no,'waters arround centre',centreNo

        if len(eachCentreAndAround)>1:
            WWCentreNo.append(centreNo)
            print len(eachCentreAndAround)
            

            centreAndAround.append(eachCentreAndAround)
    return centreAndAround

##findNearByWater can be used to find waters around each high occupied centre. The input included
##several information to be calculated. self: the input of the high occupied centre. allcentre: the
##centres for all the optimized centres. The highwaterinfor includes the water informations for all
##the high conserced waters. the output of this functions includes each high occupied water centre
##and the waters around each high occupied centre.[[[highOccupiedCentreInfor1],[aroundWater1],[around
##Water1]],[...]...]

def findNearbyCentre(self,allCentre):
    centreNearby=[]

    for centreNo in range(len(self)):
        eachCentreAndAround=[]

        eachCentre=self[centreNo]
        eachCentreArray=numpy.array(eachCentre)
        eachCentreAndAround.append(eachCentre)
        no=0
        for allCentreNo in range(len(allCentre)):
            aroundCentre=allCentre[allCentreNo]
            aroundCentreArray=numpy.array(aroundCentre)
            vectorCC=eachCentreArray-aroundCentreArray
            vectorCCModulus=numpy.sqrt((vectorCC*vectorCC).sum())
            
            if 5.0>=vectorCCModulus>=1.5:
                no+=1
                eachCentreAndAround.append(aroundCentre)
##        print no,'waters arround centre','centreNo'
        if len(eachCentreAndAround)>1:
            
            centreNearby.append(eachCentreAndAround)
##    print centreNearby
    return centreNearby

##findNearByCentre can be used to find centre around each high occupied centre. The input included
##several information to be calculated. self: the input of the high occupied centre. allcentre: the
##centres for all the optimized centres. the output of this functions includes each high occupied water centre
##and the waters centre around each high occupied centre.[[[highOccupiedCentre1],[aroundCentre2],[around
##centre2]],[...]...]

def findNearbyCentreNo(self,allCentre):
    
    WWCentreNo=[]
    for centreNo in range(len(self)):
        eachCentreAndAround=[]

        eachCentre=self[centreNo]
        eachCentreArray=numpy.array(eachCentre)
        eachCentreAndAround.append(eachCentre)
        no=0
        for allCentreNo in range(len(allCentre)):
            aroundCentre=allCentre[allCentreNo]
            aroundCentreArray=numpy.array(aroundCentre)
            vectorCC=eachCentreArray-aroundCentreArray
            vectorCCModulus=numpy.sqrt((vectorCC*vectorCC).sum())
            
            if 5.0>=vectorCCModulus>=1.5:
                no+=1
                eachCentreAndAround.append(aroundCentre)
##        print no,'waters arround centre','centreNo'
        if len(eachCentreAndAround)>1:
            
            WWCentreNo.append(centreNo)
    print WWCentreNo
    return WWCentreNo

##findNearByCentreNo can be used to find centre around each high occupied centre. The input included
##several information to be calculated. self: the input of the high occupied centre. allcentre: the
##centres for all the optimized centres. the output of this functions includes each high occupied water centre
##and the waters centre around each high occupied centre.[[[highOccupiedCentre1],[aroundCentre2],[around
##centre2]],[...]...]

def findNearbyCentreAndAroundNo(self,allCentre):
    
    WWCentreNo=[]
    for centreNo in range(len(self)):
        eachCentreNoAndAround=[]

        eachCentre=self[centreNo]
        eachCentreArray=numpy.array(eachCentre)
        eachCentreNoAndAround.append(centreNo)
        no=0
        for allCentreNo in range(len(allCentre)):
            aroundCentre=allCentre[allCentreNo]
            aroundCentreArray=numpy.array(aroundCentre)
            vectorCC=eachCentreArray-aroundCentreArray
            vectorCCModulus=numpy.sqrt((vectorCC*vectorCC).sum())
            
            if 5.0>=vectorCCModulus>=1.5:
                no+=1
                eachCentreNoAndAround.append(allCentreNo)
                print allCentreNo
        print no,'waters arround centre',centreNo
        
        if len(eachCentreNoAndAround)>1:
            
            WWCentreNo.append(eachCentreNoAndAround)
    print 'each CenterNo and around:',WWCentreNo
    return WWCentreNo

##findNearByCentreAndAroundNo can be used to find centre around each high occupied centre. The input included
##several information to be calculated. self: the input of the high occupied centre. allcentre: the
##centres for all the optimized centres. the output of this functions includes each high occupied water centre
##and the waters centre around each high occupied centre.[[[highOccupiedCentre1No,aroundCentreNo1,aroundCentreNo2],
##[...]...]



def extractWaterCoor(self):
    waterCoor=[]
    for waterSet in self:
        waterCoorSet=[]
        for waterInfor in waterSet:
            eachWaterCoor=[]
            eachWaterCoorStr=waterInfor[3:-3]
##            print eachWaterCoorStr
            for i in eachWaterCoorStr:
##                print i
                i=float(i)
                eachWaterCoor.append(i)
            waterCoorSet.append(eachWaterCoor)
        waterCoor.append(waterCoorSet)
    return waterCoor
##ExtractWaterCoor use the output of  function of removeLessOcuWaterInfor.The waterCoord were extracted
##according to each high occupied water coordidate centre. The format of the output is[[waterCoor,waterCoor]
##...]

def extractWaterCoorAround(self):
    waterCoordAround=[]
    for eachWaterCentre in waterInforAroundCentre:
        
##        print len(eachWaterCentre)
        
        
        eachWaterCoordAround=extractWaterCoor(eachWaterCentre)

        
        
        waterCoordAround.append(eachWaterCoordAround)
##    print len(waterCoordAround)

    return waterCoordAround

##ExtractWaterCoorAround use mainly use the function of ExtractWaterCoor to extract the water coordinates
##Around each centre. This function read the output of 'findNearByWater' and then the save each waterCoordnate
##according to each high occupied water coordidate centre. The highly accupied centre was saved as the first
##item in each sublist. The waters around this centre were saved as the following items
##The format of the output is[[[centreWaterCoorList],[waterCoorAroundCentre],[waterCoorAroundCentre]..]..]



def calculateGR(self,centre,frameNo):
    gr=[]

    waterDensPerA2=0.0331725
    
    constant1=(4/3.0)*pi
    for centreNo in range(len(centre)):
        grNo=range(24)
        waterCoorSet=self[centreNo]
##        frameNo=len(waterCoorSet)
        centreCoor=centre[centreNo]
        centreCoor=numpy.array(centreCoor)
        for no in range(len(grNo)):
            grNo[no]=0
##        print grNo    
        for waterCoor in waterCoorSet:
##            print waterCoor, centreCoor
            waterCoor=numpy.array(waterCoor)
            dist=numpy.linalg.norm(waterCoor-centreCoor)
##            print dist
            for number in range(24):
                nextNumber=number+1
                if nextNumber*0.05>dist>=number*0.05:
                    grNo[number]+=1
##                    print waterCoor, centreCoor, dist,   number
                elif dist==24*0.05:
                    grNo[23]+=1
##        print grNo
        for number in range(len(grNo)):
            prevNumber=number+1
            grNo[number]=round(grNo[number]/(constant1*(((prevNumber*0.05)**3)-((number*0.05)**3)))/waterDensPerA2/frameNo,3)
##        print grNo    
        gr.append(grNo)

    return gr

##calculateGR is used to calculate the g(r) distribution of each coordidate centre
##the input include : (1),self: the water coordidate sets according to the centre
##(2)Centre:each centre for calculation. (3) frameNo: total frames for the
##calculation. The frameNo information should be included to calculate the
##water density at each centre. The output of this function is list of g(r)
##According to each centre. Moreover, as we have to calculate the g(r) according
##to the distance between the centre and the oxygen coordidated in each frame.
##Each g(r) were saved as list. the form of the output looks like
##[[g(0~0.1),g(0.1~0.2)...g(1.1~1.2)]......]

def calculateGRTheta(self,centre,frameNo):
    grTheta=[]
##    thetaDensity=frameNo/20.0
    reference=numpy.array([0,0,1])
    referenceModulus=numpy.sqrt((reference*reference).sum())
    for centreNo in range(len(centre)):
        grThetaNo=range(20)
        waterCoorSet=self[centreNo]
        centreCoor=centre[centreNo]

        centreCoor=numpy.array(centreCoor)
#        frameNo=len(waterCoorSet)
        thetaDensity=frameNo/20.0
        for no in range(len(grThetaNo)):
            grThetaNo[no]=0
        for waterCoor in waterCoorSet:
##            print waterCoor, centreCoor
            waterCoor=numpy.array(waterCoor)
            orientation=waterCoor-centreCoor
            orientationModulus=numpy.sqrt((orientation*orientation).sum())
            dot=numpy.dot(orientation,reference)
            cosAngle=dot/referenceModulus/orientationModulus
            angle=numpy.arccos(cosAngle)
##            print 'the angle of Phi is',angle
            for number in range(20):
                nextNumber=number+1
                if nextNumber*pi/20.0>angle>=number*pi/20.0:
##                    print 'the angle of Phi is',waterCoor, centreCoor, angle,   number
                    grThetaNo[number]+=1
                elif angle==pi:
                    grThetaNo[19]+=1
##        print 'the angle distri bution',grThetaNo    
        for number in range(len(grThetaNo)):
            nextNumber=number+1
            grThetaNo[number]=(grThetaNo[number]/(cos(number*pi/20.0)-cos(nextNumber*pi/20.0)))/(frameNo/2.0)
        grTheta.append(grThetaNo)
##    print grTheta
    return grTheta
####????????????????????????????????????????????Maybe A question here        
##calculateGRTheta is used to calculate the g(theta) distribution of each coordidate centre
##the input include : (1),self: the water coordidate sets according to the centre
##(2)Centre:each centre for calculation. (3) frameNo: total frames for the
##calculation. The frameNo information should be included to calculate the
##water density at each centre. The output of this function is list of g(theta)
##According to each centre. Moreover, as we have to calculate the g(theta) according
##to the reference oritentation[0, 0,1 ] and the oritentation between the centre and the oxygen coordidated
##in each frame.Each g(theta) were saved as list. the form of the output looks like
##[[g(0~*pi/20.0),g(pi/20.0~2*pi/20.0)...g(pi~19*pi/20.0)]......]

def calculateGRPhi(self,centre,frameNo):
    grPhi=[]
##    phiDensity=frameNo/20.0
    referenceZ=numpy.array([0,0,1])
    referenceZModulus=numpy.sqrt((referenceZ*referenceZ).sum())
    referenceX=numpy.array([1,0,0])
    referenceXModulus=numpy.sqrt((referenceX*referenceX).sum())
    for centreNo in range(len(centre)):
        grPhiNo=range(40)
        waterCoorSet=self[centreNo]
        centreCoor=centre[centreNo]
        
        centreCoorArray=numpy.array(centreCoor)
#        frameNo=len(waterCoorSet)
        phiDensity=frameNo/40.0
        for no in range(len(grPhiNo)):
            grPhiNo[no]=0
        for waterCoor in waterCoorSet:
            waterCoorArray=numpy.array(waterCoor)
            orientation=waterCoorArray-centreCoorArray
            orientationModulus=numpy.sqrt((orientation*orientation).sum())
            dot1=numpy.dot(orientation,referenceZ)
            cosAngle1=dot1/referenceZModulus/orientationModulus
            zReflection=orientationModulus*cosAngle1

            waterCoorCopy=[]
            waterCoorCopy.extend(waterCoor)
            waterCoorCopy[2]=waterCoorCopy[2]-zReflection
            newWaterCoorArray=numpy.array(waterCoorCopy)
            
            newOrientation=newWaterCoorArray-centreCoorArray
            newOrientationModulus=numpy.sqrt((newOrientation*newOrientation).sum())
            dot2=numpy.dot(newOrientation,referenceX)
            cosAngle2=dot2/referenceXModulus/newOrientationModulus
            angle=numpy.arccos(cosAngle2)
            if newOrientation[1]<0:
                angle=2*pi-angle
                

                

            for number in range(40):
                nextNumber=number+1
                if nextNumber*pi/20.0>angle>=number*pi/20.0:
                    grPhiNo[number]+=1
                elif angle==2*pi:
                    grPhiNo[39]+=1
        print 'the angle distri bution grPhiNo',centreNo, grPhiNo

        for number in range(len(grPhiNo)):
            grPhiNo[number]=grPhiNo[number]/phiDensity
        grPhi.append(grPhiNo)
        print 'the angle distri bution grPhiNo',centreNo, grPhiNo
##    print grPhi
    return grPhi
####????????????????????????????????????????????Maybe A question here        
##calculateGRPhi is used to calculate the g(Phi) distribution of each coordidate centre
##the input include : (1),self: the water coordidate sets according to the centre
##(2)Centre:each centre for calculation. (3) frameNo: total frames for the
##calculation. The frameNo information should be included to calculate the
##water density at each centre. The output of this function is list of g(Phi)
##According to each centre.  We first calculate the cos(theta) according to [1,0,0]. According to cos(theta)
##The reflection of water Coordidate on XOY plane was found. Therefore, cos(theta) were obtained  according to
##the reference oritentation[0,0,1] and the oritentation between the centre and reflection on XOY
##in each frame.Each g(Phi) were saved as list. the form of the output looks like
##[[g(0~*pi/20.0),g(pi/20.0~2*pi/20.0)...g(pi~19*pi/20.0)]......]

def calculateGReachCentre(self,nearbyCentre,frameNo):
    grCentreAndAround=[]
    for eachWaterNo in range(len(self)):
        eachCentre=nearbyCentre[eachWaterNo]
        eachWater=self[eachWaterNo]
        eachGr=calculateGR(eachWater,eachCentre,frameNo)
        grCentreAndAround.append(eachGr)
    return  grCentreAndAround

##calculateGReachCentre is used to calculate the Gr of for all the centre waters and the arounded waters.
##The input included self: water coordinates which is saved as list of waters and arounded waters. This can be the
##output of function 'extractWaterCoorAround'. nearbycentre is all the centres for the calculation. The centre
##information should be consistant with the self(waterCoordinate). Morover, it can read the output of function
##'findNearbyCentre'. The output of this function looks like[[[centreGr],[aroundWatrer1 Gr]...]......]



def calculateGThetaeachCentre(self,nearbyCentre,frameNo):
    gthetaCentreAndAround=[]
    for eachWaterNo in range(len(self)):
        eachCentre=nearbyCentre[eachWaterNo]
        eachWater=self[eachWaterNo]
        eachGtheta=calculateGRTheta(eachWater,eachCentre,frameNo)
        gthetaCentreAndAround.append(eachGtheta)
    return  gthetaCentreAndAround
##calculateGThetaeachCentre is used to calculate the GTheta of for all the centre waters and the arounded waters.
##The input included self: water coordinates which is saved as list of waters and arounded waters. This can be the
##output of function 'extractWaterCoorAround'. nearbycentre is all the centres for the calculation. The centre
##information should be consistant with the self(waterCoordinate). Morover, it can read the output of function
##'findNearbyCentre'. The output of this function looks like[[[centreGTheta],[aroundWatrer1 GTheta]...]......]


def calculateGPhieachCentre(self,nearbyCentre,frameNo):
    gphiCentreAndAround=[]
    for eachWaterNo in range(len(self)):
        eachCentre=nearbyCentre[eachWaterNo]
        eachWater=self[eachWaterNo]
        eachGphi=calculateGRPhi(eachWater,eachCentre,frameNo)
        gphiCentreAndAround.append(eachGphi)
    return  gphiCentreAndAround

##calculateGPhieachCentre is used to calculate the Gphi of for all the centre waters and the arounded waters.
##The input included self: water coordinates which is saved as list of waters and arounded waters. This can be the
##output of function 'extractWaterCoorAround'. nearbycentre is all the centres for the calculation. The centre
##information should be consistant with the self(waterCoordinate). Morover, it can read the output of function
##'findNearbyCentre'. The output of this function looks like[[[centreGPhi],[aroundWatrer1 GPhi]...]......]

def getRdfOO(self):
    rdf=[]
    for line in self:
        line=line.split()
        for No in range(len(line)):
            eachElement=line[No]

            line[No]=float(eachElement)
        line[0]=line[0]*10
        rdf.append(line)
##    print rdf
    return rdf

##getRdfOO is used to get the RDF distribution of water with Oxy-Oxy paris.
##the input self is the rdf file and the output is the rdf distribution as lists.

def calculateGRInhTrans(self,rdfOO):
    print rdfOO
    waterDensPerA2=0.0331725
    constant1=(4/3.0)*pi

    grInhTrans=[]
    
    for setNo in range(len(self)):
        grInhTransSet=[0]
        centreSet=self[setNo]
##        waterCoordSet=waterCoordAround[setNo]
##        print waterCoordSet
##        print 'there were :',len(waterCoordSet),'sets in waterCoordSet'
        totalNo=len(centreSet)
        referenceCentre=centreSet[0]
        referenceCentreArray=numpy.array(referenceCentre)
        for aroundNo in range(totalNo):
            if aroundNo>=1:
                aroundCentre=centreSet[aroundNo]
                aroundCentreArray=numpy.array(aroundCentre)
                distanceVector=referenceCentreArray-aroundCentreArray
                distanceModulus=numpy.sqrt((distanceVector*distanceVector).sum())
##                print distanceModulus
                for line in rdfOO:
##                    print line
                    OOdistan=line[0]-distanceModulus
##                    print line[0]
                    if OOdistan<=0.015:
                        grAround=line[1]

##                waterCoorNo=len(waterCoordSet[aroundNo])
####                print waterCoorNo
##                print distanceModulus
##                grAround=round(waterCoorNo/(constant1*((distanceModulus+1.2)**3-(distanceModulus-1.2)**3))/waterDensPerA2/frameNo,3)
                print grAround
                grAround=numpy.log(grAround)*grAround-grAround+1
                grInhTransSet.append(grAround)
        grInhTrans.append(grInhTransSet)
    print grInhTrans
    return grInhTrans
##?????????????????????????Problem in determinning the 0.025
##calculateGRInhTrans is used to calculate the G(r,r) translocation redical distribution function used fir WW translocation entropy
##the input of this function includes self: the output of 'findNearbyCentre' as self: including all the reference centre and centers
##around reference. Therefore, it is a set of centres looks like[[[highOccupiedCentre1],[aroundCentre2],[around
##centre2]],[...]...]. The other input is the water coordinates which includes all the coordinates saved according to the self(centres)
##the output of the fromat is [[G(r1,r1)=0,G(r1,r2),G(r1,r3)...],...].

def calculateGRConstant(self,rdfOO):
    print rdfOO
    waterDensPerA2=0.0331725
    constant1=(4/3.0)*pi

    grInhTrans=[]
    
    for setNo in range(len(self)):
        grInhTransSet=[0]
        centreSet=self[setNo]
##        waterCoordSet=waterCoordAround[setNo]
##        print waterCoordSet
##        print 'there were :',len(waterCoordSet),'sets in waterCoordSet'
        totalNo=len(centreSet)
        referenceCentre=centreSet[0]
        referenceCentreArray=numpy.array(referenceCentre)
        for aroundNo in range(totalNo):
            if aroundNo>=1:
                aroundCentre=centreSet[aroundNo]
                aroundCentreArray=numpy.array(aroundCentre)
                distanceVector=referenceCentreArray-aroundCentreArray
                distanceModulus=numpy.sqrt((distanceVector*distanceVector).sum())
##                print distanceModulus
                for line in rdfOO:
##                    print line
                    OOdistan=line[0]-distanceModulus
##                    print line[0]
                    if OOdistan<=0.015:
                        grAround=line[1]

##                waterCoorNo=len(waterCoordSet[aroundNo])
####                print waterCoorNo
##                print distanceModulus
##                grAround=round(waterCoorNo/(constant1*((distanceModulus+1.2)**3-(distanceModulus-1.2)**3))/waterDensPerA2/frameNo,3)
                print grAround
##                grAround=numpy.log(grAround)*grAround-grAround+1
                grInhTransSet.append(grAround)
        grInhTrans.append(grInhTransSet)
    print grInhTrans
    return grInhTrans
##?????????????????????????Problem in determinning the 0.025
##calculateGRInhTrans is used to calculate the G(r,r) translocation redical distribution function used fir WW translocation entropy
##the input of this function includes self: the output of 'findNearbyCentre' as self: including all the reference centre and centers
##around reference. Therefore, it is a set of centres looks like[[[highOccupiedCentre1],[aroundCentre2],[around
##centre2]],[...]...]. The other input is the water coordinates which includes all the coordinates saved according to the self(centres)
##the output of the fromat is [[G(r1,r1)=0,G(r1,r2),G(r1,r3)...],...].



def getH2OEulerTheta(self,frameNo):
    referenceH1=numpy.array([0.79079641377315202, 0.61207926934631729, 0.0])
    referenceH2=numpy.array([-0.79079641377315202, 0.61207926934631729, 0.0])
    referenceZ=numpy.array([0.0,0.0,1.0])
    eulerTheta=[]
##    eulerThetaDensity=frameNo/20.0
    CentreNo=0
    for H2OCoordSet in self:
        eulerThetaNo=range(10)
#        frameNo=len(H2OCoordSet)
        eulerThetaDensity=frameNo/10.0
        CentreNo+=1
        for no in range(len(eulerThetaNo)):
            eulerThetaNo[no]=0

        for H2OCoord in H2OCoordSet:

            H1=numpy.array(H2OCoord[1])

            H2=numpy.array(H2OCoord[2])

            OT=numpy.array(H2OCoord[0])
            
            H1=H1-OT
            H2=H2-OT
            ZH2O=numpy.cross(H1,H2)
            ZH2OModulus=numpy.sqrt((ZH2O*ZH2O).sum())
            dot1=numpy.dot(referenceZ,ZH2O)
            cosAngleTheta=dot1/ZH2OModulus
            angleTheta=numpy.arccos(cosAngleTheta)
##            print angleTheta
            for number in range(10):
                nextNumber=number+1
                if nextNumber*pi/10.0>angleTheta>=number*pi/10.0:
                    eulerThetaNo[number]+=1
                elif angleTheta==pi:
                    eulerThetaNo[9]+=1
#        print 'the euler Theta No is:',CentreNo, eulerThetaNo
        for number in range(len(eulerThetaNo)):
            nextNumber=number+1
            eulerThetaNo[number]=(eulerThetaNo[number]/(cos(number*pi/10.0)-cos(nextNumber*pi/10.0)))/(frameNo/2.0)
        eulerTheta.append(eulerThetaNo)
#        print 'the euler Theta No is:',CentreNo,eulerThetaNo
    return eulerTheta

##getH2OEulerTheta is used to calculate the Euler distrubution of all the waters.
##The water coordinates were saved as lists with both Oxygen, and Hydeogen atoms.
##The coordinates looks like[[[OxygenCoord],[Hydgrogen1 coordinate],[Hydrogen2 Coordinate]....]...]
##the output of the Euler angle is very similar with other functions. the distribution is
##divided into 20 intergration spaces. The output is [[g(0~*pi/20.0),g(pi/20.0~2*pi/20.0)...g(pi~19*pi/20.0)]......]

def getH2OEulerPhi(self,frameNo):
    referenceH1=numpy.array([0.79079641377315202, 0.61207926934631729, 0.0])
    referenceH2=numpy.array([-0.79079641377315202, 0.61207926934631729, 0.0])
    referenceZ=numpy.array([0.0,0.0,1.0])
    referenceY=numpy.array([0.0,1.0,0.0])
    referenceX=numpy.array([1.0,0.0,0.0])
    revReferenceY=numpy.array([0.0,-1.0,0.0])
    eulerPhi=[]
##    eulerPhiDensity=frameNo/20.0
    CentreNo=0
    for H2OCoordSet in self:
        eulerPhiNo=range(10)
#        frameNo=len(H2OCoordSet)
        eulerPhiDensity=frameNo/10.0
        CentreNo+=1
        for no in range(len(eulerPhiNo)):
            eulerPhiNo[no]=0

        for H2OCoord in H2OCoordSet:

            H1=numpy.array(H2OCoord[1])

            H2=numpy.array(H2OCoord[2])

            OT=numpy.array(H2OCoord[0])
            
            H1=H1-OT
            H2=H2-OT
            ZH2O=numpy.cross(H1,H2)
           
            
            ZH2OOnXY=numpy.copy(ZH2O)
            ZH2OOnXY[2]=0.0
            ZH2OOnXYModulus=numpy.sqrt((ZH2OOnXY*ZH2OOnXY).sum())
            dot1=numpy.dot(revReferenceY,ZH2OOnXY)
            cosAnglePhi=dot1/ZH2OOnXYModulus
            anglePhi=numpy.arccos(cosAnglePhi)
            if ZH2OOnXY[0]<0:
                anglePhi=2*pi-anglePhi
            
            for number in range(10):
                nextNumber=number+1
                if nextNumber*pi/5.0>anglePhi>=number*pi/5.0:
                    eulerPhiNo[number]+=1
                elif anglePhi==2*pi:
                    eulerPhiNo[9]+=1
#        print 'the euler phi No is:',CentreNo, eulerPhiNo
        for number in range(len(eulerPhiNo)):
            eulerPhiNo[number]=eulerPhiNo[number]/eulerPhiDensity
        eulerPhi.append(eulerPhiNo)
#        print 'the euler phi No is:',CentreNo, eulerPhiNo
    return eulerPhi

##getH2OEulerPhi is used to calculate the Euler distrubution of all the waters.
##The water coordinates were saved as lists with both Oxygen, and Hydeogen atoms.
##The coordinates looks like[[[OxygenCoord],[Hydgrogen1 coordinate],[Hydrogen2 Coordinate]....]...]
##the output of the Euler angle is very similar with other functions. the distribution is
##divided into 20 intergration spaces. The output is [[g(0~*pi/20.0),g(pi/20.0~2*pi/20.0)...g(pi~19*pi/20.0)]......]

def getH2OEulerPsi(self,frameNo):
    referenceH1=numpy.array([0.79079641377315202, 0.61207926934631729, 0.0])
    referenceH2=numpy.array([-0.79079641377315202, 0.61207926934631729, 0.0])
    referenceZ=numpy.array([0.0,0.0,1.0])
    referenceY=numpy.array([0.0,1.0,0.0])
    revReferenceY=numpy.array([0.0,-1.0,0.0])
    eulerPsi=[]
##    eulerPsiDensity=frameNo/20.0
    CentreNo=0
    for H2OCoordSet in self:
#        frameNo=len(H2OCoordSet)
        eulerPsiDensity=frameNo/10.0
        eulerPsiNo=range(10)
        CentreNo+=1
        for no in range(len(eulerPsiNo)):
            eulerPsiNo[no]=0

        for H2OCoord in H2OCoordSet:

            H1=numpy.array(H2OCoord[1])

            H2=numpy.array(H2OCoord[2])

            OT=numpy.array(H2OCoord[0])
            
            H1=H1-OT
            H2=H2-OT
            YH2O=numpy.add(H1,H2)
            ZH2O=numpy.cross(H1,H2)
            XH2O=numpy.cross(YH2O,ZH2O)
            
            
            ZH2OOnXY=numpy.copy(ZH2O)

            ZH2OOnXY[2]=0.0
            crossLine=numpy.cross(referenceZ,ZH2O)
            

            
            XH2OModulus=numpy.sqrt((XH2O*XH2O).sum())
            crossLineModulus=numpy.sqrt((crossLine*crossLine).sum())
##            print crossLineModulus
            dot1=numpy.dot(crossLine,XH2O)
            
            cosAnglePsi=dot1/XH2OModulus/crossLineModulus
##            print cosAnglePsi
            anglePsi=numpy.arccos(cosAnglePsi)
##            print anglePsi
            if XH2O[2]<0:
                anglePsi=2*pi-anglePsi
            for number in range(10):
                nextNumber=number+1
                if nextNumber*pi/5.0>anglePsi>=number*pi/5.0:
                    eulerPsiNo[number]+=1
                elif anglePsi==2*pi:
                    eulerPsiNo[9]+=1
#        print 'the euler Psi No is:',CentreNo, eulerPsiNo
        for number in range(len(eulerPsiNo)):
            eulerPsiNo[number]=eulerPsiNo[number]/eulerPsiDensity
        eulerPsi.append(eulerPsiNo)
#        print 'the euler Psi No is:',CentreNo, eulerPsiNo
    return eulerPsi

##getH2OEulerPsi is used to calculate the Euler distrubution of all the waters.
##The water coordinates were saved as lists with both Oxygen, and Hydeogen atoms.
##The coordinates looks like[[[OxygenCoord],[Hydgrogen1 coordinate],[Hydrogen2 Coordinate]....]...]
##the output of the Euler angle is very similar with other functions. the distribution is
##divided into 20 intergration spaces. The output is [[g(0~*pi/20.0),g(pi/20.0~2*pi/20.0)...g(pi~19*pi/20.0)]......]

def getH2OCoordForWW(self):
    H2OCoordForWW=[]
    for H2OInforCenAndAround in self:
        H2OCoordCentre=getH2OCoord(H2OInforCenAndAround)
        H2OCoordForWW.append(H2OCoordCentre)
    return H2OCoordForWW

##self is all the H2OinforForWW. therefore, I use the function getH2OCoord to extract all the H2O coordidates. the output of the
##function looks like the output of 'H2OInforNearbyHighCentre' the output:[[[highOccupiedH2OCoord],[aroundH2OCoord2],[around
##H2OCoord2]],[...]...]. Moreover, the output centre is consistant with centre of findNearByCentre output.

def calculateEulerThetaReachCentre(self,nearbyCentre,frameNo):
    eulerThetaEacCentreAndAround=[]
    for eachH2ONo in range(len(self)):
        eachCentre=nearbyCentre[eachH2ONo]
        eachH2O=self[eachH2ONo]
        eacheulerTheta=getH2OEulerTheta(eachH2O,frameNo)
        eulerThetaEacCentreAndAround.append(eacheulerTheta)
    return  eulerThetaEacCentreAndAround
##calculateEulerThetaReachCentre is used to calculate the gEulerTheta of for all the centre H2O and the arounded H2Os.
##The input included self: H2O coordinates which is saved as list of H2O and arounded H2O. This can be the
##output of function 'getH2OCoordForWW'. nearbycentre is all the centres for the calculation. The centre
##information should be consistant with the self(H2OCoordinate). Morover, it can read the output of function
##'findNearbyCentre'. The output of this function looks like[[[centreGEularTheta],[aroundH2O GEularTheta]...]......]

def calculateEulerPhiReachCentre(self,nearbyCentre,frameNo):
    eulerPhiEacCentreAndAround=[]
    for eachH2ONo in range(len(self)):
        eachCentre=nearbyCentre[eachH2ONo]
        eachH2O=self[eachH2ONo]
        eacheulerPhi=getH2OEulerPhi(eachH2O,frameNo)
        eulerPhiEacCentreAndAround.append(eacheulerPhi)
    return  eulerPhiEacCentreAndAround
##calculateEulerPhiReachCentre is used to calculate the gEulerPhi of for all the centre H2O and the arounded H2Os.
##The input included self: H2O coordinates which is saved as list of H2O and arounded H2O. This can be the
##output of function 'getH2OCoordForWW'. nearbycentre is all the centres for the calculation. The centre
##information should be consistant with the self(H2OCoordinate). Morover, it can read the output of function
##'findNearbyCentre'. The output of this function looks like[[[centreGEularPhi],[aroundH2O GEularPhi]...]......]

def calculateEulerPsiReachCentre(self,nearbyCentre,frameNo):
    eulerPsiEacCentreAndAround=[]
    for eachH2ONo in range(len(self)):
        eachCentre=nearbyCentre[eachH2ONo]
        eachH2O=self[eachH2ONo]
        eacheulerPsi=getH2OEulerPsi(eachH2O,frameNo)
        eulerPsiEacCentreAndAround.append(eacheulerPsi)
    return  eulerPsiEacCentreAndAround
##calculateEulerPsiReachCentre is used to calculate the gEulerTheta of for all the centre H2O and the arounded H2Os.
##The input included self: H2O coordinates which is saved as list of H2O and arounded H2O. This can be the
##output of function 'getH2OCoordForWW'. nearbycentre is all the centres for the calculation. The centre
##information should be consistant with the self(H2OCoordinate). Morover, it can read the output of function
##'findNearbyCentre'. The output of this function looks like[[[centreGEularPsi],[aroundH2O GEularPsi]...]......]




def intgGr(self):
    sumGr=0
    for eachGrNo in range(len(self)):
        eachGr=self[eachGrNo]
        nextNo=eachGrNo+1
        eachGr=eachGr*(((nextNo*0.05)**3)-((eachGrNo*0.05)**3))/3.0
        sumGr+=eachGr
##    print sumGr
    return sumGr

##intgGr is used to do the integration of GR according to the integral of 0.05. The input is a list a GR
##along the distance R. The output of this function is sum(gr*dr)

def intgGtheta(self):
    sumGtheta=0
    for gthetaNo in range(len(self)):
        eachGtheta=self[gthetaNo]
        angle1=gthetaNo*pi/20.0
        angle2=(gthetaNo+1)*pi/20.0
        eachGtheta=eachGtheta*(cos(angle1)-cos(angle2))
        sumGtheta+=eachGtheta
    return sumGtheta

##intgGtheta is used to do the integration of Gtheta according to the integral of pi/20.0. The input is a list a Gtheta
##along each angle. The output of this function is sum(gtheta*sinTheta*dTheta).

def intgGangle(self):
    sumGangle=0
    for eachGangle in self:
        eachGangle=eachGangle*pi/20.0
        sumGangle+=eachGangle
    return sumGangle

##intgGangle is used to do the integration of Gangle according to the integral of pi/20.0. The input is a list a Gangle
##along each angle. The output of this function is sum(gangle*dangle)

def getWforWW(self,nearbyCentre,waterCoordAroundCentre,frameNo):
    WforWW=[] #single terms for WW entropy calculation
##    for centreSetNo in range(len(self)):
##        grOOSet=self[centreSetNo]
##        OcoordSet=waterCoordAroundCentre[centreSetNo]
##        H2OCoordSet=H2OCoordForWW[centreSetNo]
    gr=calculateGReachCentre(waterCoordAroundCentre,nearbyCentre,frameNo)
    grtheta=calculateGThetaeachCentre(waterCoordAroundCentre,nearbyCentre,frameNo)
    grphi=calculateGPhieachCentre(waterCoordAroundCentre,nearbyCentre,frameNo)
##    gEtheta=calculateEulerThetaReachCentre(H2OCoordForWW,nearbyCentre,frameNo)
##    gEphi=calculateEulerPhiReachCentre(H2OCoordForWW,nearbyCentre,frameNo)
##    gEpsi=calculateEulerPsiReachCentre(H2OCoordForWW,nearbyCentre,frameNo)
    for centreSetNo in range(len(self)):
        WforWWEachCentreAndAround=[]
        
        gRconstantSet=self[centreSetNo]
        grSet=gr[centreSetNo]
        grthetaSet=grtheta[centreSetNo]
        grphiSet=grphi[centreSetNo]
##        gEthetaSet=gEtheta[centreSetNo]
##        gEphiSet=gEphi[centreSetNo]
##        gEpsiSet=gEpsi[centreSetNo]
        
        cGr=grSet[0] ##cGr is centreGr
        cGrtheta=grthetaSet[0]
        cGrphi=grphiSet[0]
##        cGEtheta=gEthetaSet[0]
##        cGEphi=gEphiSet[0]
##        cGEpsi=gEpsiSet[0]
        
        cIntgGr=intgGr(cGr)
        cIntgGrtheta=intgGtheta(cGrtheta)
        cIntgGrphi=intgGangle(cGrphi)
        
##        cIntgGEtheta=intgGtheta(cGEtheta)
##        cIntgGEphi=intgGangle(cGEphi)
##        cIntgGEpsi=intgGangle(cGEpsi)

        for eachpairNo in range(len(gRconstantSet)):
            gRconstant=gRconstantSet[eachpairNo]
            
            aGr=grSet[eachpairNo] ##aGr is aroundGr
            aGrtheta=grthetaSet[eachpairNo]
            aGrphi=grphiSet[eachpairNo]
##            aGEtheta=gEthetaSet[eachpairNo]
##            aGEphi=gEphiSet[eachpairNo]
##            aGEpsi=gEpsiSet[eachpairNo]
        
            aIntgGr=intgGr(aGr)
            aIntgGrtheta=intgGtheta(aGrtheta)
            aIntgGrphi=intgGangle(aGrphi)
                
##            aIntgGEtheta=intgGtheta(aGEtheta)
##            aIntgGEphi=intgGangle(aGEphi)
##            aIntgGEpsi=intgGangle(aGEpsi)
#            print 'gRconstant*cIntgGr*cIntgGrtheta*cIntgGrphi*aIntgGr*aIntgGrtheta*aIntgGrphi*cIntgGEtheta*cIntgGEphi*cIntgGEpsi*aIntgGEtheta*aIntgGEphi*aIntgGEpsi' 
#            print gRconstant,cIntgGr,cIntgGrtheta,cIntgGrphi,aIntgGr,aIntgGrtheta,aIntgGrphi#,cIntgGEtheta,cIntgGEphi,cIntgGEpsi,aIntgGEtheta,aIntgGEphi,aIntgGEpsi   
            eachWforWW=(-0.5)*k*(pw**2)*mol*gRconstant*cIntgGr*cIntgGrtheta*cIntgGrphi*aIntgGr*aIntgGrtheta*aIntgGrphi#*cIntgGEtheta*cIntgGEphi*cIntgGEpsi*aIntgGEtheta*aIntgGEphi*aIntgGEpsi
	    #eachWforWW=(-0.5)*k*(pw**2)*((1.0/(8.0*pi**2))**2)*mol*gRconstant*cIntgGr*cIntgGrtheta*cIntgGrphi*aIntgGr*aIntgGrtheta*aIntgGrphi#*cIntgGEtheta*cIntgGEphi*cIntgGEpsi*aIntgGEtheta*aIntgGEphi*aIntgGEpsi
            WforWWEachCentreAndAround.append(eachWforWW)
        WforWW.append(WforWWEachCentreAndAround)
    return WforWW

##getWforWW is used to calculate the S(WW). As the S(WW) can be divided into one body terms and two body terms,this function can obtain all the one body terms
##for the S(WW) calculatoin. A lot of information is needed for this calculation. self: all the G(R) for the Oxygen-Oxygen interaction in bulk water.
##nearbyCentre: all the centres and arounded waters for calculation, it can be the output of function 'findNearbyCentre'.waterCoordAroundCentre: is all the water
##ocygen cooordidates for the calculation. it can be the output of function 'extractWaterCoorAround'. H2OCoordForWW: all the H2O coordinate including
##H1,H2, OW coordinates for the calculation. it can be the output of the function 'H2OCoordForWW'. frameNo is the total frameNo for the simulation. Pls be mention
##all the input were saved as list loops which is in accordence with each other. The output of the function is WforWW, which has the same format of self. it looks
##like [[0,centre1Around1SingleTerm,centre1Around2SingleTerm...],[]...]. Keep in mind there is an zero in WforWW, that is because the centre-centre interaction is
##zero.


def getEularAngle(H2OCoord):
    referenceH1=numpy.array([0.79079641377315202, 0.61207926934631729, 0.0])
    referenceH2=numpy.array([-0.79079641377315202, 0.61207926934631729, 0.0])
    referenceZ=numpy.array([0.0,0.0,1.0])
    referenceY=numpy.array([0.0,1.0,0.0])
    referenceX=numpy.array([1.0,0.0,0.0])
    revReferenceY=numpy.array([0.0,-1.0,0.0])
    euler=[]


    H1=numpy.array(H2OCoord[1])

    H2=numpy.array(H2OCoord[2])

    OT=numpy.array(H2OCoord[0])
            
    H1=H1-OT
    H2=H2-OT
    
    YH2O=numpy.add(H1,H2)
    ZH2O=numpy.cross(H1,H2)
    XH2O=numpy.cross(YH2O,ZH2O)
            
            
    ZH2OOnXY=numpy.copy(ZH2O)

    ZH2OOnXY[2]=0.0
    crossLine=numpy.cross(referenceZ,ZH2O)
            

    ZH2OModulus=numpy.sqrt((ZH2O*ZH2O).sum())       
    XH2OModulus=numpy.sqrt((XH2O*XH2O).sum())
    ZH2OOnXYModulus=numpy.sqrt((ZH2OOnXY*ZH2OOnXY).sum())
    crossLineModulus=numpy.sqrt((crossLine*crossLine).sum())
    
##            print crossLineModulus
    dotTheta=numpy.dot(referenceZ,ZH2O)
    dotPhi=numpy.dot(revReferenceY,ZH2OOnXY)
    dotPsi=numpy.dot(crossLine,XH2O)

    cosAngleTheta=dotTheta/ZH2OModulus
    cosAnglePhi=dot1/ZH2OOnXYModulus
    cosAnglePsi=dotPsi/XH2OModulus/crossLineModulus
##            print cosAnglePsi
    angleTheta=numpy.arccos(cosAngleTheta)
    anglePhi=numpy.arccos(cosAnglePhi)
    anglePsi=numpy.arccos(cosAnglePsi)
##            print anglePsi
    if ZH2OOnXY[0]<0:
        anglePhi=2*pi-anglePhi
    if XH2O[2]<0:
        anglePsi=2*pi-anglePsi
    euler.append(angleTheta)
    euler.append(anglePhi)
    euler.append(anglePsi)
    return euler


##getH2OEulerAngle is used to calculate the Euler angle for one water.
##The water coordinates were saved as lists with both Oxygen, and Hydeogen atoms.
##The coordinates looks like[OxygenCoord],[Hydgrogen1 coordinate],[Hydrogen2 Coordinate]
##the output of the Euler angle [Theta, Phi, Psi]

def getThetaCosDis():
    pi5=pi/5
    pi10=pi/10
    pi20=pi/20
    noList=range(10)
    cosAngle=[]
    for No in noList:
        angle=No*pi10+pi20
        cosAngle.append(cos(angle))
    return cosAngle

def getThetaSinDis():
    pi5=pi/5
    pi10=pi/10
    pi20=pi/20
    noList=range(10)
    sinAngle=[]
    for No in noList:
        angle=No*pi10+pi20
        sinAngle.append(sin(angle))
    return sinAngle

def getAngleCosDis():
    pi5=pi/5
    pi10=pi/10
    pi20=pi/20
    noList=range(10)
    cosAngle=[]
    for No in noList:
        angle=No*pi5+pi10
        cosAngle.append(cos(angle))
    return cosAngle

def getAngleSinDis():
    pi5=pi/5
    pi10=pi/10
    pi20=pi/20
    noList=range(10)
    sinAngle=[]
    for No in noList:
        angle=No*pi5+pi10
        sinAngle.append(sin(angle))
    return sinAngle

## getThetaCosDis,getThetaSinDis,getAngleCosDis,getAngleSinDis is used to calculate the cos angle and sin sinangle dictionary for the
## eykarRotation calculation.
        

def eularRotation(self,thetaCos,thetaSin,angleCos,angleSin):
    referenceH1=numpy.array([0.79079641377315202, 0.61207926934631729, 0.0])
    referenceH2=numpy.array([-0.79079641377315202, 0.61207926934631729, 0.0])
    referenceZ=numpy.array([0.0,0.0,1.0])
    referenceY=numpy.array([0.0,1.0,0.0])
    revReferenceY=numpy.array([0.0,-1.0,0.0])
##    angleTheta=self[0]
##    anglePhi=self[1]
##    anglePsi=self[2]
    angleTheta=self[0]
    anglePhi=self[1]
    anglePsi=self[2]
##    a11=cos(anglePsi)*cos(anglePhi)-cos(angleTheta)*sin(anglePhi)*sin(anglePsi)
##    a12=cos(anglePsi)*sin(anglePhi)+cos(angleTheta)*cos(anglePhi)*sin(anglePsi)
##    a13=sin(anglePsi)*sin(angleTheta)
##    a21=-sin(anglePsi)*cos(anglePhi)-cos(angleTheta)*sin(anglePhi)*cos(anglePsi)
##    a22=-sin(anglePsi)*sin(anglePhi)+cos(angleTheta)*cos(anglePhi)*cos(anglePsi)
##    a23=cos(anglePsi)*sin(angleTheta)
##    a31=sin(angleTheta)*sin(anglePhi)
##    a32=-sin(angleTheta)*cos(anglePhi)
##    a33=cos(angleTheta)
    a11=angleCos[anglePsi]*angleCos[anglePhi]-thetaCos[angleTheta]*angleSin[anglePhi]*angleSin[anglePsi]
    a12=angleCos[anglePsi]*angleSin[anglePhi]+thetaCos[angleTheta]*angleCos[anglePhi]*angleSin[anglePsi]
    a13=angleSin[anglePsi]*thetaSin[angleTheta]
    a21=-angleSin[anglePsi]*angleCos[anglePhi]-thetaCos[angleTheta]*angleSin[anglePhi]*angleCos[anglePsi]
    a22=-angleSin[anglePsi]*angleSin[anglePhi]+thetaCos[angleTheta]*angleCos[anglePhi]*angleCos[anglePsi]
    a23=angleCos[anglePsi]*thetaSin[angleTheta]
    a31=thetaSin[angleTheta]*angleSin[anglePhi]
    a32=-thetaSin[angleTheta]*angleCos[anglePhi]
    a33=thetaCos[angleTheta]
    
    rotation=numpy.matrix([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])
    
    H1=((rotation**-1)*(referenceH1.reshape(-1,1))).reshape(1,-1)
    H2=((rotation**-1)*(referenceH2.reshape(-1,1))).reshape(1,-1)
    H1=numpy.array(H1)
    H2=numpy.array(H2)
    O=numpy.array([0,0,0])
    H2O=numpy.array([O,H1[0],H2[0]])
    return H2O

##eularRotation is a function which can be used to the H2O coordinate according the the eular angle
##Theta, phi, Psi. the output is the rotated numpy.array [O,H1,H2] coordinates.


def rotationAngle(self,vector2,angle):
    epsilon=10**(-10)
    v1xy=[]
    v2xy=[]
    v1yz=[]
    v2yz=[]
    
    v1xy.extend(self[:2])
    v2xy.extend(vector2[:2])
    v1yz.extend(self[1:])
    v2yz.extend(vector2[1:])
    
    v1xy=numpy.array(v1xy)
    v2xy=numpy.array(v2xy)
    v1yz=numpy.array(v1yz)
    v2yz=numpy.array(v2yz)
    
    
    v1xyM=numpy.sqrt((v1xy*v1xy).sum())
    v2xyM=numpy.sqrt((v2xy*v2xy).sum())
    v1yzM=numpy.sqrt((v1yz*v1yz).sum())
    v2yzM=numpy.sqrt((v2yz*v2yz).sum())
    dotxy=numpy.dot(v1xy,v2xy)
    dotyz=numpy.dot(v1yz,v2yz)
    planeAngle=numpy.arccos(dotxy/v1xyM/v2xyM)

    if abs(numpy.cross(v1xy,v2xy))<epsilon:
        
        if abs(numpy.cross(v1yz,v2yz))<epsilon:
            anlge=0.0
        elif abs(numpy.cross(v1yz,v2yz))>=epsilon:
            yzCross=numpy.cross(v1yz,v2yz)
            if yzCross<=0:
                angle=2*pi-angle
##                print angle
            
    elif abs(numpy.cross(v1xy,v2xy))>=epsilon:
        xyCross=numpy.cross(v1xy,v2xy)
        if xyCross<=0:
            angle=2*pi-angle
##            print angle
    return angle

##rotationAngle is used to determine the rotation angle between two vectors. As we know that
##the angle obtained from the vector arccos is the angle between the two vectors,
##therefore it is just the value from 0-pi, if we want to determine the rotation angle
##we have to perform the other method. In this function, I use reflection the 3D vector
##to 2D vector to determine whether the angle is clockwise or anticlockwise. If the angle
##is anticlockwise, the angle becomes to 2*pi-angle. the input of the function included
##self: the first 3D vector for use. vector2 is the second vector. anlge is the
##3D angle between self and vector2. The output is the rotation angle between the
##two vectors.

def arcCosAngle(cosPhi):

    if cosPhi>=1:
        phi=0
    elif -1<cosPhi<1:
        phi=numpy.arccos(cosPhi)
    elif cosPhi<=-1:
        phi=pi-0.000001
    return phi


def getCorrelationAnlge(centreH2O,aroundH2O):
    pi10=pi/10
    pi20=pi/20
    correlationAngle=[]
    
    centreOCoord=centreH2O[0]
    centreH1Coord=centreH2O[1]
    centreH2Coord=centreH2O[2]
    
    aroundOCoord=aroundH2O[0]
    aroundH1Coord=aroundH2O[1]
    aroundH2Coord=aroundH2O[2]
    
    O1D1=(numpy.add(centreH1Coord,centreH2Coord))/2-centreOCoord
    O1H11=centreH1Coord-centreOCoord
    O1H12=centreH2Coord-centreOCoord
    O2D2=(numpy.add(aroundH1Coord,aroundH2Coord))/2-aroundOCoord
    O2H21=aroundH1Coord-aroundOCoord
    O2H22=aroundH2Coord-aroundOCoord
    O1O2=aroundOCoord-centreOCoord
    O2O1=centreOCoord-aroundOCoord
    H11H12=centreH2Coord-centreH1Coord
    H21H22=aroundH2Coord-aroundH1Coord
    crossO1O2O1D1=numpy.cross(O1O2,O1D1)
    crossO2D2O2O1=numpy.cross(O2D2,O2O1)
##  print crossO2D2O2O1
                            

    O1D1Modulus=numpy.sqrt((O1D1*O1D1).sum())
    O2D2Modulus=numpy.sqrt((O2D2*O2D2).sum())
    O1O2Modulus=numpy.sqrt((O1O2*O1O2).sum())
    O2O1Modulus=numpy.sqrt((O2O1*O2O1).sum())
    H11H12Modulus=numpy.sqrt((H11H12*H11H12).sum())
    H21H22Modulus=numpy.sqrt((H21H22*H21H22).sum())
    crossO1O2O1D1Modulus=numpy.sqrt((crossO1O2O1D1*crossO1O2O1D1).sum())
    crossO2D2O2O1Modulus=numpy.sqrt((crossO2D2O2O1*crossO2D2O2O1).sum())
                            


    dotTheta1=numpy.dot(O1D1,O1O2)
    dotTheta2=numpy.dot(O2D2,O2O1)
    dotPhi=numpy.dot(crossO1O2O1D1,crossO2D2O2O1)
##                            print dotPhi
    dotChi1=numpy.dot(H11H12,crossO1O2O1D1)
    dotChi2=numpy.dot(H21H22,crossO2D2O2O1)


    cosTheta1=dotTheta1/(O1D1Modulus*O1O2Modulus)
    cosTheta2=dotTheta2/(O2D2Modulus*O2O1Modulus)
    
    cosPhi=dotPhi/(crossO1O2O1D1Modulus*crossO2D2O2O1Modulus)
    cosChi1=dotChi1/(H11H12Modulus*crossO1O2O1D1Modulus)
    cosChi2=dotChi2/(H21H22Modulus*crossO2D2O2O1Modulus)

##    if cosTheta1>=1:
##        theta1=0
##    elif cosTheta1>=1:
##
##        theta1=numpy.arccos(cosPhi)
    theta1=numpy.arccos(cosTheta1)
    theta2=numpy.arccos(cosTheta2)

    theta1=arcCosAngle(cosTheta1)
    theta2=arcCosAngle(cosTheta2)
    phi=arcCosAngle(cosPhi)
    chi1=arcCosAngle(cosChi1)
    chi2=arcCosAngle(cosChi2)
##
##    if cosPhi>=1:
##        phi=0
##    elif -1<cosPhi<1:
##        phi=numpy.arccos(cosPhi)
##    elif cosPhi<=-1:
##        phi=pi-0.000001
    
    
##    print crossO2D2O2O1Modulus,crossO1O2O1D1Modulus,cosPhi,phi

##    phi=rotationAngle(crossO1O2O1D1,crossO2D2O2O1,phi)
##
####    chi1=numpy.arccos(cosChi1)
##    chi1=rotationAngle(H11H12,crossO1O2O1D1,chi1)
##                            
####    chi2=numpy.arccos(cosChi2)
##    chi2=rotationAngle(H21H22,crossO2D2O2O1,chi2)
##    if theta1>pi or theta1<0:
##        print 'theta1',theta1
##    if theta2>pi or theta2<0:
##        print 'theta2',theta2
##    if phi>=2*pi:
##        phi=2*pi-10*(10**(-5))
##    elif phi<=0:
##        phi=10*(10**(-5))
##    if chi1>2*pi or chi1<0:
##        print 'chi1',chi1
##    if chi2>2*pi or chi2<0:
##        print 'chi2',chi2
    

    theta1No=int(theta1/pi10)
    theta2No=int(theta2/pi10)
    phiNo=int(phi/pi20)
    chi1No=int(chi1/pi20)
    chi2No=int(chi2/pi20)
    if chi1No==20:
        chi1No=random.sample([19,0],1)
    if chi2No==20:
        chi2No=random.sample([19,0],1)
    if theta1No==10:
        theta1No=random.sample([9,0],1)
    if theta2No==10:
        theta2No=random.sample([9,0],1)
    if phiNo==20:
        phiNo=random.sample([19,0],1)
    correlationAngle.append(theta1No)
    correlationAngle.append(theta2No)
    correlationAngle.append(phiNo)
    correlationAngle.append(chi1No)
    correlationAngle.append(chi2No)
    return correlationAngle

##correlationAngle is used to correlation angle between two H2O molecules. As
##we know the eular anlge between each water. From the Eular angle, we got the
##coordinates from two H2O coordidated from the function 'eularRotation', Then,
##the Correlation anlgle can be obtained from these two H2O coordinates. The
##output of the is a aeries of angles saved as list. With these angles, we can
##read the bulk water correllation g(w,w) from file.





def strCoordToFlost(self):
    coordinate=[]
    for eachCoord in self:
        eachCoord=float(eachCoord)
        coordinate.append(eachCoord)
    return coordinate

##strCoordToFlost is used to translate all the string list coordinates into the folat format numerical
##list.

def openBulkCoorFile(self):
    correlation=[]
    dicCorrelation=[]
    correlationFile=open(self,'r')
    fiveAngle=[]
    
    t1=[0,0,0,0,0,0,0,0,0,0]
    for t1No in range(10):
        t1[t1No]=[0,0,0,0,0,0,0,0,0,0]
    for t1No in range(10):
        for t2No in range(10):
            t1[t1No][t2No]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for t1No in range(10):
        for t2No in range(10):
            for phiNo in range(20):
                t1[t1No][t2No][phiNo]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for t1No in range(10):
        for t2No in range(10):
            for phiNo in range(20):
                for c1No in range(20):
                    t1[t1No][t2No][phiNo][c1No]=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
           
                  
        
        
    for line in correlationFile:
        line=line.split()
        line=strCoordToFlost(line)
        correlation.append(line)

        t1No=int(line[1])
        t2No=int(line[2])
        phiNo=int(line[3])
        c1No=int(line[4])
        c2No=int(line[5])
        value=line[6]
##        print t1No,t2No,phiNo,c1No,c2No
        t1[t1No][t2No][phiNo][c1No][c2No]=value
        
##    for line in correlation:
##        if line[6]>=0:
##            optimizedCorrelation.append(line)
            
    
    return t1

##openBulkCoorFile is used to open the bulk water Coorelation file and save them as a list.

def getCentreAroundDistance(self):
    distance=[]
    
    for setNo in range(len(self)):
        eachDistance=[0]
##        grInhTransSet=[0]
        centreSet=self[setNo]
##        waterCoordSet=waterCoordAround[setNo]
##        print waterCoordSet
##        print 'there were :',len(waterCoordSet),'sets in waterCoordSet'
        totalNo=len(centreSet)
        referenceCentre=centreSet[0]
        referenceCentreArray=numpy.array(referenceCentre)
        for aroundNo in range(totalNo):
            if aroundNo>=1:
                aroundCentre=centreSet[aroundNo]
                aroundCentreArray=numpy.array(aroundCentre)
                distanceVector=referenceCentreArray-aroundCentreArray
                distanceModulus=numpy.sqrt((distanceVector*distanceVector).sum())
                eachDistance.append(distanceModulus)
        distance.append(eachDistance)
    return distance
##getCentreAroundDistance is used to calculate the distance between the centre and the around waters. The input is self:
##the water oxygen coordinated and the arounded waters saved as list according to each centre. This can be the output of
##'findNearbyCentre'. The output is the the distance according to between centre water and arounded water. The format of
##is also the list saved with similar format with input.


                


def getSmWW(self,waterCoordAroundCentre,nearbyCentre,H2OCoordForWW,frameNo):
#self is the distance of the centre to the arounded waters.
##    print 'there were ',len(self),'waters for calculation'
    SmWW=[]
    gEtheta=calculateEulerThetaReachCentre(H2OCoordForWW,nearbyCentre,frameNo)
    gEphi=calculateEulerPhiReachCentre(H2OCoordForWW,nearbyCentre,frameNo)
    gEpsi=calculateEulerPsiReachCentre(H2OCoordForWW,nearbyCentre,frameNo)
    
    pi5=pi/5
    pi10=pi/10
    pi20=pi/20
    spacePower=pi5**4

    thetaCos=getThetaCosDis()
    thetaSin=getThetaSinDis()
    angleCos=getAngleCosDis()
    angleSin=getAngleSinDis()
    
    
    cosSpaceValue=[]
    for i in range(10):
        spaceValue=cos(i*pi10)-cos((i+1)*pi10)
        cosSpaceValue.append(spaceValue)

##    cEAngle=[]
##    aEAngle=[]
    
    for eachCentreNo in range(len(self)):
        SmWWEachCentreAndAround=[]

        eachNearbyCentre=nearbyCentre[eachCentreNo]
        eachWaterNearbyCentre=waterCoordAroundCentre[eachCentreNo]
        eachCentreDistance=self[eachCentreNo]

        gEthetaSet=gEtheta[eachCentreNo]
        gEphiSet=gEphi[eachCentreNo]
        gEpsiSet=gEpsi[eachCentreNo]

        cWaters=eachWaterNearbyCentre[0]
        cCoord=eachNearbyCentre[0]
        cGEtheta=gEthetaSet[0]
        cGEphi=gEphiSet[0]
        cGEpsi=gEpsiSet[0]
        if len(cWaters)<=(frameNo*0.4):
            
            eachSmWW=[0]*len(eachNearbyCentre)
            SmWWEachCentreAndAround=eachSmWW
##            SmWW.append(SmWWEachCentreAndAround)
            print eachSmWW
            print eachNearbyCentre
        elif len(cWaters)>=(frameNo*0.4):
            
            
            for eachPairNo in range(len(eachCentreDistance)):
                pairDistance=eachCentreDistance[eachPairNo]
                if eachPairNo==0:
                    eachSmWW=0
                
                elif eachPairNo>=1:
                    eachSmWW=0
                    aCoord=eachNearbyCentre[eachPairNo]
                    move=numpy.array(aCoord)-numpy.array(cCoord)
                
                    for space in range(20,40):
                        space=space/10.0
                        if abs(pairDistance-space)<=0.05:
                            pairDistance=space
                    fileName='../eular2CorrSymm/distance_'+str(int((pairDistance-2)/0.1))+'.txt'
                    
                    distanceFile=openBulkCoorFile(fileName)
                    print 'the distance File includes',len(distanceFile)

                    aGEtheta=gEthetaSet[eachPairNo]
                    aGEphi=gEphiSet[eachPairNo]
                    aGEpsi=gEpsiSet[eachPairNo]
                    eachSmWW=0

                    for cGEpsiNo in range(len(cGEpsi)):
                        cGEEachPsi=cGEpsi[cGEpsiNo]
                        
                        if cGEEachPsi==0:
                            eachAngleSmWW=0
                        elif cGEEachPsi!=0:
##                            cRotaPsi=cGEpsiNo*pi5+pi10
                            for cGEphiNo in range(len(cGEphi)):
                                cGEEachPhi=cGEphi[cGEphiNo]
                                        
                                if cGEEachPhi==0:
                                    eachAngleSmWW=0
                                elif cGEEachPhi!=0:
##                                    cRotaPhi=cGEphiNo*pi5+pi10
                                    
                                    for cGEthetaNo in range(len(cGEtheta)):
                                        cGEEachTheta=cGEtheta[cGEthetaNo]
                                        if cGEEachTheta==0:
                                             eachAngleSmWW=0
                                        elif cGEEachTheta!=0:
##                                            cRotaTheta=cGEthetaNo*pi10+pi20
                                            
                                            cEangleNo=[]
                                            cEangleNo.append(cGEthetaNo)
                                            cEangleNo.append(cGEphiNo)
                                            cEangleNo.append(cGEpsiNo)
                                            cH2OR=eularRotation(cEangleNo,thetaCos,thetaSin,angleCos,angleSin)
                            
                                            for aGEpsiNo in range(len(aGEpsi)):
                                                aGEEachPsi=aGEpsi[aGEpsiNo]
                                                if aGEEachPsi==0:
                                                    eachAngleSmWW=0
                                
                                                elif aGEEachPsi!=0:
##                                                    aRotaPsi=aGEpsiNo*pi5+pi10
                                    
                                                    for aGEphiNo in range(len(aGEphi)):
                                                        aGEEachPhi=aGEphi[aGEphiNo]
                                                        if aGEEachPhi==0:
                                                           eachAngleSmWW=0 
                                                        elif aGEEachPhi!=0:
##                                                             aRotaPhi=aGEphiNo*pi5+pi10
                                                             
                                                             for aGEthetaNo in range(len(aGEtheta)):
                                                                aGEEachTheta=aGEtheta[aGEthetaNo]
                                                                if aGEEachTheta==0:
                                                                    eachAngleSmWW=0
                                                                elif aGEEachTheta!=0:
##                                                                    aRotaTheta=aGEthetaNo*pi10+pi20
                                                                                                                                    
                                                                    aEangleNo=[]

                                                                    aEangleNo.append(aGEthetaNo)
                                                                    aEangleNo.append(aGEphiNo)
                                                                    aEangleNo.append(aGEpsiNo)

                                                                    

                                                                    aH2OR=eularRotation(aEangleNo,thetaCos,thetaSin,angleCos,angleSin)
                                                                    aH2OR=aH2OR+move

                                                                    correlationAngle=getCorrelationAnlge(cH2OR,aH2OR)
                                                                    t1No=correlationAngle[0]
                                                                    t2No=correlationAngle[1]
                                                                    phiNo=correlationAngle[2]
                                                                    c1No=correlationAngle[3]
                                                                    c2No=correlationAngle[4]
                                                                    gbulk=distanceFile[t1No][t2No][phiNo][c1No][c2No]
                                                                    if gbulk!=0:
                                                                        
                                                                        eachAngleSmWW=cGEEachTheta*cGEEachPhi*cGEEachPsi*aGEEachTheta*aGEEachPhi*aGEEachPsi*(gbulk/(2*pi))*(numpy.log(gbulk/(2*pi)))
                                                                        eachAngleSmWW=eachAngleSmWW*(cosSpaceValue[cGEthetaNo])*(cosSpaceValue[aGEthetaNo])*spacePower
                                                                        eachSmWW+=eachAngleSmWW
##                                                                        if gbulk>=100:
                                                                            
##                                                                        print 'centreNo is:',eachCentreNo,cGEEachTheta,cGEEachPhi,cGEEachPsi,aGEEachTheta,aGEEachPhi,aGEEachPsi,gbulk,eachAngleSmWW
##                                                                        print int((pairDistance-2)/0.1),t1No,t2No,phiNo,c1No,c2No,gbulk,eachSmWW
                                                                    
##                                                                    for line in distanceFile:
##                                                                        if line[1:6]==correlationAngle:
##                                                                            eachAngleSmWW=cGEEachTheta*cGEEachPhi*cGEEachPsi*aGEEachTheta*aGEEachPhi*aGEEachPsi*(line[6])
##                                                                            eachAngleSmWW=eachAngleSmWW*(cosSpaceValue[cGEthetaNo])*(cosSpaceValue[aGEthetaNo])*spacePower
##                                                              
##                                                                            eachSmWW+=eachAngleSmWW
##                                                                            print eachSmWW
##                                                                            break
                  
                eachSmWW=eachSmWW/((8*pi**2)**2)
                print eachCentreNo,'with the correlation', eachSmWW
                print eachSmWW
                print eachNearbyCentre
                SmWWEachCentreAndAround.append(eachSmWW)
                print 'SmWWEachCentreAndAround is',SmWWEachCentreAndAround
        SmWW.append(SmWWEachCentreAndAround)
        print 'SmWW is',SmWW
    return SmWW



##calculateSorWW(self,WforWW) is used to calculate the final result of the Water-water orientation
##entropy. the input include self: the five anlges representing the Water-water angle, it is the output
##of function 'getWWAngle', the other input is the WforWW: which is the output of function getWforWW.
##the output is entropy of W-W orientation each centre saved as list.

def getSorWW(WWCentreNo,allSmWW,WforWW):
    SorWW=[]
##    line='##########################Final Sor Correlation####################'
##    sorOutput.write(line)
    for No in range(len(WWCentreNo)):
        eachCentreSorWW=[]	
        centreNo=WWCentreNo[No]
        eachSmWW=allSmWW[No]
        eachWforWW=WforWW[No]
        
        for eachPairNo in range(len(eachSmWW)):
            eachPairSmWW=eachSmWW[eachPairNo]
            eachPairWforWW=eachWforWW[eachPairNo]
            eachPairSorWW=eachPairSmWW*eachPairWforWW
##            line=str(centreNo)+'    '+str(eachPairSorWW)+'\n'
##            sorOutput.write(line)
            eachPairSorWW=eachPairSorWW/4.184
            eachCentreSorWW.append(eachPairSorWW)
            
##        eachCentreSorWW=eachCentreSorWW/4.184
##        line=str(centreNo)+'    '+str(eachCentreSorWW)+'    '+'Final'+'\n'
##        sorOutput.write(line)
        
##    sorOutput.close()
        SorWW.append(eachCentreSorWW)
    print WWCentreNo
    print SorWW
    return SorWW
##getSorWW is used to calculate the SorWW final result for the calculation.

def printOutWWOrien(self,centreAndAroundNo,WWOentropy):
    line='#CentreNo'+'  '+'WWOrienEntropy'+'\n'
    self.write(line)
    
    for No in range(len(centreAndAroundNo)):
        centreSetNo=centreAndAroundNo[No]

        entropySet=WWOentropy[No]
        centreNo=centreSetNo[0]
        entropy=sum(entropySet)
        line=str(centreNo)+'    '+str(entropy)+'\n'
        self.write(line)
    line='################explaination of each centre###########################'+'\n'
    self.write(line)
    for No in range(len(centreAndAroundNo)):
        
        centreSetNo=centreAndAroundNo[No]

        entropySet=WWOentropy[No]
        centreNo=centreSetNo[0]
        entropy=sum(entropySet)
        line=str(centreSetNo)[1:-1]+'    '+str(entropySet)[1:-1]+'\n'
        self.write(line)
    self.close()


centreCoord=getCentreCoord(centreFile)
print 'there were in total',len(centreCoord),'coorditates'
waters=getWaterInforWithCentre(waterFile)
print 'there were in total',len(waters),'waters'
waterInforCentre=orangeWaterInforWithCentre(waters,centreCoord)
print len(waterInforCentre)

removedCentre=removeLessOcuCentre(centreCoord,waterInforCentre,frameNo)

print len(removedCentre)
removedWaterInfor=removeLessOcuWaterInfor(centreCoord,waterInforCentre,frameNo)
print len(removedWaterInfor)

waterInforAroundCentre=findNearByWater(centreCoord,centreCoord,waterInforCentre,waterInforCentre)
waterCoordAroundCentre=extractWaterCoorAround(waterInforAroundCentre)
nearbyCentre=findNearbyCentre(centreCoord,centreCoord)
WWCentreNo=findNearbyCentreNo(centreCoord,centreCoord)

WWCentreAndAroundNo=findNearbyCentreAndAroundNo(centreCoord,centreCoord)


print 'the centre for WW_transloaction calculation are:',WWCentreNo


rdfOO=getRdfOO(rdfFile)

GRConstantEachCentre=calculateGRConstant(nearbyCentre,rdfOO)


allH2Os=getWaterInforWithCentre(allH2OInforFile)
allH2ORemovedTitle=allH2Os[1:]
allH2OInforCentre=orangeWaterInforWithCentre(allH2ORemovedTitle,centreCoord)
allOrganizedH2O=orangeH2OInfor(allH2OInforCentre)



##conservedH2Os=getWaterInforWithCentre(conservedH2OInforFile)
##conservedH2ORemovedTitle=conservedH2Os[1:]
##conservedH2OInforCentre=orangeWaterInforWithCentre(conservedH2ORemovedTitle,removedCentre)
##conservedOrganizedH2O=orangeH2OInfor(conservedH2OInforCentre)

H2OInforForWW=H2OInforNearbyHighCentre(centreCoord,centreCoord,allOrganizedH2O,allOrganizedH2O)
H2OCoordForWW=getH2OCoordForWW(H2OInforForWW)

##gEtheta=calculateEulerThetaReachCentre(H2OCoordForWW,nearbyCentre)
##gEphi=calculateEulerPhiReachCentre(H2OCoordForWW,nearbyCentre)
##gEpsi=calculateEulerPsiReachCentre(H2OCoordForWW,nearbyCentre)
##gr=calculateGReachCentre(waterCoordAroundCentre,nearbyCentre,3000)
##grtheta=calculateGThetaeachCentre(waterCoordAroundCentre,nearbyCentre)
##grphi=calculateGPhieachCentre(waterCoordAroundCentre,nearbyCentre)
##print GRConstantEachCentre

WforWW=getWforWW(GRConstantEachCentre,nearbyCentre,waterCoordAroundCentre,frameNo)
print WforWW
distance=getCentreAroundDistance(nearbyCentre)
        



##WWAngle=getWWAngle(H2OInforForWW,angleFile)
##SorWW=calculateSorWW(WWAngle,WforWW,frameNo)

print WWCentreNo

##print SorWW
print WforWW
print distance

##getSorWW(WWCentreNo,WforWW,WforWW,SorOutput)

thetaCos=getThetaCosDis()
thetaSin=getThetaSinDis()
angleCos=getAngleCosDis()
angleSin=getAngleSinDis()
print thetaCos
print thetaSin
print angleCos
print angleSin

print WWCentreNo
allSmWW=getSmWW(distance,waterCoordAroundCentre,nearbyCentre,H2OCoordForWW,frameNo)
print allSmWW

allSorWW=getSorWW(WWCentreNo,allSmWW,WforWW)

printOutWWOrien(SorOutput,WWCentreAndAroundNo,allSorWW)

