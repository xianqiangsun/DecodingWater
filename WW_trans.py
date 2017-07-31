############################################################################
##This file is used to find the waters near by the centre which we concerned
##to calculate the water-water translational entropy. The centres we concerned
##with is occupied with 95% of all the snapshots. Then the conserved waters
##around the 95% occupied waters around 5 Anstron is saved to the file. For
##for convinent. we just save the oxygen atoms for water-water translational
##entropy calculation.
##
##
##Xianqiang Sun
##TheoChem&Bio
##KTH
##2012-05-16
###########################################################################


import numpy
import sys
sys.path.append('/home/x/xiansu/pfs/program/numpy/lib/python2.6/site-packages')

from Numeric import *
from datetime import datetime


centreFile=open('optimizedCentre.txt','r')
waterFile=open('waterInforByCentre.txt','r')
H2OInforFile=open('WW_allCentre_H2O.txt','r')
rdfFile=open('rdf_HO_HO.xvg','r')
WWtransEntropyFile=open('WW_Trans.dat','w')

k=1.380648813*(10**(-23))   ## The unit of the boltzmann constant is J/K.
weiH2O=18.0154
mol=6.02214179*(10**23)
pw=0.0331725               ## The unit of this is No. of molecules in per A2.
frameNo=6000


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

##getCentreCoord is a function which use the input file of centre. Then the centres were splited
##to each centre. The output is a list of centres.

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
                print allCentreNo
##        print no,'waters arround centre',centreNo
        if len(eachCentreAndAround)>1:
            
            WWCentreNo.append(centreNo)
##    print WWCentreNo
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
                if nextNumber*pi/20>angle>=number*pi/20:
##                    print 'the angle of Phi is',waterCoor, centreCoor, angle,   number
                    grThetaNo[number]+=1
                elif angle==pi:
                    grThetaNo[19]+=1
##        print 'the angle distri bution',grThetaNo    
        for number in range(len(grThetaNo)):
            nextNumber=number+1
            grThetaNo[number]=(grThetaNo[number]/(cos(number*pi/20)-cos(nextNumber*pi/20)))/(frameNo/2.0)
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
##[[g(0~*pi/20),g(pi/20~2*pi/20)...g(pi~19*pi/20)]......]

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
                if nextNumber*pi/20>angle>=number*pi/20:
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
##[[g(0~*pi/20),g(pi/20~2*pi/20)...g(pi~19*pi/20)]......]

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

def intgGr(self):
    sumGr=0
    for eachGrNo in range(len(self)):
        eachGr=self[eachGrNo]
        nextNo=eachGrNo+1
        eachGr=eachGr*(((nextNo*0.05)**3)-((eachGrNo*0.05)**3))/3
        sumGr+=eachGr
##    print sumGr
    return sumGr

##intgGr is used to do the integration of GR according to the integral of 0.05. The input is a list a GR
##along the distance R. The output of this function is sum(gr*dr)

def intgGtheta(self):
    sumGtheta=0
    for gthetaNo in range(len(self)):
        eachGtheta=self[gthetaNo]
        angle1=gthetaNo*pi/20
        angle2=(gthetaNo+1)*pi/20
        eachGtheta=eachGtheta*(cos(angle1)-cos(angle2))
        sumGtheta+=eachGtheta
    return sumGtheta

##intgGtheta is used to do the integration of Gtheta according to the integral of pi/20. The input is a list a Gtheta
##along each angle. The output of this function is sum(gtheta*sinTheta*dTheta).



def intgGangle(self):
    sumGangle=0
    for eachGangle in self:
        eachGangle=eachGangle*pi/20
        sumGangle+=eachGangle
    return sumGangle

##intgGangle is used to do the integration of Gangle according to the integral of pi/20. The input is a list a Gangle
##along each angle. The output of this function is sum(gangle*dangle)


def intgGrln(self):
    sumGrln=0
    for eachGrlnNo in range(len(self)):
        eachGrln=self[eachGrlnNo]
        if eachGrln!=0.0:
            nextNo=eachGrlnNo+1
            eachGrln=numpy.log(eachGrln)*eachGrln*(((nextNo*0.05)**3)-((eachGrlnNo*0.05)**3))/3
            sumGrln+=eachGrln
    
    return sumGrln

##intgGrln is used to do the integration of GR according to the integral of 0.05. The input is a list a GR
##along the distance R. The output of this function is sum(gr*ln(gr)*dr)



def intgGthetaln(self):
    sumGthetaln=0
    for eachGthetalnNo in range(len(self)):
        eachGthetaln=self[eachGthetalnNo]
        if eachGthetaln!=0.0:
            angle1=eachGthetalnNo*pi/20
            angle2=(eachGthetalnNo+1)*pi/20
            eachGthetaln=numpy.log(eachGthetaln)*eachGthetaln*(cos(angle1)-cos(angle2))
            sumGthetaln+=eachGthetaln
    return sumGthetaln

##intgGthetaln is used to do the integration of Gthetaln according to the integral of pi/20. The input is a list a Gtheta
##along each angle. The output of this function is sum(gtheta*ln(gtheta)*sinTheta*dTheta).
        
    

def intgGangleln(self):
    sumGangleln=0
    for eachGangleln in self:
        if eachGangleln!=0.0:
            
            eachGangleln=numpy.log(eachGangleln)*eachGangleln*pi/20
            sumGangleln+=eachGangleln
    return sumGangleln

##intgGangleln is used to do the integration of Gangle according to the integral of pi/20. The input is a list a Gangle
##along each angle. The output of this function is sum(gangle*ln(gangle)*dangle)

def WWTrans(self,gr,gtheta,gphi):
    WWtransEntropy=[]
##    print 'constantis',self
##    print 'gr is',gr
##    print 'gtheat is',gtheta
##    print 'gphi is',gphi
    for setNo in range(len(self)):
        WWTransEachCentre=[]
        constantSet=self[setNo]
        grSet=gr[setNo]
        gthetaSet=gtheta[setNo]
        gphiSet=gphi[setNo]

        centreGIntg=intgGr(grSet[0])*intgGtheta(gthetaSet[0])*intgGangle(gphiSet[0])
##        WWTransSum=0
        for waterNo in range(len(constantSet)):
            
            eachConstant=constantSet[waterNo]
##            print eachConstant
            eachGr=grSet[waterNo]
            eachGtheta=gthetaSet[waterNo]
            eachGphi=gphiSet[waterNo]
##	    print 'centreGIntg*eachConstant*intgGr(eachGr)*intgGtheta(eachGtheta)*intgGangle(eachGphi)*mol'
##	    print centreGIntg,eachConstant,intgGr(eachGr),intgGtheta(eachGtheta),intgGangle(eachGphi)
            WWEachWater=(-0.5)*k*(pw**2)*centreGIntg*eachConstant*intgGr(eachGr)*intgGtheta(eachGtheta)*intgGangle(eachGphi)*mol/4.184
            WWTransEachCentre.append(WWEachWater)
##            WWTransSum+=WWEachWater
        WWtransEntropy.append(WWTransEachCentre)

    return WWtransEntropy

##WWTrans is used to calculate the final result of WW translocation entropy. The input incldes self: which is the constant between for g(R)
##actrually,We obtained it from calculation calculateGRInhTrans (It can also use the bulk g(R) for calculation).gr is the output of
##calculateGReachCentre. gtheta is the output of calculateGThetaeachCentre. gphi is the output of calculateGPhieachCentre. The output of
##this funciton included a list which saved the WW translocation entropy of each water centre. The output consistent with the watre centre No
##output of findNearbyCentreNo.
def printOutWWTrans(self,centreAndAroundNo,WWTentropy):
    line='#CentreNo'+'  '+'WWTransEntropy'+'\n'
    self.write(line)
    
    for No in range(len(centreAndAroundNo)):
        centreSetNo=centreAndAroundNo[No]

        entropySet=WWTentropy[No]
        centreNo=centreSetNo[0]
        entropy=sum(entropySet)
        line=str(centreNo)+'    '+str(entropy)+'\n'
        self.write(line)
    line='################explaination of each centre###########################'+'\n'
    self.write(line)
    for No in range(len(centreAndAroundNo)):
        
        centreSetNo=centreAndAroundNo[No]

        entropySet=WWTentropy[No]
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
##for i in waterCoordAroundCentre:
##    print '@@@@@@@@@@@@@@@',len(i)
##    for j in i:
##        print len(j)

nearbyCentre=findNearbyCentre(centreCoord,centreCoord)

WWCentreNo=findNearbyCentreNo(centreCoord,centreCoord)
WWCentreAndAroundNo=findNearbyCentreAndAroundNo(centreCoord,centreCoord)
print 'the centre for WW_transloaction calculation are:',WWCentreNo

rdfOO=getRdfOO(rdfFile)

##for eachWater in waterCoordAroundCentre:
##
##    eachCentre=nearbyCentre[eachCentreNo]
####    print eachCentre
##    print 'the GR distribution:',calculateGR(eachWater,eachCentre,3000)
##    print 'the Gtheta distribution:',calculateGRTheta(eachWater,eachCentre)
##    print 'the GPhi distribution:',calculateGRPhi(eachWater,eachCentre)
####    
##    eachCentreNo+=1
GrEachCentre=calculateGReachCentre(waterCoordAroundCentre,nearbyCentre,frameNo)

##for i in waterCoordAroundCentre:
##    print '@@@@@@@@@@@@@@@',len(i)
##    for j in i:
##        print len(j)

GRConstantEachCentre=calculateGRInhTrans(nearbyCentre,rdfOO)
GthetaEachCentre=calculateGThetaeachCentre(waterCoordAroundCentre,nearbyCentre,frameNo)
GphiEachCentre=calculateGPhieachCentre(waterCoordAroundCentre,nearbyCentre,frameNo)
print WWCentreNo

WWTransEntropy=WWTrans(GRConstantEachCentre,GrEachCentre,GthetaEachCentre,GphiEachCentre)

print WWTransEntropy
WWtransEntropyFile
printOutWWTrans(WWtransEntropyFile,WWCentreAndAroundNo,WWTransEntropy)
print 'print GRConstantEachCentre:',GRConstantEachCentre

    
    


                
                
            
            
        
