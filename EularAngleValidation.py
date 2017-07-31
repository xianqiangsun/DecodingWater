############################################################################
##This file is used to calculate the translational entropy with respect to
##the Protein of each water. The first is calculate the g(r) distribution.
##Then, the g(theta) and g(phi) is obtained respectively.
##
##
##Xianqiang Sun
##TheoChem&Bio
##KTH
##2012-05-16
##R2:2012-05-30
###########################################################################

import numpy
from Numeric import *
from datetime import datetime


centreFile=open('optimizedCentre.txt','r')
waterFile=open('waterInforByCentre.txt','r')
H2OInforFile=open('WW_allCentre_H2O.txt','r')
k=1.380648813*(10**(-23))   ## The unit of the boltzmann constant is J/K.
weiH2O=18.0154
mol=6.02214179*(10**23)
pw=0.0331725               ## The unit of this is No. of molecules in per A2.  
frameNo=3000 

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
##        print 'There were',len(H2OSetInfor), 'in H2OSetInfor!'
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

        if (float(ocupyNo)/float(frameNo))>=0.70:
            highOcuCentre.append(self[i])
    return highOcuCentre

def removeLessOcuWaterInfor(self,waterInforAccCen,frameNo):
    highOcuWaterInfor=[]
    for i in range(len(self)):
        ocupyNo=len(waterInforAccCen[i])
        
        if (float(ocupyNo)/float(frameNo))>=0.70:
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
            
        for waterCoor in waterCoorSet:
            waterCoor=numpy.array(waterCoor)
            dist=numpy.linalg.norm(waterCoor-centreCoor)
            
            for number in range(24):
                nextNumber=number+1
                
                if nextNumber*0.05>dist>=number*0.05:
                    grNo[number]+=1
                elif dist==24*0.05:
                    grNo[23]+=1
        print 'the distribution of GrNo',centreNo, grNo
        for number in range(len(grNo)):
            prevNumber=number+1
            grNo[number]=round(grNo[number]/(constant1*(((prevNumber*0.05)**3)-((number*0.05)**3)))/waterDensPerA2/frameNo,3)
        print 'the distribution of GrNo',centreNo, grNo    
        gr.append(grNo)
##    print gr
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
##        frameNo=len(waterCoorSet)
        print frameNo
        thetaDensity=frameNo/20.0
        for no in range(len(grThetaNo)):
            grThetaNo[no]=0
        for waterCoor in waterCoorSet:
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
                    grThetaNo[number]+=1
                elif angle==pi:
                    grThetaNo[19]+=1
        print 'the angle distri bution grThetaNo',centreNo,grThetaNo    
        for number in range(len(grThetaNo)):
            nextNumber=number+1
            grThetaNo[number]=(grThetaNo[number]/(cos(number*pi/20)-cos(nextNumber*pi/20)))/(frameNo/2.0)
        grTheta.append(grThetaNo)
        print 'the angle distri bution grThetaNo',centreNo, grThetaNo
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

def getH2OEulerTheta(self,frameNo):
    referenceH1=numpy.array([0.79079641377315202, 0.61207926934631729, 0.0])
    referenceH2=numpy.array([-0.79079641377315202, 0.61207926934631729, 0.0])
    referenceZ=numpy.array([0.0,0.0,1.0])
    eulerTheta=[]
##    eulerThetaDensity=frameNo/20.0
    CentreNo=0
    for H2OCoordSet in self:
        eulerThetaNo=range(20)
#        frameNo=len(H2OCoordSet)
        eulerThetaDensity=frameNo/20.0
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
            for number in range(20):
                nextNumber=number+1
                if nextNumber*pi/20>angleTheta>=number*pi/20:
                    eulerThetaNo[number]+=1
                elif angleTheta==pi:
                    eulerThetaNo[19]+=1
        print 'the euler Theta No is:',CentreNo, eulerThetaNo
        for number in range(len(eulerThetaNo)):
            nextNumber=number+1
            eulerThetaNo[number]=(eulerThetaNo[number]/(cos(number*pi/20)-cos(nextNumber*pi/20)))/(frameNo/2.0)
        eulerTheta.append(eulerThetaNo)
        print 'the euler Theta No is:',CentreNo,eulerThetaNo
    return eulerTheta

##getH2OEulerTheta is used to calculate the Euler distrubution of all the waters.
##The water coordinates were saved as lists with both Oxygen, and Hydeogen atoms.
##The coordinates looks like[[[OxygenCoord],[Hydgrogen1 coordinate],[Hydrogen2 Coordinate]....]...]
##the output of the Euler angle is very similar with other functions. the distribution is
##divided into 20 intergration spaces. The output is [[g(0~*pi/20),g(pi/20~2*pi/20)...g(pi~19*pi/20)]......]

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
        eulerPhiNo=range(40)
#        frameNo=len(H2OCoordSet)
        eulerPhiDensity=frameNo/40.0
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
            
            for number in range(40):
                nextNumber=number+1
                if nextNumber*pi/20>anglePhi>=number*pi/20:
                    eulerPhiNo[number]+=1
                elif anglePhi==2*pi:
                    eulerPhiNo[39]+=1
        print 'the euler phi No is:',CentreNo, eulerPhiNo
        for number in range(len(eulerPhiNo)):
            eulerPhiNo[number]=eulerPhiNo[number]/eulerPhiDensity
        eulerPhi.append(eulerPhiNo)
        print 'the euler phi No is:',CentreNo, eulerPhiNo
    return eulerPhi

##getH2OEulerPhi is used to calculate the Euler distrubution of all the waters.
##The water coordinates were saved as lists with both Oxygen, and Hydeogen atoms.
##The coordinates looks like[[[OxygenCoord],[Hydgrogen1 coordinate],[Hydrogen2 Coordinate]....]...]
##the output of the Euler angle is very similar with other functions. the distribution is
##divided into 20 intergration spaces. The output is [[g(0~*pi/20),g(pi/20~2*pi/20)...g(pi~19*pi/20)]......]

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
        eulerPsiDensity=frameNo/40.0
        eulerPsiNo=range(40)
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
            for number in range(40):
                nextNumber=number+1
                if nextNumber*pi/20>anglePsi>=number*pi/20:
                    eulerPsiNo[number]+=1
                elif anglePsi==2*pi:
                    eulerPsiNo[39]+=1
        print 'the euler Psi No is:',CentreNo, eulerPsiNo
        for number in range(len(eulerPsiNo)):
            eulerPsiNo[number]=eulerPsiNo[number]/eulerPsiDensity
        eulerPsi.append(eulerPsiNo)
        print 'the euler Psi No is:',CentreNo, eulerPsiNo
    return eulerPsi

##getH2OEulerPsi is used to calculate the Euler distrubution of all the waters.
##The water coordinates were saved as lists with both Oxygen, and Hydeogen atoms.
##The coordinates looks like[[[OxygenCoord],[Hydgrogen1 coordinate],[Hydrogen2 Coordinate]....]...]
##the output of the Euler angle is very similar with other functions. the distribution is
##divided into 20 intergration spaces. The output is [[g(0~*pi/20),g(pi/20~2*pi/20)...g(pi~19*pi/20)]......]

def intgGr(self):
    sumGr=0
    for eachGrNo in range(len(self)):
        eachGr=self[eachGrNo]
        nextNo=eachGrNo+1
        eachGr=eachGr*(((nextNo*0.05)**3)-((eachGrNo*0.05)**3))/3
        sumGr+=eachGr
    print  'the intergration of GR sum is:',sumGr
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


def calculateSWTrans(self,gtheta,gphi):
    SWTrans=[]
    grIntg=[]
    gthetaIntg=[]
    gphiIntg=[]
    grlngrIntg=[]
    gtlngtIntg=[]
    gplngpIntg=[]
        
    for grSet in self:
        grIntgSet=intgGr(grSet)
        grIntg.append(grIntgSet)
    print grIntg 
    for gthetaSet in gtheta:
        gthetaIntgSet=intgGtheta(gthetaSet)
        gthetaIntg.append(gthetaIntgSet)
##    print gthetaIntg
    for gphiSet in gphi:
        gphiIntgSet=intgGangle(gphiSet)
        gphiIntg.append(gphiIntgSet)
##    print gphiIntg
    for grSet in self:
        grlngrIntgSet=intgGrln(grSet)
        grlngrIntg.append(grlngrIntgSet)
##    print grlngrIntg
    for gthetaSet in gtheta:
        gtlngtIntgSet=intgGthetaln(gthetaSet)
        gtlngtIntg.append(gtlngtIntgSet)
##    print gtlngtIntg
    for gphiSet in gphi:
        gplngpIntgSet=intgGangleln(gphiSet)
        gplngpIntg.append(gplngpIntgSet)
##    print'the ln is', gplngpIntg
    for No in range(len(grIntg)):
        eachSWTrans=grlngrIntg[No]*gthetaIntg[No]*gphiIntg[No]+grIntg[No]*gtlngtIntg[No]*gphiIntg[No]+grIntg[No]*gthetaIntg[No]*gplngpIntg[No]
        eachSWTrans=(-1)*k*pw*eachSWTrans*mol/4.184
        SWTrans.append(eachSWTrans)
    
        
    
    return SWTrans

##calculateSWTrans is used to calculate SW translocation entropy of each water centre. The input of this function uses the output of
##calculateGr, calculateGRTheta,calculateGRPhi respecticely. With each input, the integration was performed with several
##other functions, such as intgGr,intgGtheta. The output of this function is a list which is saved trans entropy according
##to each centre( the output of extractWaterCoor).

def calculateSWOrien(self,eulerPhi,eulerPsi,waterCoord,frameNo):
    SWOrien=[]
    euTIntg=[]
    euPhIntg=[]
    euPsIntg=[]
    eTlneTIntg=[]
    ePhlnePhIntg=[]
    ePslnePsIntg=[]
    omega=8*pi*pi

    for euTset in self:
        euTIntgSet=intgGtheta(euTset)
        euTIntg.append(euTIntgSet)
##    print euTIntg
    for euPhset in eulerPhi:
        euPhIntgSet=intgGangle(euPhset)
        euPhIntg.append(euPhIntgSet)
##    print euPhIntg
    for euPsset in eulerPsi:
        euPsIntgSet=intgGangle(euPsset)
        euPsIntg.append(euPsIntgSet)
##    print euPsIntg
    for euTset in self:
        eTlneTIntgSet=intgGthetaln(euTset)
        eTlneTIntg.append(eTlneTIntgSet)
##    print eTlneTIntg
    for euPhset in eulerPhi:
        ePhlnePhIntgSet=intgGangleln(euPhset)
        ePhlnePhIntg.append(ePhlnePhIntgSet)
##    print ePhlnePhIntg
    for euPsset in eulerPsi:
        ePslnePsIntgSet=intgGangleln(euPsset)
        ePslnePsIntg.append(ePslnePsIntgSet)
##    print ePslnePsIntg
    for No in range(len(euTIntg)):
        eachSWOrien=eTlneTIntg[No]*euPhIntg[No]*euPsIntg[No]+euTIntg[No]*ePhlnePhIntg[No]*euPsIntg[No]+euTIntg[No]*euPhIntg[No]*ePslnePsIntg[No]
        
##        print NVi
        eachSWOrien=((-1)*k*len(waterCoord[No])/frameNo/omega)*eachSWOrien*mol/4.184
        SWOrien.append(eachSWOrien)
##        print eachSWOrien
    return SWOrien
##    print SWOrien
##?????????????????????????problem in determining omega
##calculateSWOrien can be used to calculate the binding entropy of SW orientation. The input includes Euler distribution of theta, Phi, and Psi.
##The output of this function is a list which is saved trans entropy according
##to each centre( the output of extractWaterCoor).


centreCoord=getCentreCoord(centreFile)
print 'there were in total',len(centreCoord),'coorditates'
waters=getWaterInforWithCentre(waterFile)
print 'there were in total',len(waters),'waters'
waterInforCentre=orangeWaterInforWithCentre(waters,centreCoord)
print len(waterInforCentre)
removedCentre=removeLessOcuCentre(centreCoord,waterInforCentre,frameNo)

print len(removedCentre)
removedWaterInfor=removeLessOcuWaterInfor(centreCoord,waterInforCentre,frameNo)
##print len(removedWaterInfor)
##for i in removedWaterInfor:
##    print len(i)

waterCoor=extractWaterCoor(waterInforCentre)
print len(waterCoor)
allWaterNo=0
##for i in waterCoor:
##    print len(i)
##    allWaterNo=allWaterNo+len(i)
##print 'there were totaly',allWaterNo,'waters'
grDist=calculateGR(waterCoor,centreCoord,frameNo)

##print grDist
##print len(grDist)

grTheta=calculateGRTheta(waterCoor,centreCoord,frameNo)
##print grTheta

grPhi=calculateGRPhi(waterCoor,centreCoord,frameNo)
##print grPhi

H2Os=getWaterInforWithCentre(H2OInforFile)
H2ORemovedTitle=H2Os[1:]
##Remove the title line for the calculation
print len(H2ORemovedTitle)

H2OInforCentre=orangeWaterInforWithCentre(H2ORemovedTitle,centreCoord)

##print len(H2OInforCentre)
##for i in H2OInforCentre:
##    print len(i)

organizedH2O=orangeH2OInfor(H2OInforCentre)
##print len(organizedH2O)
##for i in organizedH2O:
##    for j in i:
##        print 'there were',len(j),'waters in one water'
##
##    print len(i)

allH2OCoord=getH2OCoord(organizedH2O)

eulerTheta=getH2OEulerTheta(allH2OCoord,frameNo)
eulerPhi=getH2OEulerPhi(allH2OCoord,frameNo)
eulerPsi=getH2OEulerPsi(allH2OCoord,frameNo)



SWTrans=calculateSWTrans(grDist,grTheta,grPhi)
print SWTrans

SWOrien=calculateSWOrien(eulerTheta,eulerPhi,eulerPsi,waterCoor,frameNo)
print SWOrien


    
            
        

            
        
    



        
    
