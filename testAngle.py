import sys
import os
import numpy
from Numeric import *
from numpy import sin,cos
waterFile=open('testH2OCoord.txt','r')


def readFileCoord(self):
    allCoord=[]
    H2OCoord=[]
    each=[]
    
    for line in self:
        line=line.split()
        coord=line[3:-3]
##        print coord
        realCoord=[]
        for i in coord:
##            print i
            i=float(i)
            realCoord.append(i)
##        realCoord=numpy.array(realCoord)
        H2OCoord.append(realCoord)
    print len(H2OCoord)

    no=1   
    for H2ONo in range(len(H2OCoord)):
        atomCoord=H2OCoord[H2ONo]
        print atomCoord

        if no/3.0==1.0:
            each.append(atomCoord)
            allCoord.append(each)
            each=[]
            no=1
        else:
            each.append(atomCoord)
            no+=1

    return allCoord
        
    

def calculateEular(self):
    referenceH1=numpy.array([0.79079641377315202, 0.61207926934631729, 0.0])
    referenceH2=numpy.array([-0.79079641377315202, 0.61207926934631729, 0.0])
    referenceZ=numpy.array([0.0,0.0,1.0])
    referenceY=numpy.array([0.0,1.0,0.0])
    revReferenceY=numpy.array([0.0,-1.0,0.0])
    H1=numpy.array(self[1])
    
    H2=numpy.array(self[2])

    OT=numpy.array(self[0])
            
    H1=H1-OT
    H2=H2-OT
    print H1
    print H2
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


    ZH2OOnXY=numpy.copy(ZH2O)
    ZH2OOnXY[2]=0.0
    ZH2OOnXYModulus=numpy.sqrt((ZH2OOnXY*ZH2OOnXY).sum())
    dot1=numpy.dot(revReferenceY,ZH2OOnXY)
    cosAnglePhi=dot1/ZH2OOnXYModulus
    anglePhi=numpy.arccos(cosAnglePhi)
    if ZH2OOnXY[0]<0:
        anglePhi=2*pi-anglePhi

    ZH2O=numpy.cross(H1,H2)
    ZH2OModulus=numpy.sqrt((ZH2O*ZH2O).sum())
    dot1=numpy.dot(referenceZ,ZH2O)
    cosAngleTheta=dot1/ZH2OModulus
    angleTheta=numpy.arccos(cosAngleTheta)

##    print 'The H2O coordidate is :', self
    
    print 'the angle is:',angleTheta,anglePhi,anglePsi
    a11=cos(anglePsi)*cos(anglePhi)-cos(angleTheta)*sin(anglePhi)*sin(anglePsi)
    a12=cos(anglePsi)*sin(anglePhi)+cos(angleTheta)*cos(anglePhi)*sin(anglePsi)
    a13=sin(anglePsi)*sin(angleTheta)
    a21=-sin(anglePsi)*cos(anglePhi)-cos(angleTheta)*sin(anglePhi)*cos(anglePsi)
    a22=-sin(anglePsi)*sin(anglePhi)+cos(angleTheta)*cos(anglePhi)*cos(anglePsi)
    a23=cos(anglePsi)*sin(angleTheta)
    a31=sin(angleTheta)*sin(anglePhi)
    a32=-sin(angleTheta)*cos(anglePhi)
    a33=cos(angleTheta)
    rotation=numpy.matrix([[a11,a12,a13],[a21,a22,a23],[a31,a32,a33]])
    print rotation
    H1=H1.reshape((-1,1))
    H2=H2.reshape((-1,1))
    
    print 'H1',H1
    H1=numpy.matrix(H1)
##    print H1
##    print rotation*(rotation**-1)
    H1=((rotation**-1)*(referenceH1.reshape(-1,1))).reshape(1,-1)
    H1=numpy.array(H1)
    print H1[0]
    print 'H2',H2
    H2=((rotation**-1)*(referenceH2.reshape(-1,1))).reshape(1,-1)
    print H2[0]
    H2=numpy.array(H2)
    move=numpy.array([1,1,1])
    
    O=numpy.array([0,0,0])
    H2O=numpy.array([O,H1[0],H2[0]])
    print H2O+move
    
##    print numpy.dot(rotation,H2)

coord=readFileCoord(waterFile)
print len(coord)


for i in coord:
    print len(i)
    calculateEular(i)
    
