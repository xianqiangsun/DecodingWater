import sys
sys.path.append('/home/x/xiansu/pfs/program/numpy/lib/python2.6/site-packages')
from MDAnalysis import Universe, Writer
from MDAnalysis.analysis.distances import distance_array
import MDAnalysis
import numpy
from Numeric import *


top='npt.gro'
traj='md_extract1.trr'


water=Universe(top,traj)

o=water.selectAtoms('name O*')

resid=o.resids()
print resid

#resnu=o.resnums()
#resna=o.resnames()

atomInf=[]
for i in o.atoms:
    atomid= str(i).split()[2]
    atomseg=str(i).split()[-1]
    atomidandseg=[]
    atomidandseg.append(atomid)
    atomidandseg.append(atomseg)
    atomInf.append(atomidandseg)
print atomInf
##print len(waterResnu)

box = water.trajectory.ts.dimensions[:3]
print box

disNo=range(40)
anglet1=[]
anglet2=[]
angleph=[]
anglec1=[]
anglec2=[]
angle=[]


for i in range(40):
    disNo[i]=0
    angle.append([])
    anglet1.append([0,0,0,0,0,0,0,0,0,0])
    anglet2.append([0,0,0,0,0,0,0,0,0,0])
    angleph.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    anglec1.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    anglec2.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

##def rotationAngle(self,vector2,angle):
##    epsilon=10**(-5)
##    v1xy=[]
##    v2xy=[]
##    v1yz=[]
##    v2yz=[]
##    
##    v1xy.extend(self[:2])
##    v2xy.extend(vector2[:2])
##    v1yz.extend(self[1:])
##    v2yz.extend(vector2[1:])
##    
##    v1xy=numpy.array(v1xy)
##    v2xy=numpy.array(v2xy)
##    v1yz=numpy.array(v1yz)
##    v2yz=numpy.array(v2yz)
##    
##    
##    v1xyM=numpy.sqrt((v1xy*v1xy).sum())
##    v2xyM=numpy.sqrt((v2xy*v2xy).sum())
##    v1yzM=numpy.sqrt((v1yz*v1yz).sum())
##    v2yzM=numpy.sqrt((v2yz*v2yz).sum())
##    dotxy=numpy.dot(v1xy,v2xy)
##    dotyz=numpy.dot(v1yz,v2yz)
##    planeAngle=numpy.arccos(dotxy/v1xyM/v2xyM)
##
##    if abs(numpy.cross(v1xy,v2xy))<epsilon:
##        
##        if abs(numpy.cross(v1yz,v2yz))<epsilon:
##            anlge=0.0
##        elif abs(numpy.cross(v1yz,v2yz))>=epsilon:
##            yzCross=numpy.cross(v1yz,v2yz)
##            if yzCross<=0:
##                angle=2*pi-angle
####                print angle
##            
##    elif abs(numpy.cross(v1xy,v2xy))>=epsilon:
##        xyCross=numpy.cross(v1xy,v2xy)
##        if xyCross<=0:
##            angle=2*pi-angle
####            print angle
##    return angle

##rotationAngle is used to determine the rotation angle between two vectors. As we know that
##the angle obtained from the vector arccos is the angle between the two vectors,
##therefore it is just the value from 0-pi, if we want to determine the rotation angle
##we have to perform the other method. In this function, I use reflection the 3D vector
##to 2D vector to determine whether the angle is clockwise or anticlockwise. If the angle
##is anticlockwise, the angle becomes to 2*pi-angle. the input of the function included
##self: the first 3D vector for use. vector2 is the second vector. anlge is the
##3D angle between self and vector2. The output is the rotation angle between the
##two vectors.


def calculateWWangle(self,cCoord):
    angle=[]
##    print 'two coordidates are',self,cCoord
    vOCoord=numpy.array(self[0])
    vH1Coord=numpy.array(self[1])
    vH2Coord=numpy.array(self[2])
    
    cOCoord=numpy.array(cCoord[0])
    cH1Coord=numpy.array(cCoord[1])
    cH2Coord=numpy.array(cCoord[2])
    
    O1D1=(numpy.add(vH1Coord,vH2Coord))/2-vOCoord
##    print O1D1
    O1H11=vH1Coord-vOCoord
    O1H12=vH2Coord-vOCoord
    
    O2D2=(numpy.add(cH1Coord,cH2Coord))/2-cOCoord
##    print O2D2
    O2H21=cH1Coord-cOCoord
    O2H22=cH2Coord-cOCoord
    
    O1O2=cOCoord-vOCoord
##    print O1O2
    O2O1=vOCoord-cOCoord
##    print O2O1
    H11H12=vH2Coord-vH1Coord
    H21H22=cH2Coord-cH1Coord
    
    crossO1O2O1D1=numpy.cross(O1O2,O1D1)
##    print crossO1O2O1D1
    crossO2D2O2O1=numpy.cross(O2D2,O2O1)
##    print crossO1O2O1D1,crossO2D2O2O1
##                            print crossO2D2O2O1
                            

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
##    print dotPhi,crossO1O2O1D1,crossO2D2O2O1
    dotChi1=numpy.dot(H11H12,crossO1O2O1D1)
    dotChi2=numpy.dot(H21H22,crossO2D2O2O1)

    
    cosTheta1=dotTheta1/O1D1Modulus/O1O2Modulus
    cosTheta2=dotTheta2/O2D2Modulus/O2O1Modulus
    cosPhi=dotPhi/crossO1O2O1D1Modulus/crossO2D2O2O1Modulus
##    print cosPhi,dotPhi,crossO1O2O1D1Modulus,crossO2D2O2O1Modulus
    cosChi1=dotChi1/H11H12Modulus/crossO1O2O1D1Modulus
    cosChi2=dotChi2/H21H22Modulus/crossO2D2O2O1Modulus

#    print cosTheta1,cosTheta2,cosPhi,cosChi1,cosChi2
##    theta1=numpy.arccos(cosTheta1)
    if abs(cosTheta1-1.0)<=10**(-5):
        print 'theta1'
        theta1=0.0
    elif  abs(cosTheta1-1.0)>10**(-5):

        theta1=numpy.arccos(cosTheta1)
##        print theta1
    if abs(cosTheta2-1.0)<=10**(-5):
        
        print 'theta2'
        theta2=0.0
    elif abs(cosTheta2-1.0)>10**(-5):
        theta2=numpy.arccos(cosTheta2)
    
##    theta2=numpy.arccos(cosTheta2)

##    print cosPhi                        
##    phi=numpy.arccos(cosPhi)
    if abs(cosPhi-1.0)<=10**(-5):
        print 'phi',O1O2,O2O1
        phi=0
    elif abs(cosPhi-1.0)>=10**(-5):
        phi=numpy.arccos(cosPhi)
    print phi,O1O2,O2O1
##    print phi
        
        

##    phi=rotationAngle(crossO1O2O1D1,crossO2D2O2O1,phi)

    

    if abs(cosChi1-1.0)<=10**(-6):
        print 'cosChi1'
        chi1=0.0
    elif abs(cosChi1-1.0)>=10**(-6):
        chi1=numpy.arccos(cosChi1)
        
##    chi1=rotationAngle(H11H12,crossO1O2O1D1,chi1)
                            
    if abs(cosChi2-1.0)<=10**(-6):
        print 'chi2'
        chi2=0.0
    elif abs(cosChi2-1.0)>=10**(-6):
        chi2=numpy.arccos(cosChi2)

##    chi2=rotationAngle(H21H22,crossO2D2O2O1,chi2)
#    print theta1,theta2,phi,chi1,chi2
    angle.append(theta1)
    angle.append(theta2)
    angle.append(phi)
    angle.append(chi1)
    angle.append(chi2)
##    print angle
    return angle
     

def calculateGtt(self,angle):
    gttTotal=[]
 
    for disNo in range(len(self)):
        gtt=[]
        for i in range(10):
            gtt.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        distance=disNo/10.0+2
        
        waterNoDist=self[disNo]
        
        angleDist=angle[disNo]

        density=waterNoDist/4.0
        print 'the density is:', disNo, density

        if abs(density)<=10**(-8):
            gttTotal.append(gtt)
        elif abs(density)>=10**(-8):
            for pairNo in range(len(angleDist)):
                eachAngle=angleDist[pairNo]
                t1=eachAngle[0]
                t2=eachAngle[1]
                t1No=int(t1/(pi/10))
                t2No=int(t2/(pi/10))
                
                if t1No==10:
                    t1No=9
                if t2No==10:
                    t2No=9
                gtt[t1No][t2No]+=1
            for t1No in range(len(gtt)):
##            print gtt
##            print t1No
                gt1=gtt[t1No]
                for t2No in range(len(gt1)):
##                if density==0:
##                    gtt[t1No][t2No]=0.0
##                elif density!=0:

                    intgSpace=(cos(t1No*pi/10.0)-cos((t1No+1)*pi/10.0))*(cos(t2No*pi/10.0)-cos((t2No+1)*pi/10.0))
                    gtt[t1No][t2No]=gt1[t2No]/intgSpace/density
            gttTotal.append(gtt)
##        print gtt
                
        
    return gttTotal

def calculateGt2c1(self,angle):
    gt2c1Total=[]
 
    for disNo in range(len(self)):
        gt2c1=[]
        for i in range(10):
            gt2c1.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        distance=disNo/10.0+2
        
        waterNoDist=self[disNo]
        angleDist=angle[disNo]

        density=waterNoDist/2.0/pi

        if abs(density)<=10**(-8):
            gt2c1Total.append(gt2c1)
            
        elif abs(density)>=10**(-8):
            for pairNo in range(len(angleDist)):
                eachAngle=angleDist[pairNo]
                t2=eachAngle[1]
                c1=eachAngle[3]
                t2No=int(t2/(pi/10))
                c1No=int(c1/(pi/20))
                if t2No==10:
                    t2No=9
                if c1No==20:
                    c1No=19
                gt2c1[t2No][c1No]+=1
            for t2No in range(len(gt2c1)):
                gt2=gt2c1[t2No]
                for c1No in range(len(gt2)):
##                if density==0:
##                    gt2c1[t2No][c1No]=0.0
##                elif density!=0:
                    
                
                    
                    intgSpace=(cos(t2No*pi/10.0)-cos((t2No+1)*pi/10.0))*pi/20.0
                    gt2c1[t2No][c1No]=gt2[c1No]/intgSpace/density
                    
            gt2c1Total.append(gt2c1)
    return gt2c1Total


def calculateGt1c2(self,angle):
    gt1c2Total=[]
 
    for disNo in range(len(self)):
        gt1c2=[]
        for i in range(10):
            gt1c2.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        distance=disNo/10.0+2
        waterNoDist=self[disNo]
        angleDist=angle[disNo]

        density=waterNoDist/2.0/pi

        if abs(density)<=10**(-8):
            gt1c2Total.append(gt1c2)
        elif abs(density)>=10**(-8):
            for pairNo in range(len(angleDist)):
                eachAngle=angleDist[pairNo]
                t1=eachAngle[0]
                c2=eachAngle[4]
                t1No=int(t1/(pi/10))
                c2No=int(c2/(pi/20))
                if t1No==10:
                    t1No=9
                if c2No==20:
                    c2No=19
                gt1c2[t1No][c2No]+=1
            for t1No in range(len(gt1c2)):
                gt1=gt1c2[t1No]
                for c2No in range(len(gt1)):
##                if density==0:
##                    gt1c2[t1No][c2No]=0.0
##                elif density!=0:

                    intgSpace=(cos(t1No*pi/10.0)-cos((t1No+1)*pi/10.0))*pi/20.0
                    gt1c2[t1No][c2No]=gt1[c2No]/intgSpace/density
                
            gt1c2Total.append(gt1c2)
    return gt1c2Total

def calculateGc1c2(self,angle):
    gc1c2Total=[]
 
    for disNo in range(len(self)):
        gc1c2=[]
        for i in range(20):
            gc1c2.append([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        distance=disNo/10.0+2
        
        waterNoDist=self[disNo]
        angleDist=angle[disNo]

        density=waterNoDist/pi/pi

        if abs(density)<=10**(-8):
            gc1c2Total.append(gc1c2)
        elif abs(density)>=10**(-8):
            for pairNo in range(len(angleDist)):
                eachAngle=angleDist[pairNo]
                c1=eachAngle[3]
                c2=eachAngle[4]
                c1No=int(c1/(pi/20))
                c2No=int(c2/(pi/20))
                if c1No==20:
                    c1No=19
                if c2No==20:
                    c2No=19
                gc1c2[c1No][c2No]+=1
            for c1No in range(len(gc1c2)):
                gc1=gc1c2[c1No]
                for c2No in range(len(gc1)):
##                if density==0:
##                    gc1c2[c1No][c2No]=0.0
##                elif density!=0:
                
                    intgSpace=(pi/20.0)*(pi/20.0)
                    gc1c2[c1No][c2No]=gc1[c2No]/intgSpace/density
                
            gc1c2Total.append(gc1c2)
    return gc1c2Total

def calculateGt1(self,angle,position):
    gt1Total=[]
    for disNo in range(len(self)):
        gt1=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        distance=disNo/10.0+2
        waterNoDist=self[disNo]
        angleDist=angle[disNo]
        density=waterNoDist/2.0
        if abs(density)<=10**(-8):
            gt1Total.append(gt1)
        elif abs(density)>=10**(-8):
            for pairNo in range(len(angleDist)):
                eachAngle=angleDist[pairNo]
                t1=eachAngle[position]
                t1No=int(t1/(pi/10))
                if t1No==10:
                    t1No=9


                gt1[t1No]+=1
            for t1No in range(len(gt1)):
                intgSpace=cos(t1No*pi/10.0)-cos((t1No+1)*pi/10.0)
                gt1[t1No]=gt1[t1No]/intgSpace/density
            gt1Total.append(gt1)
    return gt1Total

def calculateGc(self,angle,position):
    gt1Total=[]
    for disNo in range(len(self)):
        gt1=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        distance=disNo/10.0+2
        waterNoDist=self[disNo]
        angleDist=angle[disNo]
        density=waterNoDist/pi
        if abs(density)<=10**(-8):
            gt1Total.append(gt1)
        elif abs(density)>=10**(-8):
            for pairNo in range(len(angleDist)):
                eachAngle=angleDist[pairNo]
                t1=eachAngle[position]
                t1No=int(t1/(pi/20))
##                print t1,t1No
                if t1No==20:
                    t1No=19
                    gt1[t1No]+=1
                else:

                    gt1[t1No]+=1
            for t1No in range(len(gt1)):
                intgSpace=pi/20.0
                gt1[t1No]=gt1[t1No]/intgSpace/density
            gt1Total.append(gt1)
    return gt1Total

def calculateGww(self,angle):
    gww=[]
    gttAll=calculateGtt(self,angle)
##    print gttAll
    gt1c2All=calculateGt1c2(self,angle)
    gt2c1All=calculateGt2c1(self,angle)
    gc1c2All=calculateGc1c2(self,angle)
    gt1All=calculateGt1(self,angle,0)
    gt2All=calculateGt1(self,angle,1)
    gphAll=calculateGc(self,angle,2)
    gc1All=calculateGc(self,angle,3)
    gc2All=calculateGc(self,angle,4)
##    fileName='distance'+'.txt'
##    output=open(fileName,'w')
    
    
    for disNo in range(len(self)):
        fileName='distance_'+str(disNo)+'.txt'
        output=open(fileName,'w')
##        
        distance=disNo*0.1+2.0
        waterNoDist=self[disNo]
##        angleDist=angle[disNo]
        gtt=gttAll[disNo]
##        print disNo,gtt
        gt1c2=gt1c2All[disNo]
        gt2c1=gt2c1All[disNo]
        gc1c2=gc1c2All[disNo]
        gt1=gt1All[disNo]
        gt2=gt2All[disNo]
        gph=gphAll[disNo]
        print disNo,gph
        gc1=gc1All[disNo]
        gc2=gc2All[disNo]
        eachDis=0
        for t1No in range(len(gtt)):
            sgt1=gt1[t1No] #means subgt1,1 exact value obtained
            sgtt=gtt[t1No]
            sgt1c2=gt1c2[t1No] 
            '''determines the theta 1 angle to be intergrated'''
            for t2No in range(len(sgtt)):
                sgt2=gt2[t2No] #means subgtt,2 exact value obtained
                ssgtt=sgtt[t2No] #means subsubgtt,3 exact value obtained
                sgt2c1=gt2c1[t2No]
                '''determines the theta 2 angle to be intergrated'''
                for c1No in range(len(sgt2c1)):
                    sgc1=gc1[c1No] #means subgtt,4 exact value obtained
                    ssgt2c1=sgt2c1[c1No] #means subgtt,5 exact value obtained
                    sgc1c2=gc1c2[c1No]
                    '''determines the chi 1 angle to be intergrated'''
                    for c2No in range(len(sgc1c2)):
                        sgc2=gc2[c2No] #means subgtt,6 exact value obtained
                        ssgc1c2=sgc1c2[c2No] #means sub sub gtt,7 exact value obtainedc2No
                        ssgt1c2=sgt1c2[c2No] #means sub sub gtt,8 exact value obtainedc2No
                        for phNo in range(len(gph)):
                            sgph=gph[phNo] #means sub sub gtt,9 exact value obtainedc2No, all exactvalue obtained to obtain exact WW g(R)

##                              grMult=ssgtt*ssgt1c2*ssgt2c1*ssgc1c2*sgph/(sgt1*sgt2*sgc1*sgc2)
                            if sgt1*sgt2*sgc1*sgc2==0:
                                line=str(disNo)+'   '+str(t1No)+'    '+str(t2No)+'    '+str(phNo)+'    '+str(c1No)+'    '+str(c2No)+'    '+'0'+'\n'
##                                output.write(line)
##                                outputFile.write(line)

                            elif sgt1*sgt2*sgc1*sgc2!=0:
                                grMult=ssgtt*ssgt1c2*ssgt2c1*ssgc1c2*sgph/(sgt1*sgt2*sgc1*sgc2)

                                if grMult==0:
                                    line=str(disNo)+'   '+str(t1No)+'    '+str(t2No)+'    '+str(phNo)+'    '+str(c1No)+'    '+str(c2No)+'    '+'0'+'\n'
##                                    output.write(line)
##                                    outputFile.write(line)
                                elif grMult!=0:
##				    print 'grMult is',grMult,t1No,t2No,c1No,c2No,phNo
                                    eachGrln=numpy.log(grMult)*grMult#*((pi/20.0)**3)*(cos(t1No*pi/20.0)-cos((t1No+1)*pi/20.0))*(cos(t2No*pi/20.0)-cos((t2No+1)*pi/20.0))
                                    
                                    line=str((disNo*0.1)+2)+'   '+str(t1No)+'    '+str(t2No)+'    '+str(phNo)+'    '+str(c1No)+'    '+str(c2No)+'    '+str(eachGrln)+'\n'
                                    output.write(line)
                                    
##                                    print eachGrln
##				    print 'grMult is grMult,t1No,t2No,c1No,c2No,phNo, ssgtt*ssgt1c2*ssgt2c1*ssgc1c2*sgph/(sgt1*sgt2*sgc1*sgc2), eachGrln'
##				    print grMult,t1No,t2No,c1No,c2No,phNo, ssgtt,ssgt1c2,ssgt2c1,ssgc1c2,sgph,sgt1,sgt2,sgc1,sgc2,eachGrln
##                                    eachDis+=eachGrln
        output.close()
        print eachDis
        eachDis=(1.0/(32*(pi**3)))*eachDis
        print eachDis
            
        gww.append(eachDis)
        distance=str(distance)+'	'+str(eachDis)+'\n'
  
        print eachDis,'@@@@@@@@@@@@@@@@@@@@@@@' ,distance
    print gww
##    output.write(distance)
##    output.close()
    
    return gww


for ts in water.trajectory:
    print ts
    oc1=o.coordinates()
    oc2=o.coordinates()
    d=distance_array(oc1,oc2,box)
        
    vwaterId=0
    for i in d:
        vwaterId+=1
        vresid='resid '+str(vwaterId)
        vH2O=water.selectAtoms(vresid)
        vCoord=vH2O.coordinates()
        cwaterId=0
        for j in i:
            cwaterId+=1
##            print j
            no=int((j-2.0)/0.1)
##            print no
            if 0<no<=39:
##                vresid='resid '+str(vwaterId)
                cresid='resid '+str(cwaterId)
##                vH2O=water.selectAtoms(vresid)
                cH2O=water.selectAtoms(cresid)
##                vCoord=vH2O.coordinates()
                cCoord=cH2O.coordinates()
##                print vwaterId,cwaterId,vCoord,cCoord
                pairAngle=calculateWWangle(vCoord,cCoord)
##                print no, pairAngle
                angle[no].append(pairAngle)
                


##                

                
                
##                print angle
                
                disNo[no]+=1
                

print disNo


calculateGww(disNo,angle)




            
            
        
        

    




