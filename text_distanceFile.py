
def strCoordToFlost(self):
    coordinate=[]
    for eachCoord in self:
        eachCoord=float(eachCoord)
        coordinate.append(eachCoord)
    return coordinate

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
        print line
        line=line.split()
        print line
        line=strCoordToFlost(line)
        correlation.append(line)

        t1No=int(line[1])
        t2No=int(line[2])
        phiNo=int(line[3])
        c1No=int(line[4])
        c2No=int(line[5])
        value=line[6]
        print t1No,t2No,phiNo,c1No,c2No,value
        t1[t1No][t2No][phiNo][c1No][c2No]=value
    return t1

distanceFile='./eular2corr/distance_17.txt'

angle=openBulkCoorFile(distanceFile)

print angle[3][8][17][4][17]



    
