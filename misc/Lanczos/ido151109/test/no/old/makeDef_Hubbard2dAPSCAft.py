#!/usr/bin/env python
import sys,time,os,random
from math import *

T1={(1,0): 1.0, (-1,0): 1.0
}
T2={(0,1): 1.0, (0,-1): 1.0
}
T3={(1,1): -0.3, (-1,-1): -0.3
}
T4={(1,-1): -0.3, (-1,1): -0.3
}
U=4.0
V1={}
V2={}
J1={}
J2={}

eScale=1.0

TInit = 0.0
TEnd  = 10.0
Time1 = 5.0
Time2 = 1000.0
Time3 = 0.0

NStoreO = 1
NStoreM = 0


def position(r):
    x=(r[0]+Lx)%Lx
    y=(r[1]+Ly)%Ly
    return (x,y)

def indexToPosition(i):
    x=i%Lx
    y=i/Lx
    return (x,y)

def positionToIndex(r):
    x=(r[0]+Lx)%Lx
    y=(r[1]+Ly)%Ly
    return x+Lx*y

def neighborIndex(i,dr):
    r=indexToPosition(i)
    x=r[0]+dr[0]
    y=r[1]+dr[1]
    return positionToIndex([x,y])

def direction(i,j):
    ri=indexToPosition(i)
    rj=indexToPosition(j)
    dx=(rj[0]-ri[0]+Lx)%Lx
    dy=(rj[1]-ri[1]+Ly)%Ly
    return (dx,dy)

def distance(i,j):
    ri=indexToPosition(i)
    rj=indexToPosition(j)
    dx=(rj[0]-ri[0]+Lx)%Lx
    dy=(rj[1]-ri[1]+Ly)%Ly
    if dy > Ly/2: 
      dy = Ly-dy
    if dx > Lx/2: 
      dx = Lx-dx
    return sqrt(dx*dx+dy*dy)

def locgrnIdx(i,j,s):
    return (Nsite*i+j)+s*Nsite*Nsite

def subIndex(i):
    r=indexToPosition(i)
    sx=r[0]%Sx
    sy=r[1]%Sy
    return sx+Sx*sy

def sgnAP(i,dr):
    r=indexToPosition(i)
    x=r[0]+dr[0]
    if (x >= Lx or x <= -1):
        return -1
    else:
        return 1

if len(sys.argv)<7:
    print "./makeDefFile.ph Lx Ly Sx Sy Ne U"
    sys.exit()

Lx = int(sys.argv[1])
Ly = int(sys.argv[2])
Sx = int(sys.argv[3])
Sy = int(sys.argv[4])
Nelectron = int(sys.argv[5])
U = float(sys.argv[6])
lambdaScale = 1.0
muScale = 1.0 # 1D scale

NSplitSize = 1
# NSplitSize = 4
# NSplitSize = 8

#Sx = Lx
#Sy = Ly
Nsub = Sx*Sy 
Nsite = Lx*Ly 

time.sleep(1.5)
seed = int(time.time())
random.seed()

separator = '--------------------\n'

pre='z'
fileName = []
fileName.append('xnamelist.def')
fileName.append(pre+'modpara.def')
fileName.append(pre+'locspn.def')
fileName.append(pre+'transfer.def')
fileName.append(pre+'coulombintra.def')
fileName.append(pre+'coulombinter.def')
fileName.append(pre+'hund.def')
fileName.append(pre+'pairhop.def')
fileName.append(pre+'exchange.def')
fileName.append(pre+'gutzwilleridx.def')
fileName.append(pre+'jastrowidx.def')
fileName.append(pre+'doublonholon2siteidx.def')
fileName.append(pre+'doublonholon4siteidx.def')
fileName.append(pre+'orbitalidx.def')
fileName.append(pre+'qptransidx.def')
fileName.append(pre+'cisajs.def')
fileName.append(pre+'cisajscktalt.def')
fileName.append(pre+'cisajscktaltdc.def')
fileName.append(pre+'interall.def')
fileName.append(pre+'quench.def')

### xnamelist.def ###
f = open(fileName[0],'w')
for x in fileName[1:]:
    f.write(x+"\n")
f.close()

### modpara.def ###
f = open(fileName[1],'w')
f.write(
    separator+
    "Model_Parameters  0\n"+
    separator+
    "VMC_Cal_Parameters\n"+
    separator+
    "CDataFileHead  zvo\n"+
    "CParaFileHead  zqp\n"+
    separator+
    "NVMCCalMode    0\n"+
    "NLanczosMode   0\n"+
    separator+
    "NDataIdxStart  0\n"+
    "NDataQtySmp    1\n"+
    separator+
    "Nsite          {0}\n".format(Nsite)+
    "Nelectron      {0}\n".format(Nelectron)+
    "NSPGaussLeg    1\n"+
    "NSPStot        0\n"+
    "NMPTrans       -1\n"+
    "NSROptItrStep  4000\n"+
    "NSROptItrSmp   200\n"+
    "NSROptFixSmp   1\n"+
    "DSROptRedCut   1e-7\n"+
    "DSROptStaDel   4e-2\n"+
    "DSROptStepDt   2e-2\n"+
    "NVMCWarmUp     10\n"+
    "NVMCIniterval  1\n"+
    "NVMCSample     500\n"+
    "NExUpdatePath  0\n"+
    "RndSeed        {0}\n".format(seed)+
    "NSplitSize     {0}\n".format(NSplitSize)+
    "NStoreO        {0}\n".format(NStoreO)+
    "NStoreM        {0}\n".format(NStoreM))
f.close()

### locspn.def ###
f = open(fileName[2],'w')
f.write(separator+
        "NLocalSpin\t0\n"+
        separator+
        "i_0LocSpn_1IteElc\n"+
        separator)
for i in range(Nsite):
    f.write("{0}\t1\n".format(i))
f.close()

### transfer.def ###
f = open(fileName[3],'w')
f.write(separator+
        "NTransfer\t"+str(2*4*Nsite)+"\n"+
        separator+
        "i_j_s_tijs\n"+
        separator)
paraList = []
for i in range(Nsite):
    for dr,var in T1.iteritems():
        j = neighborIndex(i,dr)
        sgn = sgnAP(i,dr)
        paraList.append([i,j,sgn*var/eScale])
    for dr,var in T2.iteritems():
        j = neighborIndex(i,dr)
        sgn = sgnAP(i,dr)
        paraList.append([i,j,sgn*muScale*var/eScale])
    #for dr,var in T3.iteritems():
    #    j = neighborIndex(i,dr)
    #    sgn = sgnAP(i,dr)
    #    paraList.append([i,j,sgn*muScale*var/eScale])
    #for dr,var in T4.iteritems():
    #    j = neighborIndex(i,dr)
    #    sgn = sgnAP(i,dr)
    #    paraList.append([i,j,sgn*muScale*var/eScale])
paraList.sort()
for s in [0,1]:
    for para in paraList:
        f.write("{0}\t{1}\t{2}\t{3}\n".format(para[0],para[1], s, para[2]))
f.close()

### coulombintra.def ###
f = open(fileName[4],'w')
f.write(separator+
        "NCoulombIntra\t"+str(Nsite)+"\n"+
        separator+
        "i_Ui\n"+
        separator)
for i in range(Nsite):
    f.write("{0}\t{1}\n".format(i,lambdaScale*U/eScale))
f.close()

## V ##
paraList = []
for i in range(Nsite):
    for dr,var in V1.iteritems():
        j = neighborIndex(i,dr)
        paraList.append([i,j,lambdaScale*var/eScale])
    for dr,var in V2.iteritems():
        j = neighborIndex(i,dr)
        paraList.append([i,j,lambdaScale*muScale*var/eScale])
paraList.sort()

### coulombinter.def ###
f = open(fileName[5],'w')
f.write(separator+
        "NCoulombInter\t"+ str(len(paraList)/2) +"\n"+
        separator+
        "i_j_Vij\n"+
        separator)
for para in paraList:
    if para[0] < para[1]:
        f.write("{0}\t{1}\t{2}\n".format(para[0],para[1],para[2]))
f.close()

## J ##
paraList = []
for i in range(Nsite):
    for dr,var in J1.iteritems():
        j = neighborIndex(i,dr)
        paraList.append([i,j,lambdaScale*var/eScale])
    for dr,var in J2.iteritems():
        j = neighborIndex(i,dr)
        paraList.append([i,j,lambdaScale*muScale*var/eScale])
paraList.sort()

### hunt.def ###
f = open(fileName[6],'w')
f.write(separator+
        "NHundCoupling\t"+ str(len(paraList)/2) +"\n"+
        separator+
        "i_j_Jij\n"+
        separator)
for para in paraList:
    if para[0] < para[1]:
        f.write("{0}\t{1}\t{2}\n".format(para[0],para[1],para[2]))
f.close()

### pairhop.def ###
f = open(fileName[7],'w')
f.write(separator+
        "NPairHopping\t"+ str(len(paraList)) +"\n"+
        separator+
        "i_j_Jij\n"+
        separator)
for para in paraList:
    f.write("{0}\t{1}\t{2}\n".format(para[0],para[1],para[2]))
f.close()

### exchange.def ###
f = open(fileName[8],'w')
f.write(separator+
        "NExchangeCoupling\t"+ str(len(paraList)/2) +"\n"+
        separator+
        "i_j_Jij\n"+
        separator)
for para in paraList:
    if para[0] < para[1]:
        f.write("{0}\t{1}\t{2}\n".format(para[0],para[1],para[2]))
f.close()

### gutzwilleridx.def ###
NGutzwiller = 1
f = open(fileName[9],'w')
f.write(separator+
        "NGutzwillerIdx\t"+ str(NGutzwiller) +"\n"+
        separator+
        "i_GutzwillerIdx\n"+
        separator)
for i in range(Nsite):
    f.write("{0}\t{1}\n".format(i,0))
f.write("{0}\t{1}\n".format(0,1)) # optimized
f.close()

### jastrow.def ###
jastrow = {}
idx=-1
#idx=0
for i in range(Nsite):
  #for j in range(Nsite):
    #if i==j:
    #  continue
    #r = distance(i,j)
    #if r in jastrow:
    #  continue;
    #else:
    #  jastrow[r] = idx
    #  idx += 1
    #
    dr = indexToPosition(i)
    jastrow[dr] = idx
    idx += 1
    #
    #dr = indexToPosition(i)
    #dr_rev = (Lx-dr[0], Ly-dr[1])
    #dr_rev2 = (Ly-dr[1], Lx-dr[0])
    #if dr_rev in jastrow:
    #    jastrow[dr] = jastrow[dr_rev]
    #elif dr_rev2 in jastrow:
    #    jastrow[dr] = jastrow[dr_rev2]
    #else:
    #    if dr[0]==0 and dr[1] > Ly/2: 
    #        jastrow[dr] = jastrow[(dr[0], Ly-dr[1])]
    #        continue
    #    elif dr[1]==0 and dr[0] > Lx/2:
    #        jastrow[dr] = jastrow[(Lx-dr[0], dr[1])]
    #        continue
    #    else:
    #        jastrow[dr] = idx
    #        idx += 1
NJastrow = idx # Ns/2+1 (Lx,Ly:even), Ns/2 (Lx:even, Ly:odd)
f = open(fileName[10],'w')
f.write(separator+
        "NJastrowIdx\t"+ str(NJastrow) +"\n"+
        separator+
        "i_j_JastrowIdx\n"+
        separator)
for i in range(Nsite):
    for j in range(Nsite):
        if i==j:
            continue
        dr=direction(i,j)
        idx = jastrow[dr]
        #r=distance(i,j)
        #idx = jastrow[r]
        f.write("{0}\t{1}\t{2}\n".format(i,j,idx))
for i in range(NJastrow):
    if i==(Lx/2)-1:
        f.write("{0}\t{1}\n".format(i,1)) #optimized
        # f.write("{0}\t{1}\n".format(i,0)) #fixed
    else:
        f.write("{0}\t{1}\n".format(i,1)) #optimized
f.close()

### doublonholon2siteidx.def ###
NDoublonHolon2 = 0
#NDoublonHolon2 = 4
f = open(fileName[11],'w')
f.write(separator+
        "NDoublonHolon2SiteIdx\t"+ str(NDoublonHolon2) +"\n"+
        separator+
        "i_xi_xi_DoublonHolon2siteIdx\n"+
        separator)
t=0 # nn_x
for i in range(Nsite):
    j1 = neighborIndex(i,[1,0])
    j2 = neighborIndex(i,[-1,0])
    f.write("{0}\t{1}\t{2}\t{3}\n".format(i,j1,j2,t))
t=1 # nn_y
for i in range(Nsite):
    j1 = neighborIndex(i,[0,1])
    j2 = neighborIndex(i,[0,-1])
    f.write("{0}\t{1}\t{2}\t{3}\n".format(i,j1,j2,t))
t=2 # nn_xy
for i in range(Nsite):
    j1 = neighborIndex(i,[1,1])
    j2 = neighborIndex(i,[-1,-1])
    f.write("{0}\t{1}\t{2}\t{3}\n".format(i,j1,j2,t))
t=3 # nn_-xy
for i in range(Nsite):
    j1 = neighborIndex(i,[1,-1])
    j2 = neighborIndex(i,[-1,1])
    f.write("{0}\t{1}\t{2}\t{3}\n".format(i,j1,j2,t))
for i in range(6*NDoublonHolon2):
    if i<4*NDoublonHolon2:
        f.write("{0}\t{1}\n".format(i,1)) # optimized
    else:
        f.write("{0}\t{1}\n".format(i,1)) # optimized
        # f.write("{0}\t{1}\n".format(i,0)) # fixed
f.close()

### doublonholon4siteidx.def ###
#NDoublonHolon4 = 0
NDoublonHolon4 = 2
f = open(fileName[12],'w')
f.write(separator+
        "NDoublonHolon4SiteIdx\t"+ str(NDoublonHolon4) +"\n"+
        separator+
        "i_xi_xi_xi_xi_DoublonHolon4siteIdx\n"+
        separator)
t=0 # nn_+
for i in range(Nsite):
    j1 = neighborIndex(i,[1,0])
    j2 = neighborIndex(i,[-1,0])
    j3 = neighborIndex(i,[0,1])
    j4 = neighborIndex(i,[0,-1])
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j1,j2,j3,j4,t))
t=1 # nn_x
for i in range(Nsite):
    j1 = neighborIndex(i,[1,1])
    j2 = neighborIndex(i,[-1,1])
    j3 = neighborIndex(i,[1,-1])
    j4 = neighborIndex(i,[-1,-1])
    f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j1,j2,j3,j4,t))
for i in range(10*NDoublonHolon4):
    if i<4*NDoublonHolon2:
        f.write("{0}\t{1}\n".format(i,1)) # optimized
    else:
        f.write("{0}\t{1}\n".format(i,1)) # optimized
        # f.write("{0}\t{1}\n".format(i,0)) # fixed
f.close()

### orbitalidx.def ###
NOrbital = Nsub*Nsite
f = open(fileName[13],'w')
f.write(separator+
        "NOrbitalIdx\t"+ str(NOrbital) +"\n"+
        separator+
        "i_j_OrbitalIdx_OrbitalSgn\n"+
        separator)
idx = 0
orbital=[[-1 for dj in range(Nsite)] for isub in range(Nsub)]
for isub in range(Nsub):
    for dj in range(Nsite):
        orbital[isub][dj] = idx
        idx += 1
for i in range(Nsite):
    for j in range(Nsite):
        #sgn = 1
        #if (i > j):
        #  sgn = -1
        isub = subIndex(i)
        dj = positionToIndex(direction(i,j))
        sgn = sgnAP(i,direction(i,j))
        #f.write("{0}\t{1}\t{2}\t{3}\n".format(i,j,orbital[isub][dj],1))
        f.write("{0}\t{1}\t{2}\t{3}\n".format(i,j,orbital[isub][dj],sgn))
for i in range(idx):
    if i==0:
        f.write("{0}\t{1}\n".format(i,1))
    else:
        f.write("{0}\t{1}\n".format(i,1))
f.close()

### qptransidx.def ###
f = open(fileName[14],'w')
f.write(separator+
        "NQPTrans\t"+ str(Nsub) +"\n"+
        separator+
        "TrIdx_TrWeight_and_TrIdx_i_xi_TrSgn_i_xi\n"+
        separator)
for i in range(Nsub):
    f.write("{0}\t{1}\n".format(i, 1.0))
for dy in range(Sy):
    for dx in range(Sx):
        idx = dx+dy*Sx
        for i in range(Nsite):
            j = neighborIndex(i,[dx,dy])
            #dr = direction(i,j)
            #ri=indexToPosition(i)
            #rj=indexToPosition(j)
            #sgn = 1
            #if ((dx == Sx-1) and (i>j)):
            #  sgn = -1
            sgn = sgnAP(i,[dx,dy])
            f.write("{0}\t{1}\t{2}\t{3}\n".format(idx, i, j, sgn))
            #f.write("{0}\t{1}\t{2}\{3}\n".format(idx, i, j, 1))
f.close()

### cisajs.def ###
f = open(fileName[15],'w')
f.write(separator+
        "NCisAjs\t"+ str(Nsite*Nsite*2) +"\n"+
        separator+
        "idx_i_j_s\n"+
        separator)
for i in range(Nsite):
    for j in range(Nsite):
        f.write("{0}\t{1}\t{2}\t{3}\n".format(locgrnIdx(i,j,0), i, j, 0))
for i in range(Nsite):
    for j in range(Nsite):
        f.write("{0}\t{1}\t{2}\t{3}\n".format(locgrnIdx(i,j,1), i, j, 1))
f.close()

### cisajscktalt.def ###
f = open(fileName[16],'w')
f.write(separator+
        "NCisAjsCktAlt\t"+ str(Nsite*Nsite*(4+2+4*4)) +"\n"+
        #"NCisAjsCktAlt\t"+ str(0) +"\n"+
        separator+
        "idxJIS_idxKLT_i_j_s_k_l_t\n"+
        separator)
for i in range(Nsite):
    for j in range(Nsite):
        for s in range(2):
            for t in range(2):
                idx1 = locgrnIdx(i,i,s)
                idx2 = locgrnIdx(j,j,t)
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,i,s,j,j,t))
for i in range(Nsite):
    for j in range(Nsite):
        idx1 = locgrnIdx(i,j,0)
        idx2 = locgrnIdx(i,j,1)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,0,j,i,1))
        idx1 = locgrnIdx(i,j,1)
        idx2 = locgrnIdx(i,j,0)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,1,j,i,0))
for i in range(Nsite):
  i1 = neighborIndex(i,[1,0])
  i2 = neighborIndex(i,[-1,0])
  i3 = neighborIndex(i,[0,1])
  i4 = neighborIndex(i,[0,-1])
  for j in range(Nsite):
    j1 = neighborIndex(j,[1,0])
    j2 = neighborIndex(j,[-1,0])
    j3 = neighborIndex(j,[0,1])
    j4 = neighborIndex(j,[0,-1])
    for s in range(1):
      for t in range(1):
        t = 1-t
        idx1 = locgrnIdx(i,j,s)
        idx2 = locgrnIdx(j1,i1,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i1,j1,t))
        idx2 = locgrnIdx(j2,i1,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i1,j2,t))
        idx2 = locgrnIdx(j3,i1,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i1,j3,t))
        idx2 = locgrnIdx(j4,i1,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i1,j4,t))
        idx2 = locgrnIdx(j1,i2,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i2,j1,t))
        idx2 = locgrnIdx(j2,i2,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i2,j2,t))
        idx2 = locgrnIdx(j3,i2,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i2,j3,t))
        idx2 = locgrnIdx(j4,i2,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i2,j4,t))
        idx2 = locgrnIdx(j1,i3,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i3,j1,t))
        idx2 = locgrnIdx(j2,i3,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i3,j2,t))
        idx2 = locgrnIdx(j3,i3,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i3,j3,t))
        idx2 = locgrnIdx(j4,i3,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i3,j4,t))
        idx2 = locgrnIdx(j1,i4,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i4,j1,t))
        idx2 = locgrnIdx(j2,i4,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i4,j2,t))
        idx2 = locgrnIdx(j3,i4,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i4,j3,t))
        idx2 = locgrnIdx(j4,i4,t)
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\n".format(idx1,idx2,i,j,s,i4,j4,t))
f.close()

### cisajscktaltdc.def ###
f = open(fileName[17],'w')
f.write(separator+
        "NCisAjsCktAltDC\t"+ str(0) +"\n"+
        #"NCisAjsCktAltDC\t"+ str(Nsite*Nsite*(4+2+4*4)) +"\n"+
        #"NCisAjsCktAltDC\t"+ str(Nsite*Nsite*(4+2)+Nsite*16*2-16) +"\n"+
        #"NCisAjsCktAltDC\t"+ str(Nsite*16*2-16) +"\n"+
        separator+
        "i_j_s_k_l_t\n"+
        separator)
for i in range(Nsite):
    for j in range(Nsite):
        for s in range(2):
            for t in range(2):
                f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,i,s,j,j,t))
for i in range(Nsite):
    for j in range(Nsite):
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,0,j,i,1))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,1,j,i,0))
for i in range(Nsite):
  i1 = neighborIndex(i,[1,0])
  i2 = neighborIndex(i,[-1,0])
  i3 = neighborIndex(i,[0,1])
  i4 = neighborIndex(i,[0,-1])
  for j in range(Nsite):
    j1 = neighborIndex(j,[1,0])
    j2 = neighborIndex(j,[-1,0])
    j3 = neighborIndex(j,[0,1])
    j4 = neighborIndex(j,[0,-1])
    #if( (i!=0) and (j!=0) ):
    #  continue
    for s in range(1):
      for t in range(1):
        t = 1-t
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i1,j1,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i1,j2,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i1,j3,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i1,j4,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i2,j1,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i2,j2,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i2,j3,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i2,j4,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i3,j1,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i3,j2,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i3,j3,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i3,j4,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i4,j1,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i4,j2,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i4,j3,t))
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(i,j,s,i4,j4,t))
f.close()

### interall.def ###
f = open(fileName[18],'w')
f.write(separator+
        "NInterAll\t0\n"+
        separator+
        "i_j_s_k_l_t_Jijsklt\n"+
        separator)
f.close()


### quench.def ###
f = open(fileName[19],'w')
f.write(separator+
        "Nquench\t"+str(Nsite)+"\n"+
        "Time\t{0}\t{1}\n".format(TInit, TEnd)+
        "Tquench\t{0}\t{1}\t{2}\n".format(Time1, Time2, Time3)+
        separator)
for i in range(Nsite):
    f.write("{0}\t{1}\t{2}\n".format(i,0.0,lambdaScale*U/eScale))
f.close()


### Physical Calculation  ###
NStoreO = 0

fileName = []
fileName.append('xnamelist_aft.def')
fileName.append(pre+'modpara_aft.def')
fileName.append(pre+'locspn.def')
fileName.append(pre+'transfer.def')
fileName.append(pre+'coulombintra.def')
fileName.append(pre+'coulombinter.def')
fileName.append(pre+'hund.def')
fileName.append(pre+'pairhop.def')
fileName.append(pre+'exchange.def')
fileName.append(pre+'gutzwilleridx.def')
fileName.append(pre+'jastrowidx.def')
fileName.append(pre+'doublonholon2siteidx.def')
fileName.append(pre+'doublonholon4siteidx.def')
fileName.append(pre+'orbitalidx.def')
fileName.append(pre+'qptransidx.def')
fileName.append(pre+'cisajs.def')
fileName.append(pre+'cisajscktalt.def')
fileName.append(pre+'cisajscktaltdc.def')
fileName.append(pre+'interall.def')
fileName.append(pre+'quench.def')

### xnamelist_aft.def ###
f = open(fileName[0],'w')
for x in fileName[1:]:
    f.write(x+"\n")
f.close()

### modpara_aft.def ###
f = open(fileName[1],'w')
f.write(
    separator+
    "Model_Parameters  0\n"+
    separator+
    "VMC_Cal_Parameters\n"+
    separator+
    "CDataFileHead  zvo_aft\n"+
    "CParaFileHead  zqp_aft\n"+
    separator+
    "NVMCCalMode    1\n"+
    "NLanczosMode   0\n"+
    separator+
    "NDataIdxStart  0\n"+
    "NDataQtySmp    10\n"+
    separator+
    "Nsite          {0}\n".format(Nsite)+
    "Nelectron      {0}\n".format(Nelectron)+
    "NSPGaussLeg    1\n"+
    "NSPStot        0\n"+
    "NMPTrans       -1\n"+
    "NSROptItrStep  4000\n"+
    "NSROptItrSmp   500\n"+
    "NSROptFixSmp   1\n"+
    "DSROptRedCut   5e-4\n"+
    "DSROptStaDel   2e-2\n"+
    "DSROptStepDt   2e-2\n"+
    "NVMCWarmUp     10\n"+
    "NVMCIniterval  1\n"+
    "NVMCSample     2000\n"+
    "NExUpdatePath  0\n"+
    "RndSeed        {0}\n".format(seed)+
    "NSplitSize     {0}\n".format(NSplitSize)+
    "NStoreO        {0}\n".format(NStoreO)+
    "NStoreM        {0}\n".format(NStoreM))
f.close()

