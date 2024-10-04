import cmath
import math
import sys

import numpy as np
import toml

import MakeInput
import qlms_lattice
import VMClocal

def main():
    #[s] tolm load
    input_file                   = sys.argv[1] 
    list_org,list_sub,input_dict = MakeInput.read_toml(input_file)
    dir_name                     = input_dict["mVMC_aft"]["directory"]
    model_type                   = input_dict["lattice"]["model_type"]
    dim_type                     = MakeInput.CalcDim(list_org)
    #[e] tolm load

    ini_cnt,max_cnt,calcmode = VMClocal.ReadModpara(input_dict)

    #[s] read Green
    G1,G2_sz,G2_ex = ReadGreen(list_org,max_cnt,dir_name)
    #[e] read Green

    #[s] calc Sq,Sz,Nq,Nk
    all_Sq,all_Sz,all_Nq,all_Nk = CalcSq_tot(list_org,G1,G2_sz,G2_ex,dir_name,max_cnt)
    OutputSqSzNq(list_org,all_Sq,all_Sz,all_Nq,dim_type)
    if model_type == "Hubbard":
        OutputNk(list_org,all_Nk,dim_type) #Nk is output only Hubbard
    #[e] calc Sq,Sz,Nq,Nk
    OutputReal(list_org,G1,max_cnt,dim_type)
    OutputSij(list_org,G1,G2_sz,G2_ex,max_cnt)

def OutputSij(list_org,G1,G2_sz,G2_ex,max_cnt):
    Lx       = list_org[0]
    Ly       = list_org[1]
    Lz       = list_org[2]
    orb_num  = list_org[3]
    All_N    = Lx*Ly*Lz*orb_num 
    tot_S    = np.zeros((max_cnt,All_N,All_N),dtype=np.float64)
    Sxy      = np.zeros((max_cnt,All_N,All_N),dtype=np.float64)
    Sz       = np.zeros((max_cnt,All_N,All_N),dtype=np.float64)
    for num_bin in range(0,max_cnt):
        for all_i in range(0,All_N):
            for all_j in range(0,All_N):
                tmp_Sz   = G2_sz[num_bin][all_i][all_j][0][0]
                tmp_Sz  += -G2_sz[num_bin][all_i][all_j][1][0]
                tmp_Sz  += -G2_sz[num_bin][all_i][all_j][0][1]
                tmp_Sz  += G2_sz[num_bin][all_i][all_j][1][1]
                tmp_Sz   = 0.25*tmp_Sz
                tmp_Sxy  = -0.5*G2_ex[num_bin][all_i][all_j][0][1]
                tmp_Sxy += -0.5*G2_ex[num_bin][all_i][all_j][1][0]
                if all_i == all_j:
                    tmp_Sxy += 0.5*G1[num_bin][all_i][all_i][0]
                    tmp_Sxy += 0.5*G1[num_bin][all_i][all_i][1]
                tot_S[num_bin][all_i][all_j] = tmp_Sxy+tmp_Sz
                Sxy[num_bin][all_i][all_j]   = tmp_Sxy
                Sz[num_bin][all_i][all_j]    = tmp_Sz
    #
    ave_tot_S  = np.mean(tot_S,axis=0)
    err_tot_S  = np.std(tot_S,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    ave_Sxy    = np.mean(Sxy,axis=0)
    err_Sxy    = np.std(Sxy,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    ave_Sz     = np.mean(Sz,axis=0)
    err_Sz     = np.std(Sz,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    with open("Sij.dat", 'w') as f:
        print(" %s " % ("# i  j tot_S err_tot_S Sxy err_Sxy Sz err_Sz "), file=f)
        for all_i in range(All_N):
            for all_j in range(All_N):
                print("%d %d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f " % (all_i,all_j,\
                ave_tot_S[all_i][all_j],err_tot_S[all_i][all_j],\
                ave_Sxy[all_i][all_j],err_Sxy[all_i][all_j],\
                ave_Sz[all_i][all_j],err_Sz[all_i][all_j]), file=f)
            print("  ",  file=f)

def OutputReal(list_org,G1,max_cnt,dim_type):
    Lx        = list_org[0]
    Ly        = list_org[1]
    Lz        = list_org[2]
    orb_num   = list_org[3]
    All_N     = Lx*Ly*Lz*orb_num 
    charge    = np.zeros((max_cnt,All_N),dtype=np.float64)
    spin      = np.zeros((max_cnt,All_N),dtype=np.float64)
    for num_bin in range(0,max_cnt):
        for all_i in range(0,All_N):
            charge[num_bin][all_i] = G1[num_bin][all_i][all_i][0]+G1[num_bin][all_i][all_i][1]
            spin[num_bin][all_i]   = G1[num_bin][all_i][all_i][0]-G1[num_bin][all_i][all_i][1]
    #
    ave_charge  = np.mean(charge,axis=0)
    err_charge  = np.std(charge,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    ave_spin    = np.mean(spin,axis=0)
    err_spin    = np.std(spin,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    with open("Real.dat", 'w') as f:
        print(" %s " % ("# x  y charge err_charge spin err_spin"), file=f)
        for i_x in range(0,list_org[0]):
            if dim_type == 1:
                i_y   = 0
                all_i =  i_x + i_y*list_org[0]
                print("%d %d %12.8f %12.8f %12.8f %12.8f" % (i_x,i_y,ave_charge[all_i],err_charge[all_i],ave_spin[all_i],err_spin[all_i]), file=f)
            elif dim_type == 2:
                for i_y in range(0,list_org[1]):
                    all_i =  i_x + i_y*list_org[0]
                    print("%d %d %12.8f %12.8f %12.8f %12.8f" % (i_x,i_y,ave_charge[all_i],err_charge[all_i],ave_spin[all_i],err_spin[all_i]), file=f)
                print(" " , file=f)


def OutputSqSzNq(list_org,all_Sq,all_Sz,all_Nq,dim_type):
    max_cnt = all_Sq.shape[0]
    #[s] Sq,Sz,Nq
    ave_Sq  = np.mean(all_Sq,axis=0)
    err_Sq  = np.std(all_Sq,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    ave_Sz  = np.mean(all_Sz,axis=0)
    err_Sz  = np.std(all_Sz,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    ave_Nq  = np.mean(all_Nq,axis=0)
    err_Nq  = np.std(all_Nq,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    max_Sq  = 0.0
    max_Nq  = 0.0
    with open("SqNq.dat", 'w') as f:
        print(" %s " % ("# kx ky Sq err_Sq Sz err_Sz Nq err_Nq"), file=f)
        for kx in range(0,list_org[0]+1):
            if dim_type==1:
                ky = 0
                if ave_Sq[kx][ky] > max_Sq:
                    max_Sq     = ave_Sq[kx][ky] 
                    max_Sq_err = err_Sq[kx][ky] 
                    max_Sq_kx  = kx 
                    max_Sq_ky  = ky
                    #
                if ave_Nq[kx][ky] > max_Nq:
                    max_Nq     = ave_Nq[kx][ky] 
                    max_Nq_err = err_Nq[kx][ky] 
                    max_Nq_kx  = kx 
                    max_Nq_ky  = ky
                print("%d %d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f" % (kx,ky,ave_Sq[kx][ky],err_Sq[kx][ky],ave_Sz[kx][ky],err_Sz[kx][ky],ave_Nq[kx][ky],err_Nq[kx][ky]), file=f)
            elif dim_type==2:
                for ky in range(0,list_org[1]+1):
                    if ave_Sq[kx][ky] > max_Sq:
                        max_Sq     = ave_Sq[kx][ky] 
                        max_Sq_err = err_Sq[kx][ky] 
                        max_Sq_kx  = kx 
                        max_Sq_ky  = ky
                        #
                    if ave_Nq[kx][ky] > max_Nq:
                        max_Nq     = ave_Nq[kx][ky] 
                        max_Nq_err = err_Nq[kx][ky] 
                        max_Nq_kx  = kx 
                        max_Nq_ky  = ky
                    print("%d %d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f" % (kx,ky,ave_Sq[kx][ky],err_Sq[kx][ky],ave_Sz[kx][ky],err_Sz[kx][ky],ave_Nq[kx][ky],err_Nq[kx][ky]), file=f)
                print(" " , file=f)
            elif dim_type==3:
                print("not implemented yet")
    with open("MaxSq.dat", 'w') as f:
        print("%12.8f %12.8f  %d %d " % (max_Sq,max_Sq_err,max_Sq_kx,max_Sq_ky), file=f)
    with open("MaxNq.dat", 'w') as f:
        print("%12.8f %12.8f  %d %d " % (max_Nq,max_Nq_err,max_Nq_kx,max_Nq_ky), file=f)
    #[e] Sq,Nq

def OutputNk(list_org,all_Nk,dim_type):
    max_cnt = all_Nk.shape[0]
    #[s] Nk
    ave_Nk  = np.mean(all_Nk,axis=0)
    err_Nk  = np.std(all_Nk,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    with open("Nk.dat", 'w') as f:
        print(" %s " % ("# kx ky Nk err_Nk "), file=f)
        for kx in range(-int(list_org[0]/2),int(list_org[0]/2)+1):
            if kx <0:
                tmp_kx = kx+list_org[0]
            else:
                tmp_kx = kx
            if dim_type==1:
                ky     = 0
                tmp_ky = 0
                print("%d %d %12.8f %12.8f " % (kx,ky,ave_Nk[tmp_kx][tmp_ky],err_Nk[tmp_kx][tmp_ky]), file=f)
            elif dim_type==2:
                for ky in range(-list_org[1],list_org[1]+1):
                    if ky <0:
                        tmp_ky = ky+list_org[1]
                    else:
                        tmp_ky = ky
                    print("%d %d %12.8f %12.8f " % (kx,ky,ave_Nk[tmp_kx][tmp_ky],err_Nk[tmp_kx][tmp_ky]), file=f)
                print(" " , file=f)
            elif dim_type==3:
                print("Not implemented yet")
        #[e] Nk

def ReadGreen(list_org,max_cnt,dir_name):
    Lx        = list_org[0]
    Ly        = list_org[1]
    Lz        = list_org[2]
    orb_num   = list_org[3]
    All_N     = Lx*Ly*Lz*orb_num 
    #[s] allocate
    G1        = np.zeros((max_cnt,All_N,All_N,2),dtype=np.float64)
    G2_ex     = np.full((max_cnt,All_N,All_N,2,2),-100,dtype=np.float64)
    G2_sz     = np.full((max_cnt,All_N,All_N,2,2),-100,dtype=np.float64)
    #[e] allocate
    for num_bin in range(max_cnt):
        file_name ="%s/zvo_cisajs_00%d.dat" %(dir_name,num_bin)

        with open(file_name) as f:
            tmp_G1      = f.read()
            tmp_G1      = tmp_G1.split("\n")
        #[s] count not empty elements
        cnt = 0
        for i in range(0,len(tmp_G1)):
            if tmp_G1[i]: # if data[i] is not empty
                cnt += 1
        #print(cnt)
        cnt_max = cnt
        #[e] count not empty elements
        for cnt in range(0,cnt_max):
            tmp                   = tmp_G1[cnt].split()  
            all_i                 = int(tmp[0]) 
            all_j                 = int(tmp[2]) 
            spn                   = int(tmp[1])
            #print(tmp)
            G1[num_bin][all_i][all_j][spn] = float(tmp[4])
    
        file_name ="%s/zvo_cisajscktalt_00%d.dat" %(dir_name,num_bin)
        with open(file_name) as f:
            data      = f.read()
            data      = data.split("\n")
        #[s] count not empty elements
        cnt = 0
        for i in range(0,len(data)):
            if data[i]: # if data[i] is not empty
                cnt += 1
        #print(cnt)
        cnt_max = cnt
        #[e] count not empty elements
        for cnt in range(0,cnt_max):
            tmp                   = data[cnt].split()  
            all_i                 = int(tmp[0]) 
            all_j                 = int(tmp[2]) 
            all_k                 = int(tmp[4]) 
            all_l                 = int(tmp[6]) 
            spn_0                 = int(tmp[1])
            spn_1                 = int(tmp[5])
            if all_i == all_j and all_k == all_l and all_i==all_k:
                G2_sz[num_bin][all_i][all_j][spn_0][spn_1] = float(tmp[8])
                G2_ex[num_bin][all_i][all_j][spn_0][spn_1] = float(tmp[8])
            elif all_i == all_j and all_k == all_l:
                G2_sz[num_bin][all_i][all_k][spn_0][spn_1] = float(tmp[8])
            elif all_i == all_l and all_j == all_k:
                G2_ex[num_bin][all_i][all_j][spn_0][spn_1] = float(tmp[8])
            else:
                print("fatal error in 2-body Green functions")
    return G1,G2_sz,G2_ex

def CalcSq_tot(list_org,G1,G2_sz,G2_ex,dir_name,max_cnt):
    Lx        = list_org[0]
    Ly        = list_org[1]
    Lz        = list_org[2]
    orb_num   = list_org[3]
    All_N     = Lx*Ly*Lz*orb_num 
    All_site  = Lx*Ly*Lz
    all_Sq    = np.full((max_cnt,list_org[0]+1,list_org[1]+1),-100,dtype=np.float64)
    all_Sz    = np.full((max_cnt,list_org[0]+1,list_org[1]+1),-100,dtype=np.float64)
    all_Nq    = np.full((max_cnt,list_org[0]+1,list_org[1]+1),-100,dtype=np.float64)
    all_Nk    = np.full((max_cnt,list_org[0]+1,list_org[1]+1),-100,dtype=np.float64)

    for num_bin in range(max_cnt):
        with open("%s/total_SqNq_%d.dat" % (dir_name,num_bin), 'w') as f:
            for kx in range(0,Lx+1):
                for ky in range(0,Ly+1):
                    tmp_Sq  = 0.0
                    tmp_Sz  = 0.0
                    tmp_Nq  = 0.0
                    tmp_Nk  = 0.0
                    Ncond   = 0
                    for all_i in range(0,All_N):
                        Ncond+=  G1[num_bin][all_i][all_i][0]
                        Ncond+=  G1[num_bin][all_i][all_i][1]
                        for all_j in range(0,All_N):
                            list_site       = qlms_lattice.get_site(all_i,list_org)
                            i_x             = list_site[0]
                            i_y             = list_site[1]
                            #
                            list_site       = qlms_lattice.get_site(all_j,list_org)
                            j_x             = list_site[0]
                            j_y             = list_site[1]
            
                            theta           = 2*math.pi*kx*(i_x-j_x)/Lx+2*math.pi*ky*(i_y-j_y)/Ly
                            #
                            tmp_uu  = G2_sz[num_bin][all_i][all_j][0][0]
                            tmp_ud  = G2_sz[num_bin][all_i][all_j][0][1]
                            tmp_du  = G2_sz[num_bin][all_i][all_j][1][0]
                            tmp_dd  = G2_sz[num_bin][all_i][all_j][1][1]

                            #
                            tmp_Nk += G1[num_bin][all_i][all_j][0]*math.cos(theta)
                            tmp_Nk += G1[num_bin][all_i][all_j][1]*math.cos(theta)
                            #
                            tmp_Nq += tmp_uu*math.cos(theta)
                            tmp_Nq += tmp_dd*math.cos(theta)
                            tmp_Sq += 0.25*tmp_uu*math.cos(theta)
                            tmp_Sq += 0.25*tmp_dd*math.cos(theta)
                            tmp_Sz += 0.25*tmp_uu*math.cos(theta)
                            tmp_Sz += 0.25*tmp_dd*math.cos(theta)
                            #
                            tmp_Nq += tmp_ud*math.cos(theta)
                            tmp_Nq += tmp_du*math.cos(theta)
                            tmp_Sq += -0.25*tmp_ud*math.cos(theta)
                            tmp_Sq += -0.25*tmp_du*math.cos(theta)
                            tmp_Sz += -0.25*tmp_ud*math.cos(theta)
                            tmp_Sz += -0.25*tmp_du*math.cos(theta)
                            #
                            tmp_Sq += -0.5*G2_ex[num_bin][all_i][all_j][0][1]*math.cos(theta)
                            tmp_Sq += -0.5*G2_ex[num_bin][all_i][all_j][1][0]*math.cos(theta)
                    tmp_Sq += Ncond/2
                    if kx%Lx == 0 and ky%Ly ==0:
                       tmp_Nq = tmp_Nq- Ncond**2
                    tmp_Sq  = tmp_Sq/(All_site)
                    tmp_Sz  = tmp_Sz/(All_site)
                    tmp_Nq  = tmp_Nq/(All_site)
                    tmp_Nk  = tmp_Nk/(All_site)
                    all_Sq[num_bin][kx][ky] = tmp_Sq
                    #print(num_bin,kx,ky)
                    all_Sz[num_bin][kx][ky] = tmp_Sz
                    all_Nq[num_bin][kx][ky] = tmp_Nq
                    all_Nk[num_bin][kx][ky] = tmp_Nk
                    #print(kx,ky,tmp_Sq,tmp_Sz,tmp_Nq)
                    print("%d %d %12.8f %12.8f %12.8f " % (kx,ky,tmp_Sq,tmp_Sz,tmp_Nq), file=f)
                print(" " , file=f)
                #print(" ")
    return all_Sq,all_Sz,all_Nq,all_Nk

if __name__ == "__main__":
    main()
