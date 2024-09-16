import cmath
import math
import sys

import numpy as np
import toml

import MakeInput
import qlms_lattice


def main():
    #[s] tolm load
    input_file                   = sys.argv[1] 
    list_org,list_sub,input_dict = MakeInput.read_toml(input_file)

    ini_cnt,max_cnt,calcmode = ReadModpara(input_dict)
    All_N    = list_org[0]*list_org[1]*list_org[2]*list_org[3]
    All_site = list_org[0]*list_org[1]*list_org[2]
    orb_num  = list_org[3]
    dir_name = input_dict["mVMC_aft"]["directory"]

    tot_Ene    = np.zeros([max_cnt], dtype=np.float64)
    tot_occ    = np.zeros([max_cnt,orb_num], dtype=np.float64)
    tot_AF     = np.zeros([max_cnt,orb_num], dtype=np.float64)
    for i_smp in range(ini_cnt,ini_cnt+max_cnt):
        file_name = "%s/zvo_cisajs_00%d.dat" % (dir_name,i_smp)
        occ,AF    = ReadG1(file_name,list_org,i_smp)
        for orb_i in range(orb_num):
            tot_occ[i_smp][orb_i] = occ[orb_i] 
            tot_AF[i_smp][orb_i]  = AF[orb_i] 
        #
        file_name = "%s/zvo_out_00%d.dat" % (dir_name,i_smp)
        tot_Ene[i_smp] = ReadEne(file_name,list_org,i_smp)

    Ave_Ene = np.mean(tot_Ene,axis=0)
    Err_Ene = np.std(tot_Ene,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    # 
    Ave_occ = np.mean(tot_occ,axis=0)
    Err_occ = np.std(tot_occ,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    Ave_AF  = np.mean(tot_AF,axis=0)
    Err_AF  = np.std(tot_AF,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    
    with open("Ene.dat", 'w') as f:
        print("%s         " % ("# Ene err_Ene Ene/(All_site) err_Ene/(All_site)"),file=f)
        print("%f %f %f %f" % (Ave_Ene,Err_Ene,Ave_Ene/(All_site),Err_Ene/(All_site)),file=f)
    with open("occ.dat", 'w') as f:
        print("%s         " % ("# occ err_occ AF err_AF"),file=f)
        for orb_i in range(orb_num):
            print("%f %f %f %f" % (Ave_occ[orb_i],Err_occ[orb_i],Ave_AF[orb_i],Err_AF[orb_i]),end="",file=f)
        print(" " ,file=f)

def ReadModpara(input_dict):
    file_name = input_dict["mVMC_aft"]["modpara"]
    with open(file_name) as f:
        data      = f.read()
        data      = data.split("\n")
        for i in range(0,len(data)):
            if data[i]: # if data[i] is not empty
                tmp = data[i].split()
                if tmp[0] == "NDataIdxStart":
                    ini_cnt = int(tmp[1])
                if tmp[0] == "NDataQtySmp":
                    max_cnt = int(tmp[1])
                if tmp[0] == "NVMCCalMode":
                    calcmode = int(tmp[1])
    return ini_cnt,max_cnt,calcmode
 

def ReadEne(file_name,list_org,i_smp):
    with open(file_name) as f:
        data      = f.read()
        data      = data.split("\n")
        #print(len(data))
    tmp     = data[0].split()
    tmp_Ene = tmp[0]
    return tmp_Ene
 

def ReadG1(file_name,list_org,i_smp):
    Lx       = list_org[0]
    Ly       = list_org[1]
    Lz       = list_org[2]
    orb_num  = list_org[3]
    #file_name = "output/zvo_aft_cisajs_001.dat" 
    with open(file_name) as f:
        data      = f.read()
        data      = data.split("\n")
        #print(len(data))
    #[s] count not empty elements
    occ = np.zeros([orb_num], dtype=np.float64)
    AF  = np.zeros([orb_num], dtype=np.float64)
    for i in range(0,len(data)):
        if data[i]: # if data[i] is not empty
            tmp = data[i].split()
            if tmp[0] == tmp[2]:
                all_i     = int(tmp[0])
                list_site = qlms_lattice.get_site(all_i,list_org)
                x_i       = list_site[0]
                y_i       = list_site[1]
                orb_i     = list_site[3]
                sgn       = math.cos(math.pi*x_i+math.pi*y_i)   
                if tmp[1] == tmp[3] and int(tmp[1]) == 0: 
                    occ[orb_i]   += float(tmp[4])   
                    AF[orb_i]    += sgn*float(tmp[4]) 
                if tmp[1] == tmp[3] and int(tmp[1]) == 1: 
                    occ[orb_i]   += float(tmp[4])   
                    AF[orb_i]    += -sgn*float(tmp[4]) 
    #[e] count not empty elements
    occ   = occ/(Lx*Ly*Lz)
    AF    = AF/(Lx*Ly*Lz)
    return occ,AF


if __name__ == "__main__":
    main()
