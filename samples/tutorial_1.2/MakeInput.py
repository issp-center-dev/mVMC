import numpy as np
import math
import cmath
import toml
import sys

def main():
    #[s] tolm load
    input_file                   = sys.argv[1] 
    list_org,list_sub,input_dict = read_toml(input_file)
    #[e] tolm load
    #[s] output StdFace files
    OutputStdFace(list_org,list_sub,input_dict)
    #[s] output StdFace files
    #[s] output StdFace files
    OutputGreen(list_org)
    #[s] output StdFace files

def read_toml(input_file):
    input_dict  = toml.load(input_file)
    #[e] tolm load
    #[s]define constants
    Lx            = int(input_dict["lattice"]["Lx"])
    Ly            = int(input_dict["lattice"]["Ly"])
    Lz            = int(input_dict["lattice"]["Lz"])
    orb_num       = int(input_dict["lattice"]["orb_num"])
    sub_x         = int(input_dict["mVMC"]["sub_x"])
    sub_y         = int(input_dict["mVMC"]["sub_y"])
    sub_z         = int(input_dict["mVMC"]["sub_z"])
    #[e]define constants
    All_N = Lx*Ly*Lz*orb_num
    print('Lx      = ',Lx)
    print('Ly      = ',Ly)
    print('Ly      = ',Lz)
    print('orb_num = ',orb_num)
    print('sub_x   = ',sub_x)
    print('sub_y   = ',sub_y)
    print('sub_z   = ',sub_z)

    #[s] initialize
    list_org   = [Lx,Ly,Lz,orb_num]
    list_sub   = [sub_x,sub_y,sub_z]
    #[e] initialize

    return list_org,list_sub,input_dict

def CalcDim(list_org):
    Lx      = list_org[0]
    Ly      = list_org[1]
    Lz      = list_org[2]
    orb_num = list_org[3]

    #[s] calc dim_type
    if Ly == 1 and Lz == 1: 
        dim_type = 1
    elif Lz == 1:  
        dim_type = 2
    elif Lx>1 and Ly > 1 and Lz>1:  
        dim_type = 3
    else :
        print(" Possible error in Lx, Ly, Lz")
        print(" This script supports chain,square, and cubic in the following way: ")
        print(" Lx>1 Ly=Lz=1   -> chain")
        print(" Lx>1 Ly>1 Lz=1 -> square")
        print(" Lx>1 Ly>1 Lz>1 -> cubic")
        print(" Input Lx=%d Ly=%d  Lz=%d do not meet above conditions." %(Lx,Ly,Lz))
    #[e] calc dim_type
    return dim_type

def OutputStdFace(list_org,list_sub,input_dict):
    Lx          = list_org[0]
    Ly          = list_org[1]
    Lz          = list_org[2]
    orb_num     = list_org[3]
    Ncond       = Lx*Ly*Lz*orb_num
    sub_x       = list_sub[0]
    sub_y       = list_sub[1]
    sub_z       = list_sub[2]
    model_type  = input_dict["lattice"]["model_type"]

    dim_type = CalcDim(list_org)

    #[s] make stan_com
    stan_com = []
    if dim_type == 1:
        stan_com.append("L             = %d       "%(Lx))
        stan_com.append("Lsub          = %d       "%(sub_x))
        stan_com.append("lattice       = \"chain\"")
    elif dim_type == 2:
        stan_com.append("W             = %d       "%(Lx))
        stan_com.append("Wsub          = %d       "%(sub_x))
        stan_com.append("L             = %d       "%(Ly))
        stan_com.append("Lsub          = %d       "%(sub_y))
        stan_com.append("lattice       = \"square\"")
    elif dim_type == 3:
        stan_com.append("W             = %d       "%(Lx))
        stan_com.append("Wsub          = %d       "%(sub_x))
        stan_com.append("L             = %d       "%(Ly))
        stan_com.append("Lsub          = %d       "%(sub_y))
        stan_com.append("Height        = %d       "%(Lz))
        stan_com.append("Hsub          = %d       "%(sub_z))
        stan_com.append("lattice       = \"cubic\"")

    if model_type == "Spin":
        stan_com.append("model         = \"%s\"   "%(model_type))
        stan_com.append("J             = 1.0      ")
    elif model_type == "Hubbard":
        stan_com.append("model         = \"%s\"   "%(model_type))
        stan_com.append("t             = 1.0      ")
        stan_com.append("U             = 4.0      ")
        stan_com.append("ncond         = %d       " %(Ncond))
    elif model_type == "Kondo":
        stan_com.append("model         = \"%s\"   "%(model_type))
        stan_com.append("t             = 1.0      ")
        stan_com.append("J             = 1.0      ")
        stan_com.append("ncond         = %d       " %(int(Ncond/2)))
    else: 
        print("This scropt only support Spin, Hubbard, Kondo")

    stan_com.append("2Sz           = 0        ")
    stan_com.append("NVMCSample    = 200      ")
    stan_com.append("NSROptItrStep = 600      ")
    stan_com.append("NMPTrans      = 1        ")
    stan_com.append("NSPStot       = 0        ")
    #[e] make stan_com

    with open("stan_opt.in", 'w') as f:
        for cnt_std in stan_com:
            print(cnt_std,file=f)

    with open("stan_aft.in", 'w') as f:
        for cnt_std in stan_com:
            print(cnt_std,file=f)
        print("NVMCCalMode   = 1        ",file=f)
        print("NDataIdxStart = 0        ",file=f)
        print("NDataQtySmp   = 10        ",file=f)

 

def OutputGreen(list_org): 
    All_N = list_org[0]*list_org[1]*list_org[2]*list_org[3]

    with open("green1", 'w') as f:
        print("==================", file=f)
        print("onebody %d "%(2*All_N**2), file=f)
        print("==================", file=f)
        print("==================", file=f)
        print("==================", file=f)
        for all_i in range(0,All_N):
            for all_j in range(0,All_N):
                print(" %d %d %d %d "% (all_i,0,all_j,0), file=f)
                print(" %d %d %d %d "% (all_i,1,all_j,1), file=f)

    with open("green2", 'w') as f:
        print("==================", file=f)
        print("twobody %d "%(6*All_N**2), file=f)
        print("==================", file=f)
        print("==================", file=f)
        print("==================", file=f)
        for all_i in range(0,All_N):
            for all_j in range(0,All_N):
                print(" %d %d %d %d %d %d %d %d"% (all_i,0,all_i,0,all_j,0,all_j,0), file=f)
                print(" %d %d %d %d %d %d %d %d"% (all_i,0,all_i,0,all_j,1,all_j,1), file=f)
                print(" %d %d %d %d %d %d %d %d"% (all_i,1,all_i,1,all_j,0,all_j,0), file=f)
                print(" %d %d %d %d %d %d %d %d"% (all_i,1,all_i,1,all_j,1,all_j,1), file=f)
                #
                print(" %d %d %d %d %d %d %d %d"% (all_i,0,all_j,0,all_j,1,all_i,1), file=f)
                print(" %d %d %d %d %d %d %d %d"% (all_i,1,all_j,1,all_j,0,all_i,0), file=f)
 




if __name__ == "__main__":
    main()
