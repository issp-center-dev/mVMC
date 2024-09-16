import numpy as np
import math
import cmath
import toml
import sys

import MakeInput

def main():
    #[s] tolm load
    #[s] tolm load
    input_file                   = sys.argv[1] 
    list_org,list_sub,input_dict = MakeInput.read_toml(input_file)
    #[e] tolm load

    trans_1x,trans_1y,trans_2x,trans_2y,trans_3x,trans_3y = MakeFormFactor()
    #[s] output green files
    Ini_site = MakeIniSite(list_org,list_sub,"full")
    Output_SC(list_org,Ini_site,trans_1x,trans_1y,"swave","SC_1swave")
    #Output_SC(list_org,Ini_site,trans_2x,trans_2y,"diagonal","SC_2s2d.def")  for 2s-wave and 2d-wave
    #[e] output green files
 
def MakeIniSite(list_org,list_sub,type_ini):
    Lx      = list_org[0]
    Ly      = list_org[1]
    Lz      = list_org[2]
    if type_ini == "full":
        site_max = Lx*Ly*Lz
        Ini_site  = np.zeros(site_max, dtype=np.int64)
        for cnt_i in range(site_max):
            Ini_site[cnt_i] = cnt_i
    elif type_ini == "reduce":
        sub_Lx    = list_sub[0]
        sub_Ly    = list_sub[1]
        site_max  = sub_Lx*sub_Ly
        Ini_site  = np.zeros(site_max, dtype=np.int64)
        for cnt_i in range(site_max):
            tmp_x_i = cnt_i%sub_Lx
            tmp_y_i = int((cnt_i-tmp_x_i)/sub_Lx)
            all_i   = tmp_x_i+tmp_y_i*Lx
            Ini_site[cnt_i] = all_i
    return Ini_site
        

def MakeFormFactor():
    trans_1x  = np.zeros([1], dtype=np.int64)
    trans_1y  = np.zeros([1], dtype=np.int64)
    # 
    trans_2x  = np.zeros([4], dtype=np.int64)
    trans_2y  = np.zeros([4], dtype=np.int64)
    #
    trans_3x  = np.zeros([4], dtype=np.int64)
    trans_3y  = np.zeros([4], dtype=np.int64)
    #
    trans_1x[0] = 0
    trans_1y[0] = 0
    #
    trans_2x[0] = 1
    trans_2y[0] = 0
    #
    trans_2x[1] = -1
    trans_2y[1] = 0
    # 
    trans_2x[2] = 0
    trans_2y[2] = 1
    #
    trans_2x[3] = 0
    trans_2y[3] = -1
    #
    trans_3x[0] = 1
    trans_3y[0] = 1
    #
    trans_3x[1] = -1
    trans_3y[1] = -1
    # 
    trans_3x[2] = 1
    trans_3y[2] = -1
    #
    trans_3x[3] = -1
    trans_3y[3] = 1

    return trans_1x,trans_1y,trans_2x,trans_2y,trans_3x,trans_3y
 
   
def Output_SC(list_org,Ini_site,trans_x,trans_y,type_of_orbital,name_file):
    Lx      = list_org[0]
    Ly      = list_org[1]
    Lz      = list_org[2]
    orb_num  = list_org[3]
    ini_site_max = len(Ini_site)
    site_max = Lx*Ly*Lz
    num_trans = len(trans_x)
    if type_of_orbital == "diagonal":
        tot_orb = int(orb_num)
        cnt_max = 8*(site_max*ini_site_max)*(num_trans**2)*tot_orb
    elif type_of_orbital == "swave":
        tot_orb = int(orb_num)
        cnt_max = 2*(site_max*ini_site_max)*tot_orb*tot_orb
    else:
        print("fatal error type_or_orbital should be diagonal or swave")
    #print(' cnt_V = ',cnt_V)
    with open("%s" % (name_file), 'w') as f:
        print("====== ", file=f)
        print("%s %d " % ('N', cnt_max), file=f)
        print("====== ", file=f)
        print("====== ", file=f)
        print("====== ", file=f)
        if type_of_orbital == "swave":
            for tmp_site_i in range(ini_site_max):
                site_i = Ini_site[tmp_site_i]
                x_i = site_i%Lx
                y_i = int((site_i-x_i)/Lx)
                for orb_i in range(orb_num):
                    all_i     = orb_i + site_i*orb_num
                    for all_j in range(site_max*orb_num):
                        print(" %8d %8d %8d %8d %8d %8d %8d %8d " % (all_i,0,all_j,0,all_i,1,all_j,1), file=f)
                        print(" %8d %8d %8d %8d %8d %8d %8d %8d " % (all_j,0,all_i,0,all_j,1,all_i,1), file=f)
 
        if type_of_orbital == "diagonal":
            for tmp_site_i in range(ini_site_max):
                site_i = Ini_site[tmp_site_i]
                x_i = site_i%Lx
                y_i = int((site_i-x_i)/Lx)
                for site_j in range(site_max):
                    x_j = site_j%Lx
                    y_j = int((site_j-x_j)/Lx)
                    for cnt_i in range(len(trans_x)):
                        til_site_i = ((x_i+trans_x[cnt_i]+Lx)%Lx)+((y_i+trans_y[cnt_i]+Ly)%Ly)*Lx
                        for cnt_j in range(len(trans_x)):
                            til_site_j = ((x_j+trans_x[cnt_j]+Lx)%Lx)+((y_j+trans_y[cnt_j]+Ly)%Ly)*Lx
                            for orb_i in range(orb_num):
                                all_i     = orb_i + site_i*orb_num
                                til_all_i = orb_i + til_site_i*orb_num
                                all_j     = orb_i + site_j*orb_num
                                til_all_j = orb_i + til_site_j*orb_num
                                print(" %8d %8d %8d %8d %8d %8d %8d %8d " % (all_i,0,all_j,0,til_all_i,1,til_all_j,1), file=f)
                                print(" %8d %8d %8d %8d %8d %8d %8d %8d " % (all_j,0,all_i,0,til_all_j,1,til_all_i,1), file=f)
                                #
                                print(" %8d %8d %8d %8d %8d %8d %8d %8d " % (all_i,0,til_all_j,0,til_all_i,1,all_j,1), file=f)
                                print(" %8d %8d %8d %8d %8d %8d %8d %8d " % (til_all_j,0,all_i,0,all_j,1,til_all_i,1), file=f)
                                #
                                print(" %8d %8d %8d %8d %8d %8d %8d %8d " % (til_all_i,0,all_j,0,all_i,1,til_all_j,1), file=f)
                                print(" %8d %8d %8d %8d %8d %8d %8d %8d " % (all_j,0,til_all_i,0,til_all_j,1,all_i,1), file=f)
                                #
                                print(" %8d %8d %8d %8d %8d %8d %8d %8d " % (til_all_i,0,til_all_j,0,all_i,1,all_j,1), file=f)
                                print(" %8d %8d %8d %8d %8d %8d %8d %8d " % (til_all_j,0,til_all_i,0,all_j,1,all_i,1), file=f)


if __name__ == "__main__":
    main()
