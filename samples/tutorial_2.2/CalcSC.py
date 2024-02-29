import numpy as np
import math
import cmath
import toml
import sys
import SCGreen
import MakeInput
import VMClocal

def main():

    #[s] tolm load
    input_file                   = sys.argv[1] 
    list_org,list_sub,input_dict = MakeInput.read_toml(input_file)
    ini_cnt,max_cnt,calcmode     = VMClocal.ReadModpara(input_dict)
    dir_name = input_dict["mVMC_aft"]["directory"]
    All_N    = list_org[0]*list_org[1]*list_org[2]*list_org[3]
    #[e] tolm load

    trans_1x,trans_1y,trans_2x,trans_2y,trans_3x,trans_3y = SCGreen.MakeFormFactor()
    print(trans_2x)
    print(trans_2y)
    #[s] output misc files
    Ini_site = SCGreen.MakeIniSite(list_org,list_sub,"full")
    print(Ini_site)

    G1 =  read_G1(All_N,max_cnt,dir_name)
    #
    G2_swave =  read_G2_swave(All_N,max_cnt,dir_name,G1,trans_1x,trans_1y)
    ave_SC_swave,err_SC_swave = ave_swave_G2(list_org,max_cnt,dir_name,G2_swave)
    Output_SC(list_org,Ini_site,ave_SC_swave,err_SC_swave,"Result_1swave.dat")
    #

def ave_swave_G2(list_org,max_cnt,dir_name,G2):
    Lx           = list_org[0]
    Ly           = list_org[1]
    Lz           = list_org[2]
    orb_num      = list_org[3]
    All_site     = Lx*Ly*Lz
 
    SC_swave = np.zeros((max_cnt,All_site),dtype=np.float64)
    for num_bin in range(0,max_cnt):
        for site_i in range(0,All_site):
            vec_x_i = site_i %Lx
            vec_y_i = int((site_i-vec_x_i)/Lx)
            all_i = site_i
            tmp   = 0.0
            for site_j in range(0,All_site):
                ini_x_j = site_j %Lx
                ini_y_j = int((site_j-ini_x_j)/Lx)
                ini_all_j = site_j
                x_j = (ini_x_j+vec_x_i+Lx)%Lx
                y_j = (ini_y_j+vec_y_i+Ly)%Ly
                all_j   = x_j+y_j*Lx
                #if all_i == 0 :
                #    print(num_bin,all_i,all_j, G2[num_bin][all_j][all_i][0][0], G2[num_bin][all_j][all_i][0][1])
                tmp    += G2[num_bin][ini_all_j][all_j][0][0]+G2[num_bin][ini_all_j][all_j][0][1]
            SC_swave[num_bin][site_i] = tmp/(4.0*All_site)
    print(SC_swave[0][0])
    ave_G2  = np.mean(SC_swave,axis=0)
    err_G2  = np.std(SC_swave,axis=0,ddof=1)/math.sqrt(1.0*max_cnt)
    #for all_i in range(0,All_N,2):
    #    for all_j in range(0,All_N,2):
    #        print(all_i,all_j,ave_G2[all_i][all_j])
    return ave_G2,err_G2
    #print(ave_G2)

def read_G2_swave(All_N,max_cnt,dir_name,G1,trans_x,trans_y):
    num_neighbor = len(trans_x)*len(trans_y)
    print("num_neighbor",num_neighbor)
    G2 = np.zeros((max_cnt,All_N,All_N,num_neighbor,2),dtype=np.float64)
    for num_bin in range(0,max_cnt):
        file_name ="{}".format(dir_name)+"/zvo_cisajscktalt_00"+"{0:1d}".format(num_bin)+".dat"

        with open(file_name) as f:
            data      = f.read()
            data      = data.split("\n")
            print(len(data))
        #[s] count not empty elements
        cnt = 0
        for i in range(0,len(data)):
            if data[i]: # if data[i] is not empty
                cnt += 1
        #print(cnt)
        cnt_max = cnt
        print("G2",cnt_max)
        #[e] count not empty elements
        for cnt_neighbor in range(num_neighbor):
            for cnt in range(0,cnt_max,2):
                tmp                   = data[cnt].split()  
                all_i                 = int(tmp[0]) 
                all_j                 = int(tmp[2]) 
                til_all_i             = int(tmp[4]) 
                til_all_j             = int(tmp[6]) 
                delta_0               = 0.0
                delta_1               = 0.0
                delta_2               = 0.0
                delta_3               = 0.0
                tmp              = data[cnt].split()  
                G2[num_bin][all_i][all_j][cnt_neighbor][0] = float(tmp[8])

                if all_i == all_j:
                    delta_0 = 1.0
                if til_all_i == til_all_j:
                    delta_1 = 1.0
                if all_i == til_all_j:
                    delta_2 = 1.0
                if til_all_i == all_j:
                    delta_3 = 1.0

                tmp_2     = data[cnt+1].split()  
                tmp_G     = delta_0*delta_1-delta_1*G1[num_bin][all_j][all_i][0]-delta_0*G1[num_bin][til_all_j][til_all_i][1]
                #print(all_i,all_j,til_all_i,til_all_j,1,tmp_G)
                G2[num_bin][all_i][all_j][cnt_neighbor][1] = float(tmp_2[8])+tmp_G
    return G2

def read_G1(All_N,max_cnt,dir_name):
    G1        = np.zeros((max_cnt,All_N,All_N,2),dtype=np.float64)
    for num_bin in range(0,max_cnt):
        file_name ="{}".format(dir_name)+"/zvo_cisajs_00"+"{0:1d}".format(num_bin)+".dat"
        with open(file_name) as f:
            tmp_G1      = f.read()
            tmp_G1      = tmp_G1.split("\n")
            print(len(tmp_G1))
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
    return G1
   
def Output_SC(list_org,Ini_site,ave_SC,err_SC,name_file):
    Lx           = list_org[0]
    Ly           = list_org[1]
    Lz           = list_org[2]
    orb_num      = list_org[3]
    ini_site_max = len(Ini_site)
    site_max     = Lx*Ly*Lz

    max_num_dis  = (int(Lx/2)+1)**2+(int(Ly/2)+1)**2
    ave_SC_dis   = np.zeros(max_num_dis,dtype=np.float64)
    err_SC_dis   = np.zeros(max_num_dis,dtype=np.float64)
    for site_j in range(site_max):
        x_j = site_j%Lx
        y_j = int((site_j-x_j)/Lx)
        diff_x = abs(x_j)
        if diff_x > Lx/2:
            diff_x = diff_x - Lx
        diff_y = abs(y_j)
        if diff_y > Ly/2:
            diff_y = diff_y-Ly
        tmp_dis  = diff_x**2+diff_y**2
        #print(x_i,y_i,diff_x,x_j,y_j,diff_y,tmp_dis,max_cnt)
        if abs(ave_SC[site_j]) > abs(ave_SC_dis[tmp_dis]):
            ave_SC_dis[tmp_dis] = ave_SC[site_j]
            err_SC_dis[tmp_dis] = err_SC[site_j]

    with open("%s" % (name_file) , 'w') as f:
        for cnt in range(max_num_dis):
            if abs(ave_SC_dis[cnt])>1e-12:
                print(" %f %f %f %d " % (math.sqrt(cnt),abs(ave_SC_dis[cnt]),err_SC_dis[cnt],cnt) , file=f)

if __name__ == "__main__":
    main()
