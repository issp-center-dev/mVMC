import subprocess
import filecmp
import numpy as np
import os

def read_out(filename):
    # drop the first two columns
    array = np.loadtxt(filename, dtype='float').astype('float')
    return array
        
if __name__ == '__main__':
    cdir = os.getcwd()
#    subprocess.call("sh "+cdir+"/make_ref.sh", shell = True)
    array = read_out("output1/zqp_opt.dat")
    for i in range(2, 11):
        filename = "output{0}/zqp_opt.dat".format(i) 
        array = np.append(array, read_out(filename))
    array = array.reshape(10, array.shape[0]//10)
    std = np.std( array, axis = 0 )
    ave = np.mean( array, axis = 0 )
    np.savetxt("ref_std.dat", X = std)
    np.savetxt("ref_mean.dat", X = ave)
    #subprocess.call("rm -rf "+cdir+"/output*", shell = True)
