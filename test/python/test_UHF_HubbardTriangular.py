import unittest
import os
import filecmp
import numpy as np
import shutil
import glob
import math
import subprocess

dir="data/UHF_HubbardTriangular"

def read_opt_def(filename, total_site):
    rf = open(filename, "r")
    _arr = np.zeros(total_site, dtype=complex)    
    for line in rf.readlines()[5:]:
        line1 = line.split()
        _arr[int(line1[0])]=float(line1[1])+1j*float(line1[2])
    rf.close()
    return _arr

def read_out(filename):
    array = np.loadtxt(filename, dtype='float').astype('float')
    return array


def read_trans_def(filename, total_site):
    rf = open(filename, "r")
    _arr = np.zeros(int(math.pow(total_site*2, 2))*2,dtype=complex).\
                 reshape((total_site, 2, total_site, 2, 2))    
    for line in rf.readlines()[5:]:
        line1 = line.split()
        _arr[int(line1[0])][int(line1[1])][int(line1[2])][int(line1[3])]=float(line1[4])+1j*float(line1[5])
    rf.close()
    return _arr

class TestAtomic(unittest.TestCase):

    def test_HubbardTri(self):

        # run
        for _file in glob.glob("%s/*.def" %dir):
            file_name = os.path.split(_file)
            shutil.copyfile(_file, "./%s" %file_name[1])
        
        # run
        self.assertIs(0, os.system("../../src/ComplexUHF/UHF %s/namelist.def" %dir))        
        # get results
        array_calc = read_out("./zvo_eigen.dat")
        ref_ave = read_out("%s/ref/zvo_eigen.dat" %dir)
        testArray = (abs(array_calc-ref_ave) < 1e-8)
        for _test in testArray:
            self.assertTrue(_test.all())

        # get results
        #array_calc = read_opt_def("./zqp_APOrbital_opt.dat", 64)
        #ref_ave = read_opt_def("%s/ref/zqp_APOrbital_opt.dat" %dir, 64)
        #testArray = (abs(array_calc-ref_ave) < 1e-8)
        #for _test in testArray:
        #    self.assertTrue(_test.all())

        # get results
        array_calc = read_trans_def("./zvo_UHF_cisajs.dat", 8)
        ref_ave = read_trans_def("%s/ref/zvo_UHF_cisajs.dat" %dir, 8)
        testArray = (abs(array_calc-ref_ave) < 1e-8)
        for _test in testArray:
            self.assertTrue(_test.all())
            
         #clean directory
        subprocess.call("rm ./*.dat", shell = True)
        subprocess.call("rm ./*.def", shell = True)
        
if __name__ == '__main__':
    unittest.main()
