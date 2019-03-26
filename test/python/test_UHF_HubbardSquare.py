import unittest
import os
import filecmp
import numpy as np
import shutil
import glob
import math
import subprocess

dir="data/UHF_HubbardSquare"

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

class TestAtomic(unittest.TestCase):

    def test_HeisenbergChain(self):

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
        array_calc = read_opt_def("./zqp_APOrbital_opt.dat", 32)
        ref_ave = read_opt_def("%s/ref/zqp_APOrbital_opt.dat" %dir, 32)
        testArray = (abs(array_calc-ref_ave) < 1e-8)
        for _test in testArray:
            self.assertTrue(_test.all())

         #clean directory
        subprocess.call("rm ./*.dat", shell = True)
        subprocess.call("rm ./*.def", shell = True)
        
if __name__ == '__main__':
    unittest.main()
