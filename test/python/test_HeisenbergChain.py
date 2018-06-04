import unittest
import os
import filecmp
import numpy as np

dir="data/HeisenbergChain"

def read_out(filename):
    # drop the first two columns
    array = np.loadtxt(filename, dtype='float').astype('float')
    return array

class TestAtomic(unittest.TestCase):
    def test_HeisenbergChain(self):
        # run
        #self.assertIs(0, os.system("../../src/mVMC/vmc.out -s %s/StdFace.def" %dir))
        
        # get results
        array_calc = read_out("./output/zqp_opt.dat")[0:2]
        ref_ave = read_out("%s/ref/ref_mean.dat" %dir)[0:2]
        ref_std = read_out("%s/ref/ref_std.dat" %dir)[0:2]
        testFalse = (ref_std < abs(array_calc - ref_ave))
        print(abs(array_calc[testFalse]) < 1e-8)
        print(abs(array_calc[testFalse]), abs(ref_ave[testFalse]), abs(ref_std[testFalse]))

        
if __name__ == '__main__':
    unittest.main()
