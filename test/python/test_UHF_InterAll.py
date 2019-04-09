import unittest
import os
import filecmp
import numpy as np
import shutil
import glob
import math
import subprocess
import itertools

dir_diag = "data/UHF_InterAll_Diagonal"    
dir_exchange =  "data/UHF_InterAll_Exchange"

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

def change_coulombIntra_to_Interall(filename, total_int):
    _arr_int_all = []
    with open(filename, "r") as rf:
        for line in rf.readlines()[5:]:
            line1 = line.split()
            _arr_int_all.append((int(line1[0]), 0, int(line1[0]), 0, int(line1[0]), 1, int(line1[0]), 1, float(line1[1]), 0.0))
    return _arr_int_all

def change_coulombInter_to_Interall(filename, total_int):
    _arr_int_all = []
    with open(filename, "r") as rf:
        for line in rf.readlines()[5:]:
            line1 = line.split()
            for i_spin, j_spin in itertools.product(range(2), range(2)):
                _arr_int_all.append((int(line1[0]), i_spin, int(line1[0]), i_spin,
                                              int(line1[1]), j_spin, int(line1[1]), j_spin, float(line1[2]), 0.0))
    return _arr_int_all

def change_exchange_to_Interall(filename, total_int):
    _arr_int_all = []
    with open(filename, "r") as rf:
        for line in rf.readlines()[5:]:
            line1 = line.split()
            _arr_int_all.append((int(line1[0]), 0, int(line1[1]), 0, int(line1[1]), 1, int(line1[0]), 1, float(line1[2]), 0.0))
            _arr_int_all.append((int(line1[0]), 1, int(line1[1]), 1, int(line1[1]), 0, int(line1[0]), 0, float(line1[2]), 0.0))
    return _arr_int_all


def MakeHeader(file, comment, number):
    file.write("========================\n")
    file.write("TotalNumber "+str(number)+"\n")
    file.write("Comment: "+str(comment)+"\n")
    file.write("========================\n")
    file.write("========================\n")
    return

def write_Interall(filename, arr_int):
    with open(filename, "w") as wf:
        MakeHeader(wf, "interall", len(arr_int))
        for list1 in arr_int:
            for str_list in list1:
                wf.write(str(str_list)+" ")
            wf.write("\n")
            
class TestAtomic(unittest.TestCase):

    def test_Diagonal(self):
        # run
        self.assertIs(0, os.system("../../src/mVMC/vmcdry.out %s/stan.in" %dir_diag))        
        for _file in glob.glob("%s/*.def" %dir_diag):
            file_name = os.path.split(_file)
            shutil.copyfile(_file, "./%s" %file_name[1])

        arr_int = change_coulombIntra_to_Interall("coulombintra.def", 4)
        arr_int += change_coulombInter_to_Interall("coulombinter.def", 4)
        write_Interall("interall.def", arr_int)
            
        # run
        self.assertIs(0, os.system("../../src/ComplexUHF/UHF %s/namelist_all.def" %dir_diag))        
        # get results
        array_calc = read_out("./zvo_eigen.dat")
        ref_ave = read_out("%s/ref/zvo_eigen.dat" %dir_diag)
        testArray = (abs(array_calc-ref_ave) < 1e-8)
        for _test in testArray:
            self.assertTrue(_test.all())

        # get results
        array_calc = read_trans_def("./zvo_UHF_cisajs.dat", 8)
        ref_ave = read_trans_def("%s/ref/zvo_UHF_cisajs.dat" %dir_diag, 8)
        testArray = (abs(array_calc-ref_ave) < 1e-8)
        for _test in testArray:
            self.assertTrue(_test.all())
        
         #clean directory
        subprocess.call("rm ./*.dat", shell = True)
        subprocess.call("rm ./*.def", shell = True)

    def test_Exchange(self):
        # copy
        for _file in glob.glob("%s/*.def" %dir_exchange):
            file_name = os.path.split(_file)
            shutil.copyfile(_file, "./%s" %file_name[1])

        arr_int = change_coulombIntra_to_Interall("coulombintra.def", 4)
        arr_int += change_exchange_to_Interall("exchange.def", 4)

        
        write_Interall("interall.def", arr_int)
            
        # run
        self.assertIs(0, os.system("../../src/ComplexUHF/UHF %s/namelist_all.def" %dir_exchange))        
        # get results
        array_calc = read_out("./zvo_eigen.dat")
        ref_ave = read_out("%s/ref/zvo_eigen.dat" %dir_exchange)
        testArray = (abs(array_calc-ref_ave) < 1e-8)
        for _test in testArray:
            self.assertTrue(_test.all())

        # get results
        array_calc = read_trans_def("./zvo_UHF_cisajs.dat", 8)
        ref_ave = read_trans_def("%s/ref/zvo_UHF_cisajs.dat" %dir_exchange, 8)
        testArray = (abs(array_calc-ref_ave) < 1e-8)
        for _test in testArray:
            self.assertTrue(_test.all())
        
         #clean directory
        subprocess.call("rm ./*.dat", shell = True)
        subprocess.call("rm ./*.def", shell = True)

        
if __name__ == '__main__':
    unittest.main()
