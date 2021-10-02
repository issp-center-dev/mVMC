import numpy as np
import os
import copy
import math
import cmath
import random

with open('zqp_opt.dat') as f:
    data = f.read()
    data = data.split()
    #print(data)
    #print(len(data))
cnt_max = len(data)
var_org = np.zeros([cnt_max], dtype=np.float)
var_new = np.zeros([cnt_max], dtype=np.float)

random.seed(0)
for i in range(0,cnt_max):
    var_org[i] = float(data[i])
    if(abs(float(data[i]))>1e-8) :
         var_new[i] = var_org[i]+random.uniform(-0.02,0.02)
    #print(var_org[i],var_new[i])
#[s] count not empty elements

with open('initial.def','w') as f:
    for i in range(0,cnt_max):
        print("%.12lf " % (var_new[i]),end='',file=f)
    print("  ",file=f)
