# python3

import scipy.odr as odr
import numpy as np

import sys
sys.path.append("D:\BIO\PhD\ditalia\zebrafish\github\ditalia-zebrafish\code")
print(sys.path)

import mymod
f = mymod.fit_func_py

xdata = [1,2,3,4,5]
ydata = [1,4,9,16,22]
power_py = odr.Model(mymod.fit_func_py)
mydata = odr.Data(xdata, ydata)
myodr = odr.ODR(mydata, power_py, beta0=np.array([1, 1]))
myoutput = myodr.run();
myoutput.pprint()