# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 16:15:48 2018

@author: Jens Reichwein
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 13:23:59 2018

@author: Jens Reichwein
"""
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

def interpolation():
    
    path       = 'AWT/cp_lambda'
    file = 'cp_lambda_original'
    linestyles = ['-','-.']
    
    with open(path +'/'+ file +'.txt', 'r') as fp:
        header = ''
        for j in range(4):
            header += next(fp)

        data = fp.read()
        TSR         = []
        c_P         = []
        c_P2        = []
                    
        for line in data.splitlines():
            
            temp    = line.split()
            if temp:      
                cp = (eval(temp[1]))
                if cp >=0:
                    TSR.append(float(eval(temp[0])))
                    c_P.append(float(cp))
#        for j in c_P:
#            c_p = savgol_filter(j, 51, 3)
#        c_P2.append(c_p)
        plt.title('Power Coefficient vs Tip Speed Ratio')
        plt.plot(TSR, c_P, color = 'k')
#    plt.legend()
    plt.minorticks_on()
    plt.grid(True,which='major')
    plt.grid(True,which='minor', linestyle='--')
    plt.xlabel('TSR')
    plt.ylabel(r'$c_P$')
    plt.xlim(0)
    plt.ylim(0)
    plt.show()
                           
interpolation()
