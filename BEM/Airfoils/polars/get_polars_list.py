# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 11:44:37 2018

@author: Jens Reichwein
"""

import matplotlib.pyplot as plt
import numpy as np

def airfoil_polars():
   
    Profile = 'NACA_2418'
    Re      = [55]#,55,60]#,65,70,75,80,85,90,95,100]

    max_cl      = [];           max_cd      = [];       max_L2D     = []
    AoA         = [];           Re_list     = []
    for r in Re:
        Alpha       = []
        c_L         = [];           c_D         = [];       L2D         = []         
        link = str(Profile) +'/Montgomerie/'+ '/'+ str(Profile) + '_T1_Re' + str('{:.3f}'.format((int(r)/1000))) + '_M0.00_N9.0_360_M.dat'
        with open(link, 'r') as fp:
            for i in range (14):
                next(fp)
            data = fp.read()

            for line in data.splitlines():
                temp    = line.split()
                if temp:                
                    alpha   = eval(temp[0]);    cl      = eval(temp[1]);    cd      = eval(temp[2]);    l2d = cl/cd
                    Alpha.append(alpha);        c_L.append(cl);             c_D.append(cd);             L2D.append(l2d)                       
                else:
                    break
        position    = L2D.index(max(L2D[50:]))
        aoa = Alpha[position];  cl = c_L[position];   cd = c_D[position];   max_l2d     = L2D[position]
        AoA.append(aoa);        max_cl.append(cl);    max_cd.append(cd);    max_L2D.append(max_l2d) 
        
        start = 220
        end= 400
    
#        plt.figure(1)
#        plt.plot(Alpha[start:end], c_L[start:end])
#        plt.figure(2)
#        plt.plot(Alpha[start:end], c_D[start:end])        
        plt.figure(3)
        plt.ylim(8, 15)
        plt.plot(Alpha[start:end], L2D[start:end])
        plt.minorticks_on()
        plt.grid(True,which='major')
        plt.grid(True,which='minor', linestyle='--')            
airfoil_polars()