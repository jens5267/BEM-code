#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:11:28 2019

@author: jere
"""
import numpy as np
import pandas as pd
import csv




def interpolate():
    Reynolds = []
    for i in range(2):
        Reynolds.append(50+i*5)
    for re in Reynolds:
        Re = 'Re'+'{:.3f}'.format(re/1000)
        print(Re)
    with open ('Airfoils/polars/NACA_2418/test/NACA_2418_T1_'+Re+'_M0.00_N9.0_360_M.dat','r') as fp:    
        ref = []
        position = 0
        for i in range(170):
            ref.append(-180.0+i)
        for j in range(300):
            ref.append(-10.0+(j/10))
        for k in range(161):
            ref.append(20.0+k)
        ref.append(180)

        header = ''
        for i in range (14):
            header += (next(fp)) 
        data = fp.read()
        alpha = []
        c_L = []
        c_D = []
        for line in data.splitlines():
            temp = line.split()
            if temp:
                alpha.append(eval(temp[0]))
                c_L.append(eval(temp[1]))
                c_D.append(eval(temp[2]))    
            else: 
                break
    for i in range(len(ref)):
        
        if alpha[i] == (ref[i]) or abs(alpha[i]-ref[i])<0.09:
            continue
        else:
            alpha.insert(i,float('{:.1f}'.format(ref[i])))
            c_L.insert(i, np.nan)
            c_D.insert(i,np.nan)

    c_L = pd.Series(c_L).interpolate().values.ravel().tolist()
    c_D = pd.Series(c_D).interpolate().values.ravel().tolist()
    
    
    with open ('Airfoils/polars/NACA_2418/test/NACA_2418_interpolated.dat','w') as f:
        print(header, file = f)
        for j in range(len(alpha)):
            print('{:10.2f}{:10.2f}{:10.2f}'.format(alpha[j], c_L[j], c_D[j]), file = f)
#
#, file = f
interpolate()
