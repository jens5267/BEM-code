# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 09:47:40 2018

@author: Jens Reichwein
"""
import matplotlib.pyplot as plt
import numpy as np

def airfoil_geometry():
    airfoil     = 'S822'
    Reynolds = [50000,75000,100000]
    
    Reynolds1    = 50000
    Reynolds2    = 75000
    Reynolds3    = 100000
    
#    Reynolds_new = 65000
        
    path        = 'polars/' + airfoil 
    polar1       = airfoil + '_' + str(int(Reynolds1/1000)) + 'k' 
    polar2       = airfoil + '_' + str(int(Reynolds2/1000)) + 'k' 
    polar3       = airfoil + '_' + str(int(Reynolds3/1000)) + 'k'     

    with open(path + '/'+ polar1 + '.txt', 'r') as fp:
        for i in range (12):
            next(fp)
        data = fp.read()
        Alpha1   = []
        c_L1     = []
        c_D1     = []
        L2D1     = []
        rev      = []
        count   = 0
        
        for line in data.splitlines():
            temp    = line.split()
            #for i in range (5):
            if temp:
                alpha1   = eval(temp[0])
                cl1      = eval(temp[1])
                cd1      = eval(temp[2])
                
                c_L1.append(temp[1])
                c_D1.append(temp[2])
                Alpha1.append(alpha1)
                
                l2d1 = cl1/cd1
                L2D1.append(l2d1)
            else:
                break
        position1 = L2D1.index(max(L2D1))
        angle1 = Alpha1[position1]



###############################################################################
    with open(path + '/'+ polar2 + '.txt', 'r') as fp:
        for i in range (12):
            next(fp)
        data = fp.read()
        Alpha2   = []
        c_L2     = []
        c_D2     = []
        L2D2     = []
        for line in data.splitlines():
            temp    = line.split()
            alpha2   = eval(temp[0])
            cl2      = eval(temp[1])
            cd2      = eval(temp[2])
            
            c_L2.append(temp[1])
            c_D2.append(temp[2])
            Alpha2.append(alpha2)
            
            l2d2 = cl2/cd2
            L2D2.append(l2d2)
        
        position2 = L2D2.index(max(L2D2))
        angle2 = Alpha2[position2]        
        
###############################################################################
    with open(path + '/'+ polar3 + '.txt', 'r') as fp:
        for i in range (12):
            next(fp)
        data = fp.read()
        Alpha3   = []
        c_L3     = []
        c_D3     = []
        L2D3     = []
        for line in data.splitlines():
            temp    = line.split()
            alpha3   = eval(temp[0])
            cl3      = eval(temp[1])
            cd3      = eval(temp[2])
            
            c_L3.append(temp[1])
            c_D3.append(temp[2])
            Alpha3.append(alpha3)
            
            l2d3 = cl3/cd3
            L2D3.append(l2d3)
        
        position3 = L2D3.index(max(L2D3))
        angle3 = Alpha3[position3]        

 
        
        
    
    
    
###############################################################################
    print('Airfoil: {:s}'.format(airfoil))
    print('{:>8s} {:>8s} {:>8s}'.format('max L/D','alpha', 'Re'))
    print('{:8.3f} {:8.3f} {:8d}'.format(max(L2D1),angle1, Reynolds1))
    print('{:8.3f} {:8.3f} {:8d}'.format(max(L2D2),angle2, Reynolds2))
    print('{:8.3f} {:8.3f} {:8d}'.format(max(L2D3),angle3, Reynolds3))        

    
airfoil_geometry()