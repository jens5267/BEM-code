# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 12:23:04 2018

@author: Jens Reichwein
"""

#def machines(screws, a_ax=1):
#    '''
#
#    ''' 
#    machine = 0
#    count = 0
#    
#    while screws > 54:
#        machine += 1
#        screws -= 54
#        count += 1
#        a_ax -= 0.1
#    
#    return machine, screws, count, a_ax
import numpy as np
import pandas as pd

airfoil = 'SG6043'
Reynolds = 100000
tests = []
nacas   = [2408,2412,2415,2418,4408,4412,4415,4418]
sg      = [6040,6041,6042,6043]
s       = [822,1210,2091,4061,4180,4320,]
sd      = [7034] 

 
reynolds = [50,75,100]

for naca in nacas:
    for r in reynolds:
        tests.append(str('NACA_'+str(naca)+'_'+str(r))+'k.txt')
        
for i in sg:
    for r in reynolds:
        tests.append(str('SG'+str(i)+'_'+str(r))+'k.txt')

for i in s:
    for r in reynolds:
        tests.append(str('S'+str(i)+'_'+str(r))+'k.txt')
        
for i in sd:
    for r in reynolds:
        tests.append(str('SD'+str(i)+'_'+str(r))+'k.txt')
#print(test)

#with open(airfoil+'_'+str(int(Reynolds/1000))+'k.txt', 'r') as fp:
for test in tests:
    with open(test, 'r') as fp:    
        for i in range (12):
            next(fp)
        data = fp.read()
        Alpha   = []
        c_L     = []
        c_D     = []
        L2D     = []
        ref     = []  
        values  = []
        
        for line in data.splitlines():
            temp    = line.split()
            if temp:                
                alpha   = eval(temp[0])
                cl      = eval(temp[1])
                cd      = eval(temp[2])
                
                c_L.append(float(temp[1]))
                c_D.append(float(temp[2]))
                Alpha.append(float(alpha))
    
        for angle in range(301):
            ref.append(-100+angle)
            ref[angle]=ref[angle]/10
        for element in ref:
            if element not in Alpha:
                values.append(element)
        for i in range(len(ref)):
            if Alpha[i]!=ref[i]:
                Alpha.insert(i,ref[i])
                c_L.insert(i,np.NaN)
                c_D.insert(i,np.NaN)
        for i in range(len(c_L)):
            if Alpha[i] == min(Alpha) and np.isnan(c_L[i]):
                c_L[i] = c_L[i+1]
                c_D[i] = c_D[i+1]
            
            elif np.isnan(c_L[i]):
                j = 0
                while np.isnan(c_L[i+j]):
                    CL = c_L[i-1]+(c_L[i+j+1]-c_L[i-1])/(Alpha[i+j+1]-Alpha[i-1])*(Alpha[i]-Alpha[i-1])
                    CD = c_D[i-1]+(c_D[i+j+1]-c_D[i-1])/(Alpha[i+j+1]-Alpha[i-1])*(Alpha[i]-Alpha[i-1])
                    j += 1
                c_L[i] = CL
                c_D[i] = CD
#
#for i in range(len(ref)):
#    print('{:6.2f},{:8.3f},{:8.3f}'.format(Alpha[i], c_L[i], c_D[i]))      
with open(airfoil+'_'+str(int(Reynolds/1000))+'k_interpolated.txt','w') as f:
    for i in range(len(ref)):
        f.write('{:6.2f}\t{:7.3f}\t{:7.3f} \n'.format(Alpha[i], c_L[i], c_D[i]))