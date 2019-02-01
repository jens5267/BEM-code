# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 09:28:04 2019

@author: Jens Reichwein
"""
import matplotlib.pyplot as plt
import numpy as np

def compare_airfoils():
    airfoils = ['NACA_2415','NACA_2418','S833']#,NACA_2415 
    Linestyle = ['-','--', ':']
    l = -1
    Re = []
    for r in range(31):
        Re.append(50+r*5)
    for r in range(30):
        Re.append(210+r*10)
    
#    for airfoil in airfoils:
#        max_L2D     = []
#        max_Alpha   = []
#        l += 1
#        for re in Re:
#            data = airfoil+'_T1_Re'+str('{:.3f}'.format(re/1000))+'_M0.00_N9.0_360_M.dat'
#            link = 'Airfoils/polars/'+ airfoil + '/Montgomerie/'+ data #str(p) + '/'+ s+'/'+ str(p) + '_T1_Re' + str('{:.3f}'.format((int(r)/1000))) + '_M0.00_N9.0_360'+abb+'.dat'
#            header = ''
#            with open(link, 'r') as fp:
#                for j in range(14):
#                    header += next(fp)
#                data = fp.read()
#                Alpha       = []
#                c_L         = []
#                c_D         = []
#                L2D         = []
#
#                
#                for line in data.splitlines():
#                    temp    = line.split()
#                    if temp:                
#                        Alpha.append(eval(temp[0]))
#                        c_L.append(eval(temp[1]))
#                        c_D.append(eval(temp[2]))
#                
#                for i in range(len(Alpha)):
#                    L2D.append(c_L[i]/c_D[i])
#                max_L2D.append(max(L2D[10:]))
#                max_Alpha.append(Alpha[L2D.index(max(L2D[50:]))])
#                
#        mean = sum(max_Alpha)/len(max_Alpha)
#        s = []
#        for xi in max_Alpha:
#            s.append(((xi-mean)**2))
#        std = np.sqrt((1/len(max_Alpha))*sum(s))
#    l = -1   
    for airfoil in airfoils:
        l += 1
        x       = []
        y       = []            
        link = 'Airfoils/geometry/'+ airfoil + '/'+  airfoil+'.dat'           
        header = ''
        with open(link, 'r') as fp:
            data = fp.read()
            for line in data.splitlines():
                temp    = line.split()
                if temp: 
                    x.append(eval(temp[0]))
                    y.append(eval(temp[1]))
                
#        plt.minorticks_on()
#        plt.grid(True,which='major')
#        plt.grid(True,which='minor', linestyle='--')
#        plt.xlim(50, 500)
                     
#        plt.figure(1)
#        plt.plot(Re, max_Alpha, label = airfoil, color = 'k', linestyle = Linestyle[l])
#        plt.plot([50,500], [sum(max_Alpha)/len(max_Alpha),sum(max_Alpha)/len(max_Alpha)], color = 'grey', linestyle = Linestyle[l])
#        plt.ylim(0,12.5)
#        plt.xlabel('Reynolds number')
#        plt.ylabel('AoA at maximum L2D')
#        plt.legend()
#        
#        plt.figure(2)
#        plt.plot(Re, max_L2D, label = airfoil, color = 'k', linestyle = Linestyle[l])
#        plt.plot([50,500], [sum(max_L2D)/len(max_L2D),sum(max_L2D)/len(max_L2D)], color = 'grey', linestyle = Linestyle[l])
#        plt.ylim(0,100)
#        plt.xlabel('Reynolds number')
#        plt.ylabel('Maximum Lift to Drag ratio L2D')
#        plt.legend()
#        if airfoil == 'NACA_2418':
        plt.plot(x,y, 'k',linestyle = Linestyle[l], label = airfoil)
        plt.ylim(-0.5,0.5)
        plt.xlim(0,1)
        plt.axis('off') 
        plt.legend()
        plt.show()
#        print('standard deviation: ', std)
#        print('average L2D: ', sum(max_L2D)/len(max_L2D), airfoil)
#        print('average AoA: ', sum(max_Alpha)/len(max_Alpha), airfoil)
            
compare_airfoils()