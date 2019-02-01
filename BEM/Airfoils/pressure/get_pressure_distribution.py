# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 12:55:13 2018

@author: Jens Reichwein
"""

import numpy as np
import matplotlib.pyplot as plt

def pressure():
    NACA        = [2418]#,2412,2415,2418,4408,4412,4415,4418]
    S           = [833]#822,1210,2091,4061,4180,4320]
    
    reynolds    = []
    paths       = []
    
    for n in range(2):
        reynolds.append(50+n*50) #50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120
    
    for naca in NACA:
        paths.append(str('NACA_'+str(naca)))
        
    #for sg in SG:
    #    paths.append(str('SG'+str(sg)))
        
    #for s in S:
    #    paths.append(str('S'+str(s)))
        
    #for sd in SD:
    #    paths.append(('SD'+str(sd)))
    
    
    #x_list = []
    #reynolds_list = []
    AoA = [-10,0,10]
    linestyles  = [':', '-', '-.']# , ':'
    f = 0

    for r in reynolds:
        f += 1

        for p in paths:
            airfoil = p
            for aoa in AoA:
                header = ''
                if aoa == -10:
                    add = '_m10'
                    alpha = -10
                    l = 0
                elif aoa == 0:
                    add = '_0'
                    alpha = 0
                    l = 1
                elif aoa == 10:
                    add = '_p10'
                    alpha = 10
                    l = 2
                      
                link = str(p) + '/'+ str(p) + '_Re' + str('{:.3f}'.format((int(r)/1000))) + add+'.txt'
                with open(link, 'r') as fp:
                    for j in range(6):
                        header += next(fp)
                    data = fp.read()
                    X           = []
                    c_pi        = []
                    c_pv        = []
    
                    for line in data.splitlines():
                        temp    = line.split()
                        if temp:                
                            x       = float(eval(temp[0]))
                            cpi     = float(eval(temp[1]))
                            cpv     = float(eval(temp[2]))
                            
                            c_pi.append(cpi)
                            c_pv.append(cpv)
                            X.append(x)
                
                
                plt.figure(f)
                plt.title(airfoil +' '+ 'Re' + ' = '+ str('{:,.0f}'.format(reynolds[f-1]*1000)))
                plt.plot(X, c_pv,color = 'k', linestyle = linestyles[l], label= r'$\alpha$'+ ' = '+ '{:>6s}'.format(str(alpha)+'°') )
                plt.gca().invert_yaxis()
                plt.legend()
                plt.xlim(-0.05, 1.05)
                plt.ylim(1.2, -4.2)
                plt.minorticks_on()
                plt.grid(True,which='major')
                plt.grid(True,which='minor', linestyle='--')  
                plt.axhline(0, color='k', lw=1)
                plt.show()

pressure()