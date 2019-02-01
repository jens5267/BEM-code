# -*- coding: utf-8 -*-
"""
Created on Fri Nov 23 10:56:57 2018

@author: Jens Reichwein
"""

import numpy as np
import matplotlib.pyplot as plt
def interpolation():
    plt.close('all')
    airfoil     = []
    tests       = []
    paths       = []
    NACA    = [2418]#2408,2412,2415,4408,4412,4415,4418
    SG      = [6042]#6040,6041,,6043
    S       = [822,1210,2091,4061,4180,4320]
    SD      = [7034] 
    
    reynolds    = [50,75,100]
    a1 = []
    a2 = []
    cl1 = []
    cl2 = []
    cd1 = []
    cd2 = []
    
    for naca in NACA:
        for r in reynolds:
            tests.append(str('NACA_'+str(naca)+'_'+str(r))+'k.txt')
        paths.append(str('NACA_'+str(naca)))
        
    #for sg in SG:
    #    for r in reynolds:
    #        tests.append(str('SG'+str(sg)+'_'+str(r))+'k.txt')
    #    paths.append(str('SG'+str(sg)))
    
    #for s in S:
    #    for r in reynolds:
    #        tests.append(str('S'+str(s)+'_'+str(r))+'k.txt')
    #    paths.append(str('S'+str(s)))
        
#    for sd in SD:
#        for r in reynolds:
#            tests.append(str('SD'+str(sd)+'_'+str(r))+'k.txt')
#        paths.append(('SD'+str(sd)))
    
    states      = ['lift', 'drag', 'L2D']    
    linestyles  = [':', '--', '-.']# , ':'
    markers     = ['.','+','x']
    colors      = ['k', 'dimgrey','lightgrey'] 
    interpolation = True
    if interpolation == True:
        text = 'k_interpolated.txt'
        skip_data = 0
        fig = 1
    elif interpolation == False:
        text = 'k.txt'
        skip_data = 12
        fig = 4
    l = -1
    c = -1
    s = -1

    for r in reynolds:
        Alpha   = []
        c_L     = []
        c_D     = []
        L2D     = []
        s += 1
        Re = r
        l += 1
        m = -1
        c += 1
        for p in paths:
            m += 1

            with open(str(p) + '/'+ str(p) + '_' + str(r) + text , 'r') as fp:
                for i in range(skip_data):
                    next(fp)
                data = fp.read()

    
                for line in data.splitlines():
                    temp    = line.split()
                    if temp:                
                        alpha   = eval(temp[0])
                        cl      = eval(temp[1])
                        cd      = eval(temp[2])
                        if cd == 0:
                            cd = 0.00001
                        c_L.append(float(temp[1]))
                        c_D.append(float(temp[2]))
                        
                        Alpha.append(float(alpha))
                        L2D.append(cl/cd)


            Label = ('{:5s} {:>10s}  {:<10s}'.format('Re =',(str(Re)+',000'), p))                
            for state in states:
                if state == 'drag':
                    plt.figure(fig)
                    plt.plot(Alpha, c_D, label = Label, linestyle = linestyles[l], marker = markers[m], color=colors[c])            
                    plt.ylabel('Drag coefficient c$_D$ [-]')
    
                elif state == 'lift':
                    plt.figure(fig+1)
                    plt.plot(Alpha, c_L, label = Label, linestyle = linestyles[l], marker = markers[m], color=colors[c])   
#                    plt.plot(alpha_new, cL_new)
                    plt.ylabel('Lift coefficient c$_L$ [-]')
    
                elif state == 'L2D':
                    #ax = plt.figure(3)
                    #ax.set_title('Test')
                    plt.figure(fig+2)
                    plt.plot(Alpha, L2D, label = Label, linestyle = linestyles[l], marker = markers[m], color=colors[c])
                    plt.ylim(-30,70)
                    plt.ylabel('Lift to drag ratio L2D [-]')
                
                plt.grid(True,which='major')
                plt.minorticks_on()
                plt.grid(True,which='minor', linestyle='--')
                plt.xlim(-10,20)
                plt.xlabel('Angle of Attack ' + r'$\alpha$' +'[°]')
                plt.legend()

    #            mng = plt.get_current_fig_manager()
                
    #            plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
                #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,ncol=4, mode="expand", borderaxespad=0.)
            #print(Label)
    
    
    #plt.show()
    #mng.window.setGeometry(0,30, 420,350)
   
interpolation()