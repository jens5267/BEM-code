# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 13:23:59 2018

@author: Jens Reichwein
"""
import numpy as np
import matplotlib.pyplot as plt

def interpolation():
    NACA        = [2418]#,2412,2415,2418,4408,4412,4415,4418]
    S           = [833]#822,1210,2091,4061,4180,4320]
    
    reynolds    = []
    paths       = []
    L2D         = []
    for n in range(1):
        reynolds.append(100+n*50) #50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120
    
    for naca in NACA:
        paths.append(str('NACA_'+str(naca)))
        
    #for sg in SG:
    #    paths.append(str('SG'+str(sg)))
        
#    for s in S:
#        paths.append(str('S'+str(s)))
        
#    for sd in SD:
#        paths.append(('SD'+str(sd)))
    
    methods = ['Viterna', 'Montgomerie']
    linestyles = ['-','-.']
    diagrams = 'lift'#, 'drag', 'L2D']
    
    for r in reynolds:
        l = -1
        for s in methods:
            
            start = 170#0#
            end = 470#630#
            if s == 'Montgomerie':
                abb = '_M'
                
            elif s == 'Viterna':
                abb = '_V'
            for p in paths:
                l += 1
                header = ''
                link = str(p) + '/'+ s+'/'+ str(p) + '_T1_Re' + str('{:.3f}'.format((int(r)/1000))) + '_M0.00_N9.0_360'+abb+'.dat'
                with open(link, 'r') as fp:
                    for j in range(14):
                        header += next(fp)
                    data = fp.read()
                    Alpha       = []
                    c_L         = []
                    c_D         = []
                    ref         = []  
                    L2D         = []
                    max_L2D     = []
                    max_Alpha   = []
                    
                    for line in data.splitlines():
                        temp    = line.split()
                        if temp:                
                            alpha   = eval(temp[0])
                            cl      = eval(temp[1])
                            cd      = eval(temp[2])
                            
                            c_L.append(float(cl))
                            c_D.append(float(cd))
                            Alpha.append(float(alpha))
                            
                    for angle in range(170):
                        ref.append(int(-180.0+angle))
                    for angle2 in range(301):
                        ref.append(-10.0+(angle2/10))
                    for angle3 in range(160):
                        ref.append(21.0+angle3)
                    ref.append(180.0)
                            
    # =============================================================================
    #                 for i in range(len(ref)):
    #                     if abs(Alpha[i]-ref[i])>0.0001:
    #                         Alpha.insert(i,ref[i])
    #                         c_L.insert(i,np.NaN)
    #                         c_D.insert(i,np.NaN)
    #                 for i in range(len(c_L)):
    #                     if Alpha[i] == min(Alpha) and np.isnan(c_L[i]) and np.isnan(c_D[i]):
    #                         c_L[i] = c_L[i+1]
    #                         c_D[i] = c_D[i+1]
    # 
    #                     elif np.isnan(c_L[i]):
    #                         j = 0
    #                         while np.isnan(c_L[i+j]):
    #                             CL = c_L[i-1]+(c_L[i+j+1]-c_L[i-1])/(Alpha[i+j+1]-Alpha[i-1])*(Alpha[i]-Alpha[i-1])
    #                             CD = c_D[i-1]+(c_D[i+j+1]-c_D[i-1])/(Alpha[i+j+1]-Alpha[i-1])*(Alpha[i]-Alpha[i-1])
    #                             j += 1
    #                         c_L[i] = CL
    #                         c_D[i] = CD
    #                 
    # =============================================================================
#                    for aoa in range(len(Alpha)):
#                            l2d = c_L[aoa]/c_D[aoa]
#                            L2D.append(l2d)

#                    if c_D[i] == 0:
#                        pass
#                    else:    
#                        L2D = [ c_l/c_D[i] for i,c_l in enumerate(c_L) ]       
#                        alpha_pos = L2D.index(max(L2D[20:]))     
#                        max_alpha = Alpha[alpha_pos] 
#                        alpha_list.append(max_alpha)
#                        reynolds_list.append(r)
    #                max_alpha_pos = L2D.index(max(L2D))
    #                max_alpha = Alpha[int(max_alpha_pos)]
    #                max_Alpha.append(max_alpha)                
                if diagrams == 'lift':
                    ylabel = 'Lift Coefficient $c_L$'
                    coefficient = c_L
                    ylim = (-1,1.5)
                elif diagrams == 'drag':
                    ylabel = 'Drag Coefficient $c_D$'
                    coefficient = c_D
                    ylim = (0, 0.25)
                elif diagrams == 'L2D':
                    ylabel = 'Glide Ratio $L2D$' 
                    coefficient = np.array(np.array(c_L)/np.array(c_D))
                    ylim = (-30, 50)
            #plt.title(p)
            #plt.figure(methods.index(s))
            plt.plot(Alpha[start:end], coefficient[start:end], label = 'Re = '+str('{:,.0f}'.format(r*1000)) +' | '+s, color = 'k', linestyle = linestyles[l])
    plt.legend()
    plt.minorticks_on()
    plt.grid(True,which='major')
    plt.grid(True,which='minor', linestyle='--')
    plt.axhline(0, color = 'k', lw = 1)
    plt.axvline(0, color = 'k', lw = 1)
    plt.xlabel('Angle of Attack ' + r'$\alpha$')
    plt.ylabel(ylabel)
    plt.xlim(-10,20)
    plt.ylim(ylim)
#    plt.show()

    ALPHA = []
    lift_coeffcient = []
    drag_coefficient = []
    start = -10
    for i in range(300):
        ALPHA.append(start)
        lift_coeffcient.append(2*np.sin(start*np.pi/180)*np.cos(start*np.pi/180))
        drag_coefficient.append(2*np.sin(start*np.pi/180)*np.sin(start*np.pi/180))   
        start += 0.1
        if drag_coefficient[i] != 0 and lift_coeffcient[i] != 0:
            if (lift_coeffcient[i]/drag_coefficient[i]) <= 0:
                L2D.append(0)
            elif (lift_coeffcient[i]/drag_coefficient[i]) > 100:
                L2D.append(100)  
            else:                
                L2D.append(lift_coeffcient[i]/drag_coefficient[i])                  
    if diagrams == 'lift':
        plt.plot(ALPHA, lift_coeffcient, ':', color = 'k', label = '$c_L$ = 2$\sin$' r'($\alpha$) ' '$cos$' r'($\alpha$)')
    elif diagrams == 'drag':
        plt.plot(ALPHA, drag_coefficient, ':',color = 'k', label = '$c_D$ = 2$sin^2$' r'($\alpha$)')
    #elif diagrams 
    plt.legend()
#    plt.xlim(-10,20)
#    plt.minorticks_on()
#    plt.grid(True,which='major')
#    plt.grid(True,which='minor', linestyle='--')
#    plt.axvline(x=0, color='k', lw = 1)
#    plt.axhline(y=0, color='k', lw = 1)
#    plt.xlabel('Angle of Attack ' +r'$\alpha$')
#    plt.ylabel('Lift and Drag coefficients $c_L$ & $c_D$')

                           
#            link2 = str(p) + '/'+ str(p) + '_T1_Re' + str('{:.3f}'.format((int(r)/1000))) + '_M0.00_N9.0_360_M.dat'
#            with open(link2,'w') as f:
#                f.write(header)
#                for i in range(len(ref)):
#                    f.write('{:10.2f}\t{:10.4f}\t{:10.4f} \n'.format(Alpha[i], c_L[i], c_D[i]))      
#            for ltd in max_L2D:
#                print('{:6.2f} {:8,.0f} {:8.2f}'.format(ltd, r*1000, max_alpha))
#    print(max_alpha)
interpolation()

#for r in reynolds:
#    for p in paths:                    
#        link3 = str(p) + '/'+ str(p) + '_T1_Re' + str('{:.3f}'.format((int(r)/1000))) + '_M0.00_N9.0_360_M.dat'
#        with open(link3,'r') as fd:
#            for i in range (14):
#                next(fd)
#            data = fd.read()
#            Alpha   = []
#            c_L     = []
#            c_D     = []
#            L2D     = []
#            AoA     = []
#            max_L2D = []
#
#            for line in data.splitlines():
#                temp    = line.split()
#                if temp:                
#                    alpha   = eval(temp[0])
#                    cl      = eval(temp[1])
#                    cd      = eval(temp[2])
#                    
#                    c_L.append(cl)
#                    c_D.append(cd)
#                    Alpha.append(alpha)
#                    
#                    l2d = cl/cd
#                    L2D.append(l2d)
#                        
#                else:
#                    break
#                position = L2D.index(max(L2D))
#            AoA.append(Alpha[position])                                        
#            max_l2d = max(L2D)
#            max_L2D.append(max_l2d)
#        print((AoA), link3)