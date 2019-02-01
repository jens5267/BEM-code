# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 09:47:49 2019

@author: Jens Reichwein
"""
import numpy as np
import matplotlib.pyplot as plt
def airfoil_selection():
    Profiles = ['AS_5048','Dornier_A5','Eppler_560','Eppler_1233','FX_83_W_160','GOE_508','GOE_513','GOE_535','GOE_626','LS_417','MS_317','NACA_2418','S_825','S833']#,'S833',
    Reynolds = [50,75,100,125]
    Colors = ['r', 'b', 'g', 'purple', 'k', 'y']
    Markers = ['.','+','x','^','o', '1','2','3']
    av_alpha    = []
    av_l2d      = []
    std         = []
    var         = []
    L2D_max     = []
    Alpha_max   = []
    L_pm_5      = []; A_pm_5 = []
    for Profile in Profiles:                                                    # reads in data from XFoil for different profiles
        alpha_max   = []
        l2d_max     = []
        av_a = []
        av_l = []
        
        for re in Reynolds: 
            l_pm_5  = []
            a_pm_5  = []
            Re = 'Re'+'{:.3f}'.format((re)/1000) 
#            inp = 'Airfoils/polars/'+Profile+'/Montgomerie/'+Profile+'_T1_'+Re + '_M0.00_N9.0_360_M.dat'
            inp = 'Airfoils/polars/airfoil_data/'+Profile+'_T1_'+Re + '_M0.00_N9.0_360_M.dat'
            with open(inp, 'r') as fp:
                header = ''
                for i in range (14):
                    header += next(fp)            
                data = fp.read()
                
                alpha   = []
                c_L     = []
                c_D     = []
                L2D     = []
                for line in data.splitlines():
                    temp = line.split()
                    if temp:
                        alpha.append(eval(temp[0]))
                        c_L.append(eval(temp[1]))
                        c_D.append(eval(temp[2]))  
                    else: 
                        break
                for i in range(len(c_L)):
                    L2D.append(c_L[i]/c_D[i])
            position = L2D.index(max(L2D[120:-120]))
            Alpha = alpha[position]
            l2d = L2D[position]
            alpha_max.append(Alpha)
            l2d_max.append(l2d)
            for pm in range(11):
                l_pm_5.append(L2D[position-5+pm])
                a_pm_5.append(alpha[position-5+pm])
        L_pm_5.append(l_pm_5)
        A_pm_5.append(a_pm_5)
        
        Alpha_max.append(alpha_max)
        L2D_max.append(l2d_max)
        av_a = (sum(alpha_max)/len(alpha_max))
        av_l = sum(l2d_max)/len(l2d_max)
        av_alpha.append(av_a)
        av_l2d.append(av_l)
        s2 = []
        for m in range(len(alpha_max)):
            s2.append((alpha_max[m]-av_a)**2)
        std.append(sum(s2)/len(s2))
    for s in range(len(std)):
        var.append(np.sqrt(std[s]))
    print('{:12s} {:>8s} {:>8s} {:>8s} {:>8s} {:s}'.format('Profile', 'av AoA', 'av L2D', 'std dev', 'variance', 'L2D_max'))        
    c = -1
    for l in range(len(Profiles)):
        stage1 = av_l2d[l]>sum(av_l2d)/len(av_l2d)
        stage2 = var[l] <=1
        stage3 = True#(abs(max(L_pm_5[l])-min(L_pm_5[l]))) < 3
        if stage1 and stage2 and stage3:                    # check which airfoil is worse than average and if variance is <=1 
            c += 1
            print('{:12s} {:8.2f} {:8.2f} {:8.2f} {:8.2f}'.format(Profiles[l],av_alpha[l],av_l2d[l], std[l], var[l]), L2D_max[l])
            plt.plot(A_pm_5[l],L_pm_5[l], label= Profiles[l], color = 'k', marker = Markers[c])
            plt.legend()
            plt.xlabel('-0.5° < max AoA < +0.5°')
            plt.ylabel('max glide ratio')
            plt.xlim(0,10)
            plt.ylim(30,55)
    print('average AoA:', '{:8.2f}'.format(sum(av_alpha)/len(av_alpha)),'average L2D', '{:8.2f}'.format(sum(av_l2d)/len(av_l2d)))     
airfoil_selection()