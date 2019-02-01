# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 14:44:49 2018

@author: Jens Reichwein

Shows a GUI with input parameters. The parameters can be plotted, saved and the blade can be seen in 3D
The aerodynamic values of the blades are stored in .txt files for Re- numbers of 50.000. They can be generated using XFoil.
So far, only 4-digit NACA airfoils can be used.
A correction for number of blades could not be succesfully implemented, as the values for the tip showed impossible values.

Input: design wind speed, design rotations, design rotor radius, .txt files of aerodynamic values
Output: bladed design with chord and twist distribution. 
"""
import numpy as np
import matplotlib.pyplot as plt
import os

def aerodynamics():        
    Profile             = 'NACA_2418'
    Re                  = 100000
    rotor_diameter      = 0.4;     hub_diameter        = 0.056    
    rotations           = int(input('rotations: '))
    v_des               = 15
    profile             = 'NACA_2418';    method              = 'Glauert'
        
    sections            = 151;     B                   = 3
    efficiency          = 0.8
    
    frequency           = rotations /60;     Omega               = 2*np.pi*frequency
    
    nu                  = 15.11*10**-6
    rho                 = 1.225

    Radius              = rotor_diameter/2;     hub_radius          = hub_diameter/2
    dr                  = (Radius- hub_radius)/(sections-1)
    
    TSR                 = Omega * Radius / v_des
    v_tip               = Omega*Radius 

    A                   = np.pi * (Radius**2 -  hub_radius**2)
    P_w                 = 0.5* rho * v_des**3 * A;              P_t             = P_w * 16/27

    a_ax_list           = [];       a_rot_list          = [];       TSR_list        = []
    chord_list          = [];       Re_list             = [];       radius_list     = []
    c_P_list            = [];       c_T_list            = []     
    cl_list             = [];       cd_list             = [];       L2D_list        = []
    f_tip_list          = [];       f_hub_list          = [];       F_list          = []
    
    dT_list             = [];       dQ_list             = [];       dP_list         = []
    
    alpha_list          = [];       phi_deg_list        = [];       twist_list      = []; phi_rad_list = []
    v_app_list          = [];       iteration_list      = []

    for i in range(sections-1):  
        iteration = 0                                                                     
        radius          = hub_radius + (i*dr+(i+1)*dr)/2
        Re = 100000
        cl = 1; cd = 0.1; alpha = 0
        a_ax_init       = 0.0 # can be any value below 0.3, must be smaller than a_ax
        delta_a         = 1.0
        a_ax            = 0.0
        a_rot           = 0.001
        solidity        = 0.01
        
        chord           = 0
        count           = 0
        convergence_error = 10e-8
    
        while abs(delta_a)> convergence_error:
            iteration += 1
            v_rot          = radius*2*np.pi*frequency
            v_app          = np.sqrt((v_des*(1-a_ax))**2+(v_rot*(1+a_rot))**2) 
            
            TSR_local      = TSR * radius / Radius
    
            #phi_rad        = np.arctan(((1-a_ax)*v_des)/((1+a_rot)*(Omega*radius)))
            phi_rad = np.arctan(2/3/(TSR*radius/Radius*(1+2/(3*TSR**2*(radius/Radius)**2))))
            phi_deg        = np.rad2deg(phi_rad)
            
###############        Prandtl Correction             ##########################                 
            f_tip           = 2/np.pi*(np.arccos(np.exp(-(B/2 *((Radius-radius)/(radius * np.sin(phi_rad)))))))
            f_hub           = 2/np.pi*(np.arccos(np.exp(-(B/2 *((radius-hub_radius)/(radius * np.sin(phi_rad)))))))
            F               = f_tip*f_hub
            #F = 1
            chord           = (8*np.pi*radius*np.sin(phi_rad)/(3*B*cl*TSR_local))
            #chord          =  8*np.pi*radius/(B*cl)*(1-np.cos(phi_rad))
            #chord = ((2*np.pi)/(9*B*TSR*np.sqrt((4/9)+TSR**2*(radius/Radius)**2*((1+(2/(9*TSR**2*(radius/Radius)**2)))**2))))/cl
            cn             = cl*np.cos(phi_rad) + cd*np.sin(phi_rad)
            cr             = cl*np.sin(phi_rad) - cd*np.cos(phi_rad)        
            
            a_rot          = 1/((4*F*(np.sin(phi_rad)*np.cos(phi_rad)))/(cr*solidity)-1)   
            if a_rot <= 0:
                a_rot = 0
                       
            solidity       = B*chord/(2*np.pi*radius)  
            a_ax           = 1/(4*F*np.sin(phi_rad)*np.sin(phi_rad)/(solidity*cn)+1) 
###############        Correction Methods            ##########################            
            c_T1   = 1.816
            if a_ax >= 1/3:
                c_T     = 8/9+(4*F-20*c_T1/9)*a_ax +(25*c_T1/9-4*F)* a_ax**2
                a_ax    = (18*F-20-3*np.sqrt(c_T*(50-36*F)+12*F*(3*F-4)))/(36*F-50)
            else:
                c_T     = 4*F*a_ax*(1-a_ax)            
                        
            c_P = 4*a_ax*(1-a_ax)**2
            count += 1
            if abs(delta_a) >= 10**(-10):                              
                a_ax_init = (a_ax_init+a_ax)/2
            if count > 100:
                break    
            delta_a        = abs(a_ax_init - a_ax)
             
        conv = 0
        j = 0
        while conv <= 1:
            j += 1
            Alpha       = []
            c_L         = [];           c_D         = [];       L2D         = []         
            link = 'Airfoils/polars/'+str(Profile) +'/Montgomerie/'+ '/'+ str(Profile) + '_T1_Re' + str('{:.3f}'.format((int(Re)/10e5))) + '_M0.00_N9.0_360_M.dat'
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
                        
                if Re >200000:
                    base = 10000
                else:
                    base = 5000
                re = int(base * round(float((v_app*chord/nu))/(base)))                
                if re < 50000:
                    re = 50000 
                elif re >200000:
                    re = 200000
                Re = re
                
     
            if Re == re:
                conv += 1   
            
            dT = 4*F*(1-a_ax)*0.5*rho*v_des**2*2*np.pi*radius*dr
            dQ = abs(4*F*a_rot*(1-a_ax)*0.5*v_des*Omega*radius**2*2*np.pi*radius*dr)
            dP = Omega * dQ
            
            twist          = phi_deg-aoa
        cl_list.append(cl);    cd_list.append(cd);    L2D_list.append(max_l2d)
        Re_list.append(Re)            
        iteration_list.append(iteration)
        a_ax_list.append(a_ax),         a_rot_list.append(a_rot),       TSR_list.append(TSR_local)
        chord_list.append(chord),       radius_list.append(radius)
        c_P_list.append(c_P),           c_T_list.append(c_T)
        
        F_list.append(F),               f_tip_list.append(f_tip),       f_hub_list.append(f_hub)
        dT_list.append(dT),             dQ_list.append(dQ),             dP_list.append(dP)
        
        twist_list.append(twist),       phi_deg_list.append(phi_deg),   alpha_list.append(aoa),     phi_rad_list.append(phi_rad)
        v_app_list.append(v_app)
            
#    print('{:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}'.format('Section', "Radius", 'Phi deg','Phi rad', "a","a'","chord", 'v_app','Re', 'f_hub', 'f_tip', 'cl', 'TSR_l')) 
#    for l in range(len(a_ax_list)):
#        print('{:>8.0f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8,.0f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}'.format(l+1, radius_list[l], phi_deg_list[l], phi_rad_list[l], a_ax_list[l], a_rot_list[l], chord_list[l], v_app_list[l], Re_list[l], f_hub_list[l], f_tip_list[l], cl_list[l], TSR_list[l]))
  
    M_tot = sum(dQ_list)
    P_tot = sum(dP_list)*efficiency
    T_tot = sum(dT_list)

###############################################################################
###############          print           ######################################
#    print(50 * '--')
#    print('{:<15s}: {:>10.3f} m | {:<15s}: {:>10.0f}   | {:<15s}: {:>10.4f} cm'.format      ('Rot Diameter', rotor_diameter,'Sections', sections -1, 'delta r', dr*100))
#    print('{:<15s}: {:>10s}   | {:<15s}: {:>10.3f} m/s'.format          ('Profile', profile, 'des. wind speed', v_des))
#    print('{:<15s}: {:>10.0f}   | {:<15s}: {:>10.3f} rpm'.format          ('Blades', B, 'Rotations', rotations))
#    if v_tip > 80:
#        print('{:<15s}: {:>10.3f} m/s -> above 80 m/s'.format   ('Tip velocity', v_tip, ))  
#    else:
#        print('{:<15s}: {:>10.3f} m/s -> below 80 m/s'.format   ('Tip velocity', v_tip))
#    print('{:<15s}: {:>10.3f}   '.format       ('TSR', TSR, ))
#    print('{:<15s}: {:>10.3f} W | {:<15s}: {:>10.3f} W | {:<15s}: {:>10.3f} W | {:<15s}: {:>10.3f} %'.format('Power wind', P_w, 'Power Betz', P_t, 'Power WT', P_tot, 'Efficiency η', 100*P_tot/P_w))
#    print('{:<15s}: {:>10.3f} Nm| {:<15s}: {:>10.3f} N'.format('Torque WT', M_tot, 'Thrust WT', T_tot))  
#    print(50*'--')
#    plt.figure(1)
#    plt.plot(radius_list, phi_deg_list,   label = 'inflow')
#    plt.plot(radius_list, twist_list, label = 'Twist')
#    plt.plot(radius_list, alpha_list, label = 'Angle of attack')
#    plt.legend()
    plt.figure(2)
    plt.plot(radius_list, f_hub_list, 'k', linestyle = ':', label = 'f_tip', marker = '.')
    plt.plot(radius_list, f_tip_list, 'k', linestyle = '--',label = 'f_hub', marker = 'x')
    plt.plot(radius_list, F_list, 'k', linestyle = '-', label = 'F')
    plt.ylim(0,1.1)


#    plt.plot(iteration_list)
    #plt.xlim(0.025,0.2)
    plt.minorticks_on()
    plt.grid(True,which='major')
    plt.grid(True,which='minor', linestyle='--')    
    #plt.ylim(0,0.6)


aerodynamics()
