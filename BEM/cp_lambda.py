#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 13:17:33 2019

@author: jere
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 14:44:49 2018

@author: Jens Reichwein

Shows a GUI with input parameters. The parameters can be plotted, saved and the blade can be seen in 3D
The aerodynamic values of the blades are stored in .dat files for Re- numbers from 50,000 to 200,000 in steps of 5,000. They are generated using XFoil.

Input: design wind speed, design rotations, design rotor radius, .txt files of aerodynamic values
Output: bladed design with chord and twist distribution. 
"""
from mpl_toolkits.mplot3d import Axes3D
from tkinter import *
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import matplotlib.ticker as mtick

     
###############################################################################
def calculation(): 
    a_ax_list          = [];       a_rot_list         = [];       TSR_list           = [];       solidity_list      = []
    Re_list            = [];       re_list            = [] 
    
    max_alpha_list     = [];       max_cl_list        = [];       max_cd_list        = [];       max_l2d_list       = []
    c_T_list           = [];       c_Q_list           = [];       c_P_list           = []
    f_tip_list         = [];       f_hub_list         = [];       F_list             = []

    dT_list            = [];       dQ_list            = [];       dP_list            = []
    dT_list            = [];       dQ_list            = [];       P_list             = []

    alpha_list         = [];       phi_deg_list       = [];       phi_rad_list       = [];       twist_list         = []
    v_rot_list         = [];       v_app_list         = [];       iteration_list     = []  
###############################################################################
    
    chords              = [0.06893015616415156, 0.059761724798064124, 0.05080310987235474, 0.04349865668179045, 0.03774119286126487, 0.033191391487197054, 0.02954840695366159, 0.026585665528526488, 0.024139026188495713, 0.022089934897953255]
    twists              = [12.648078752835502, 10.728967262825156, 8.408273286836527, 6.260587042364667, 4.032691133394809, 1.8883273448830824, 0.6253317695802434, -0.3023779418871911, -0.5358723617040535, -0.47655194492204345]
    radii               = [0.0347, 0.05210000000000001, 0.0695, 0.0869, 0.1043, 0.1217, 0.13910000000000003, 0.1565, 0.17390000000000003, 0.1913]       
    sections            = [0,1,2,3,4,5,6,7,8,9]
    a_list = [0.49400907913701875, 0.38875099266262136, 0.35344941084770265, 0.33763069428058956, 0.33018121337645545, 0.3276320353529738, 0.3297332963308741, 0.3393118214206339, 0.36596182108452213, 0.45470111403226865]
    ap_list = [0.43316464840752167, 0.18816475514167863, 0.10460382006735092, 0.06643042271088116, 0.04594482488113761, 0.033759046095433114, 0.026035450872421706, 0.02103941135055198, 0.01806595482159363, 0.017692880251115557]
    
    rotor_diameter     = 0.4
    hub_diameter       = 0.056
    v_1                = 15
    
    Profile            = 'NACA_2418';                           
    method             = 'Wilson-Walker'

    efficiency         = 0.8
         
    nu                 = 15.11*10**-6;                            rho                = 1.225

    Radius             = rotor_diameter/2;                        hub_radius         = hub_diameter/2;         dr                 = (Radius- hub_radius)/(len(sections)-1)
    tsr = 0.1
    for i in range(10):
        frequency           = tsr * v_1/(2*np.pi*Radius)
        rotations           = frequency * 60                          
        Omega               = 2*np.pi*frequency  
        TSR                 = Omega * Radius / v_1
        v_tip               = Omega*Radius 

        A                  = np.pi * (Radius**2 - hub_radius**2)
        P_w                = 0.5* rho * v_1**3 * A
        P_t                = P_w * 16/27
        convergence_error  = 10e-6
        B                  = 3
        tsr += 0.5
    #################################################################################       
        for section in sections:
            radius             = radii[section]
            chord              = chords[section]
            twist              = twists[section] 

            dA                 = 2*np.pi*radius*dr
            v_rot              = radius*2*np.pi*frequency
            
            a_ax                = a_list[section]
            a_rot               = ap_list[section] 
            solidity           = B*chord/(2*np.pi*radius)
            
            Re_conv                = 0;        count              = 0
            Re                 = 100000;           re                     = 0
            
            cl                 = 1.0;              cd                     = 0.01
            while Re_conv != 1: 
                
                v_app          = np.sqrt((v_1*(1-a_ax))**2+(v_rot*(1+a_rot))**2) 
                #print(a_ax, a_rot,v_app)
                phi_rad        = abs(np.arctan(((1-a_ax)*v_1)/((1+a_rot)*(Omega*radius))))
                phi_deg        = np.rad2deg(phi_rad)
                
###############        Prandtl Correction             ##########################                 
                f_tip              = 2/np.pi*(np.arccos(np.exp(-(B/2 *((Radius-radius)/(radius * np.sin(phi_rad)))))))
                f_hub              = 2/np.pi*(np.arccos(np.exp(-(B/2 *((radius-hub_radius)/(radius * np.sin(phi_rad)))))))
                F                  = f_tip*f_hub
                
                cn             = cl*np.cos(phi_rad) + cd*np.sin(phi_rad)
                ct             = cl*np.sin(phi_rad) - cd*np.cos(phi_rad)  
                
                a_ax           = 1/((4*F*(np.sin(phi_rad)*np.sin(phi_rad)))/(cn*solidity)+1)   
                
                Re             = v_app*chord/nu  
                            
                if Re < 60000:
                    Re = 60000  
                elif Re > 500000:
                    Re = 500000
                if Re <= 200000:
                    base = 5000
                elif Re > 200000:
                    base = 10000
                     
                re             = int(base * round(float(Re)/base))  
#################################################################################                      
                c_T             = solidity*(1-a_ax)**2*cn/(np.sin(phi_rad)*np.sin(phi_rad))                 

###############        Correction Methods            ########################## 
                a_c       = 0.2
                if a_ax <= a_c:                     
                    c_T    = 4*F*a_ax*(1-a_ax)                        
                else:
                    K      = 4*F*np.sin(phi_rad)*np.sin(phi_rad)/(solidity*cn)
                    a_ax   = 1/2*(2+K*(1-2*a_c)-np.sqrt((K*(1-2*a_c)+2)**2+4*(K*a_c**2-1)))
                    c_T    = solidity*(1-a_ax)**2*cn/(np.sin(phi_rad)*np.sin(phi_rad))

                a_rot          = 1/((4*F*(np.sin(phi_rad)*np.cos(phi_rad)))/(ct*solidity)-1) 
               
                count          += 1                        
          
                link = 'Airfoils/polars/' + str(Profile) + '/Montgomerie/' + str(Profile) + '_T1_Re'+str('{:.3f}'.format(((re)/10e5))) + '_M0.00_N9.0_360_M.dat'
                with open(link, 'r') as fp:
                    Alpha = [];    c_L    = [];       c_D    = [];       L2D    = []  
                    for i in range (15):
                        next(fp)            
                    data = fp.read()
                    for line in data.splitlines():
                        temp = line.split()
                        if temp:
                            alpha = eval(temp[0]);          cl = eval(temp[1]);         cd = eval(temp[2]);     l2d = cl/cd
                            Alpha.append(alpha);       c_L.append(cl);        c_D.append(cd);    L2D.append(l2d)
                        else: 
                            break                    
                    aoa = (int((phi_deg-twist)*10)/10)
                    if aoa < -10 or aoa > 20:
                        aoa = int(aoa)
                    position        = Alpha.index(aoa)
                    alpha      = Alpha[position];      cl = c_L[position];      cd = c_D[position];      l2d = L2D[position]
        
                    delta_Re   = abs(Re-re)
                    if count > 2000:
                        print('Iteration loop exceeded', a_ax)
                        break  
                    if delta_Re <base:
                        Re_conv    = 1
            print(tsr, Re, cl, cd)
        print('new section', radius)
            #print(chord, radius, v_app, Re)
    ####################################################                       
        c_P            = 4*a_ax*(1-a_ax)**2
        d_c_P          = 2*(1-a_ax)*a_rot *TSR**2*(radius/Radius)**4        
         
        dT             = 4*F*a_ax*(1-a_ax)*rho*v_1**2*np.pi*radius*dr
        dQ             = 4*F*a_rot*(1-a_ax)*v_1*Omega*radius**2*np.pi*radius*dr*rho
        dP             = 2*np.pi*dA*v_1**3*a_ax*(1-a_ax)**2#dQ *Omega#0.5*rho*v_1**3*radius**2*np.pi*c_P
        
        dT_MT          = 4*F*a_ax*(1-a_ax)*rho*v_1**2*np.pi*radius*dr
        dQ_MT          = 4*F*a_rot*(1-a_ax)*v_1*Omega*radius**2*np.pi*radius*dr*rho

        dT_BET         = B*cn *0.5*rho*chord*(v_app**2)*dr
        dQ_BET         = B*ct *0.5*rho*chord*(v_app**2)*radius*dr            
         
        max_alpha_list.append(alpha);     max_cl_list.append(cl);           max_cd_list.append(cd);       max_l2d_list.append(l2d)
        Re_list.append(Re);               re_list.append(re)               
        phi_deg_list.append(phi_deg);     twist_list.append(twist)   #alpha_list.append(alpha);          
        a_ax_list.append(a_ax);           a_rot_list.append(a_rot);         TSR_list.append(TSR);         solidity_list.append(solidity)

        f_tip_list.append(f_tip);         f_hub_list.append(f_hub);         F_list.append(F)
        c_T_list.append(c_T);             c_P_list.append(d_c_P)
        dT_list.append(dT);               dQ_list.append(dQ);               P_list.append(dP)        
        v_app_list.append(v_app);         v_rot_list.append(v_rot)

    twist_list_smoothed = []
    #print(len(twist_list))
    for t in range(len(twist_list)):
        twist_smoothed = sum(twist_list[t:t+5])/(5)
        twist_list_smoothed.append(twist_smoothed)

    #print('{:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}'.format('Section', "Radius",'Twist','Phi',"Chord", 'v_app', 'Re', 're','alpha', 'c_L', 'c_D', 'L2D', "a","a'",'f_hub', 'f_tip', 'F', 'solidity'))         

    

    M_tot = sum(dQ_list)
    P_tot = sum(P_list)*efficiency#
    T_tot = sum(dT_list)
    F_av = sum(F_list)/len(F_list)


calculation()