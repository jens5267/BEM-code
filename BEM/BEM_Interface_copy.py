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

try:
    import IPython
    shell = IPython.get_ipython()
    shell.enable_matplotlib(gui='qt')
except:
    print('failed') 
    
#########
###############          entries          #####################################
    def BEM(self):
        self.rotor_diameter         = 0.2                                       # rotor diameter
        self.hub_diameter           = 0.056                                     # hub diameter       
        self.rotations              = 3000                                      # rotations 
        self.v_des                  = 15                                        # wind speed
        self.Blades                 = 3                                         # blades 
        self.sections               = 15                                        # sections
        self.efficiency             = 80                                        # efficiency
        self.airfoil                = 'NACA_2418'                               # airfoil
###############################################################################
###############          functions           ##################################
        self.aerodynamics()
        self.fix_blades()
###############################################################################
    def aerodynamics(self):        
        self.Profile            = self.airfoil
        self.Re_init            = 100000
        try:
            self.Re
        except AttributeError:                
            link = 'Airfoils/polars/' + str(self.Profile) + '/Montgomerie/' + str(self.Profile) + '_T1_Re'+str('{:.3f}'.format(((self.Re_init)/10e5))) + '_M0.00_N9.0_360_M.dat'
        else:
            self.base = 5000
            self.Re = int(self.base * round(float(self.Re)/self.base))

            link = 'Airfoils/polars/' + str(self.Profile) + '/Montgomerie/' + str(self.Profile) + '_T1_Re'+str('{:.3f}'.format(((self.Re)/10e5))) + '_M0.00_N9.0_360_M.dat'

        with open(link, 'r') as fp:
            for i in range (15):
                next(fp)            
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
                
        L2D = [ x/c_D[i] for i,x in enumerate(c_L) ]       
        alpha_pos = L2D.index(max(L2D[10:]))     
        self.alpha = alpha[alpha_pos] 
        self.cl = c_L[alpha_pos]
        self.cd = c_D[alpha_pos]           

    def calculation(self):
        '''
        Input:  rotor diameter, hub diameter, rotations, design wind velocity, Airfoil, amount of sections, efficiency, if desired: amount of blades
        Requirements: txt. data of the aerodynamic values of the used airfoils in the same folder
        Output: BEM code calculated parameters: local radii, chord-/ twist-distribution, thrust, torque, power (-coefficients)   
        '''
        self.aerodynamics()
        
        self.frequency          = self.rotations /60
        self.Omega              = 2*np.pi*self.frequency
        
        self.nu                 = 15.11*10**-6
        self.rho                = 1.225

        self.Radius             = self.rotor_diameter/2
        self.hub_radius         = self.hub_diameter/2
        self.d_r                = (self.Radius- self.hub_radius)/(self.sections-1)
        
        self.TSR                = self.Omega * self.Radius / self.v_des
        self.v_tip              = self.Omega*self.Radius 

        self.A                  = np.pi * (self.Radius**2 -  self.hub_radius**2)
        self.P_w                = 0.5* self.rho * self.v_des**3 * self.A
        self.P_t                = self.P_w * 16/27
        self.Method             = 'mixed'#'Glauert' #'straight Line' 'Wilson-Walker' 'mixed'

        #################################################################################
        print('{:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}'.format("a'",'a',"Radius","ax. ind.","chord", 'c/R', 'Twist', 'Re', 'solidity','c_P', 'F', 'f_tip', 'f_tip')) 
        
        self.a_ax_list          = []
        self.a_rot_list         = []
        
        self.chord_list         = []
        self.Re_list            = []
        self.twist_list         = []
        self.v_app_list         = []

        self.F_list             = []
        
        self.c_T_list           = []
        self.c_P_list           = []
        
        self.T_list             = []
        self.M_list             = []
        self.P_list             = []
        
        self.radius_list        = []


        self.a_ax_init       = 0.0 
        self.a_rot_init      = 0.001
        self.delta_a         = 1.0
        for i in range(self.sections-1):                                                                       
            self.radius          = self.hub_radius + (i*self.d_r+(i+1)*self.d_r)/2
            self.radius_list.append(self.radius)
            
            while abs(self.delta_a)> 10e-8:#self.a_ax_init < self.a_ax:#abs(a_ax - a_ax_init)>0:#a_ax > a_ax_init:#
                
                self.v_rot          = self.radius*2*np.pi*self.frequency
                self.v_app          = np.sqrt((self.v_des*(1-self.a_ax))**2+(self.v_rot*(1+self.a_rot))**2) 
                
                self.TSR_local      = self.TSR * self.radius / self.Radius
        
                self.phi_rad        = np.arctan(((1-self.a_ax)*self.v_des)/((1+self.a_rot)*(self.Omega*self.radius)))
                self.phi_deg        = np.rad2deg(self.phi_rad)
                
                self.twist          = self.phi_deg-self.alpha
###############        Prandtl Correction             ##########################                 
                self.f_tip          = 2/np.pi*(np.arccos(np.exp(-(self.B/2 *((self.Radius-self.radius)/(self.radius * np.sin(self.phi_rad)))))))
                self.f_hub          = 2/np.pi*(np.arccos(np.exp(-(self.B/2 *((self.radius-self.hub_radius)/(self.radius * np.sin(self.phi_rad)))))))
                self.F              = self.f_tip*self.f_hub

                self.chord          = 8*np.pi*self.radius*np.sin(self.phi_rad)/(3*self.B*self.cl*self.TSR_local)
                
                self.d_L            = 0.5*self.rho*self.chord*self.v_app**2*self.cl*self.d_r#*radius
                self.d_D            = 0.5*self.rho*self.chord*self.v_app**2*self.cd*self.d_r#*radius
                
                self.cn             = self.cl*np.cos(self.phi_rad) + self.cd*np.sin(self.phi_rad)
                self.cr             = self.cl*np.sin(self.phi_rad) - self.cd*np.cos(self.phi_rad)        
                
                self.a_rot          = 1/((4*self.F*(np.sin(self.phi_rad)*np.cos(self.phi_rad)))/(self.cr*self.solidity)-1)   
                
                self.d_T_MOM        = 4*np.pi*self.radius*self.rho*self.v_des**2*(1-self.a_ax)*self.a_ax*self.d_r*self.F
                self.d_M_MOM        = 4*np.pi*self.radius*self.rho*self.v_des*(self.Omega*self.radius**2*(1-self.a_ax))*self.a_rot*self.d_r*self.F
        
                self.d_T_BE         = self.B*self.cn *0.5*self.rho*self.chord*(self.v_app**2)*self.d_r
                self.d_M_BE         = self.B*self.cr *0.5*self.rho*self.chord*(self.v_app**2)*self.radius*self.d_r 
                
                self.solidity       = self.B*self.chord/(2*np.pi*self.radius)        

                if self.Method == 'Glauert':
                    if self.a_ax <= 0.4:
                        self.a_ax = 1/self.F*(0.143 + np.sqrt(0.0203-0.6427*(0.899-self.c_T)))
                        self.c_T = self.F*self.solidity*(1-self.a_ax)**2*(self.cl*self.cn/(np.sin(self.phi_rad)*np.sin(self.phi_rad)))
                        
                
                


                self.count       += 1
                self.delta_a     = self.a_ax_init - self.a_ax
               
#                self.B_EP = self.a_ax/(1+self.a_rot)*4*np.sin(self.phi_rad)
#                self.chord2 = self.B_EP *2*np.pi*self.v_des/(self.B*self.cl*self.Omega)

                self.Re = self.v_app*self.chord/self.nu
                
                self.epsilon = self.cd/self.cl
                self.eta = (1-self.epsilon*(1/np.tan(self.phi_rad)*(1/np.tan(self.phi_rad))/((1/np.tan(self.phi_rad)+self.epsilon))))
                
                self.c_P = 4*self.a_ax*(1-self.a_ax)**2
                self.d_c_P = 2*(1-self.a_ax)*self.a_rot *self.TSR**2*(self.radius/self.Radius)**4
                
                self.d_a = abs(self.a_ax_init-self.a_ax)

                if self.a_ax_init+self.a_ax >0:
                    self.a_ax_init = (self.a_ax_init+self.a_ax)/2
                
###############        Correction Methods            ##########################                
                self.c_T = (1-self.a_ax)**2*self.solidity*self.cn/(np.sin(self.phi_rad)*np.sin(self.phi_rad))

                self.T = 0.5*self.rho*self.c_T*np.pi*self.radius**2*self.v_des**2
                self.P = self.Omega*self.d_M_BE*self.efficiency
                

            self.a_ax_list.append(self.a_ax),       self.a_rot_list.append(self.a_rot)
            self.chord_list.append(self.chord)
            self.c_P_list.append(self.c_P)
            self.c_Thrust_list.append(self.c_T)
            self.F_list.append(self.F)
            self.M_list.append(self.d_M_BE)
            self.P_list.append(self.P)

            self.Re_list.append(self.Re)#,          self.Re2_list.append(self.rey2)
            #self.F_list.append(self.F)
            #self.Thrust_list.append(self.d_T_BE)
            self.Thrust_list.append(self.T)
            #self.T2_list.append(self.T2)
            self.twist_list.append(self.twist)
            self.v_app_list.append(self.v_app), self.T_list.append(self.d_T_BE)
            #self.T_list.append(self.d_T_MOM)
            print('{:>8,.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8,.0f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}'.format(self.a_rot, self.a_ax_init, self.radius, self.a_ax_init, self.chord, self.chord/self.Radius, self.twist, self.Re, self.solidity, self.c_P, self.F, self.f_tip, self.f_hub, self.phi_rad))
        self.M_tot = sum(self.M_list)
        self.P_tot = sum(self.P_list)*self.efficiency

###############################################################################
###############          print           ######################################
        print(50 * '--')
        print('{:<15s}: {:>10.3f} m | {:<15s}: {:>10.0f}   | {:<15s}: {:>10.4f} m'.format      ('Rot Diameter', self.rotor_diameter,'Sections', self.sections -1, 'delta r', self.d_r))
        print('{:<15s}: {:>10s}   | {:<15s}: {:>10.3f} m/s'.format          ('Profile', self.profile, 'des. wind speed', self.v_des))
        print('{:<15s}: {:>10.0f}   | {:<15s}: {:>10.3f} rpm'.format          ('Blades', self.B, 'Rotations', self.rotations))
        if self.v_tip > 80:
            print('{:<15s}: {:>10.3f} m/s -> above 80 m/s'.format   ('Tip velocity', self.v_tip, ))  
        else:
            print('{:<15s}: {:>10.3f} m/s -> below 80 m/s'.format   ('Tip velocity', self.v_tip))
        print('{:<15s}: {:>10.3f}   '.format       ('TSR', self.TSR, ))
        print('{:<15s}: {:>10.3f} W | {:<15s}: {:>10.3f} W | {:<15s}: {:>10.3f} W | {:<15s}: {:>10.3f} %'.format('Power wind', self.P_w, 'Power Betz', self.P_t, 'Power WT', self.P_tot, 'Efficiency η', 100*self.P_tot/self.P_w))
        print('{:<15s}: {:>10.3f}   | {:<15s}: {:>10.3f}   | {:<15s}: {:>10.3f} ° | {:<15s}: {:>10.3f} '.format('lift coeff.', self.cl, 'drag coeff.', self.cd, 'AoA', self.alpha, 'L2D', self.cl/self.cd))
        print('{:<15s}: {:>10.3f} Nm'.format('Torque WT', self.M_tot))
        print(50*'--')
###############################################################################
###############          plot           #######################################
    def plotting(self):   
        if self.plot.get() == 1:              
            
            plt.figure('plots', figsize=(10,5))
            plt.title(self.profile)
            ax1 = plt.subplot(221)
            plt.xlabel('Radius [m]')
            plt.ylabel('Twist [°]')
            plt.xlim(0,self.Radius*1.1)
            plt.ylim(0.0, 36.0)
            plt.minorticks_on()
            plt.grid(True,which='major')
            plt.grid(True,which='minor', linestyle='--')
            ax1.plot(self.radius_list,self.twist_list, color = 'k', marker = '.')
            
            ax2 = plt.subplot(222)
            plt.xlabel('Radius [m]')
            plt.ylabel('Chord [m]')
            plt.xlim(0,self.Radius*1.1)
            plt.ylim(0.0, 0.12)
            plt.minorticks_on()
            plt.grid(True,which='major')
            plt.grid(True,which='minor', linestyle='--')
            ax2.plot(self.radius_list,self.chord_list, color = 'k', marker = '.')
            
            ax3 = plt.subplot(223)
            plt.xlabel('Radius [m]')
            plt.ylabel('Re [-]')
            plt.xlim(0,self.Radius*1.1)
            plt.ylim(0,100000)
            plt.minorticks_on()
            plt.grid(True,which='major')
            plt.grid(True,which='minor', linestyle='--')
            ax3.plot(self.radius_list,self.Re_list, color = 'k', marker = '.')            
            
            ax4 = plt.subplot(224)
            plt.xlabel('Radius [m]')
            plt.ylabel('apparent wind velocity [m/s]')
            plt.xlim(0,self.Radius*1.1)
            plt.ylim(0,80)
            plt.minorticks_on()
            plt.grid(True,which='major')
            plt.grid(True,which='minor', linestyle='--')
            ax4.plot(self.radius_list,self.v_app_list, color = 'k', marker = '.')               
            
            plt.show()

###############################################################################
#self.F = 1
#chord       = 16*np.pi/(self.B) * radius/self.cl * (np.sin(1/3 * np.arctan(self.Radius/radius/self.TSR)))**2
#chord       = 8 * F*radius*np.pi*np.sin(phi_rad)/(self.B*self.cl)*((1-self.TSR_local*np.tan(phi_rad))/(np.tan(phi_rad+self.TSR_local)))
#self.chord       = 8*np.pi*self.radius/(self.B*self.cl)*(1-np.cos(self.phi_rad))
            
#                if self.Method == 'straight Line':
#                    self.c_T1 = 1.816
#                    self.c_T = self.c_T1 -4 *(np.sqrt(self.c_T1)-1)*(1-self.a_ax)
#                
#                if self.Method == 'Wilson-Walker':
#                    self.a_ac = 0.2
#                    self.K = 4*self.F*np.sin(self.phi_rad)*np.sin(self.phi_rad)/(self.solidity*self.cn)
#                    if self.a_ax <= self.a_ac:
#                        self.c_T = 4*self.F*self.a_ax*(1-self.a_ax)
#                    elif self.a_ax >0.2:
#                        self.a_ax = 1/2*(2+self.K*(1-2*self.a_ac))-np.sqrt((self.K*(1-2*self.a_ac)**2+4*(self.K*self.a_ax**2-1)))
#                        self.c_T = 4*self.F*(self.a_ac**2+self.a_ax*(1-2*self.a_ac))
#                if self.Method == 'this one':
#                    self.a_c = 0.2
#                    if self.a_ax <= self.a_c:
#                        self.K           = 4*self.F*np.sin(self.phi_rad)*np.sin(self.phi_rad)/(self.solidity*self.cn)
#                        self.a_ax       = 1/(self.K+1)
#                    else:
#                    #    self.a_ax    = 1/((4*self.F*(np.sin(self.phi_rad)*np.sin(self.phi_rad)))/(self.cn*self.solidity)+1)
#                        self.a_ax    = 1/2*(2+self.K*(1-2*self.a_c)-np.sqrt((self.K*(1-2*self.a_c)+2)**2+4*(self.K*self.a_c**2-1)))            
            
            
            #                self.a_c = 0.2
#                
#                if self.a_ax <= self.a_c:
#                    self.K          = 4*self.F*np.sin(self.phi_rad)*np.sin(self.phi_rad)/(self.solidity*self.cn)
#                    self.a_ax       = 1/(self.K+1)
#                else:
#                #    self.a_ax    = 1/((4*self.F*(np.sin(self.phi_rad)*np.sin(self.phi_rad)))/(self.cn*self.solidity)+1)
#                    self.a_ax    = 1/2*(2+self.K*(1-2*self.a_c)-np.sqrt((self.K*(1-2*self.a_c)+2)**2+4*(self.K*self.a_c**2-1)))