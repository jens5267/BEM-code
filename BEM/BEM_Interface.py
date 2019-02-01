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
from mpl_toolkits.mplot3d import Axes3D
from tkinter import *
import numpy as np
import matplotlib.pyplot as plt
import os
import json
import matplotlib.ticker as mtick


try:
    import IPython
    shell = IPython.get_ipython()
    shell.enable_matplotlib(gui='qt')
except:
    print('failed') 
    
class Bladedesign(object):
    def __init__(self, master):
        self.master = master
        master.title("Blade Design")
        
        window_size                 = root.geometry('400x270+100+50')
        window_title                = root.title('Blade design')
###############################################################################
###############          entries          #####################################
# rotor diameter
        self.rotor_diameter         = StringVar()
        self.label_rotor_diameter   = Label(root, text = 'Rotor Diameter [m]')
        self.label_rotor_diameter.grid(row = 0, column = 0, sticky = E)        
        self.rotor_diameter_entry   = Entry(root, textvariable = self.rotor_diameter)
        self.rotor_diameter_entry.delete(0, END)
        self.rotor_diameter_entry.insert(0, "0.4")
        self.rotor_diameter_entry.grid(   row = 0, column = 1, sticky = W)
       
###############################################################################        
# hub diameter
        self.hub_diameter           = StringVar()
        self.label_hub_diameter     = Label(root, text = 'Hub Diameter [m]')
        self.label_hub_diameter.grid(row = 1, column = 0, sticky = E)
        self.hub_diameter_entry     = Entry(root, textvariable = self.hub_diameter)
        self.hub_diameter_entry.delete(0, END)
        self.hub_diameter_entry.insert(0, "0.052")        
        self.hub_diameter_entry.grid(     row = 1, column = 1, sticky = W)
###############################################################################
# rotations        
        self.rotations              = StringVar()
        self.label_rotations        = Label(root, text = 'Rotations [rpm]:')
        self.label_rotations.grid(  row = 2, column = 0, sticky = E)
        self.rotations_entry        = Entry(root, textvariable = self.rotations)
        self.rotations_entry.delete(0, END)
        self.rotations_entry.insert(0, "3000")        
        self.rotations_entry.grid(  row = 2, column = 1, sticky = W)
###############################################################################
# wind speed
        self.v_1                  = StringVar()
        self.label_v_1            = Label(root, text = 'Design Wind Speed [m/s]:')
        self.label_v_1.grid(      row = 3, column = 0, sticky = E)
        self.v_1_entry                  = Entry(root, textvariable = self.v_1)
        self.v_1_entry.delete(0, END)
        self.v_1_entry.insert(0, "15")        
        self.v_1_entry.grid(            row = 3, column = 1, sticky = W)        

###############################################################################
# blades
        self.fixed_blades           = IntVar()
        self.label_blades           = Label(self.master, text="Number of blades:")
        self.label_blades.grid(row  = 5, column = 0, sticky = E)
        self.blades_entry           = Entry(self.master)
        self.blades_entry.grid(row  = 5, column = 1, sticky = W)
        self.blades_entry.config(state=DISABLED)
        self.plot                   = IntVar()
        self.save                   = IntVar()
###############################################################################
# efficiency
        self.efficiency             = StringVar()
        self.label_efficiency       = Label(master, text = 'Efficiency: [%]')
        self.label_efficiency.grid( row = 6, column = 0, sticky = E)     
        self.efficiency_entry             = Entry(master, textvariable = self.efficiency)
        self.efficiency_entry.delete(0, END)
        self.efficiency_entry.insert(0, "80")
        self.efficiency_entry.grid(       row = 6, column = 1, sticky = W) 
###############################################################################
#sections
        self.sections               = StringVar()
        self.label_sections         = Label(root, text = 'Sections:')
        self.label_sections.grid(   row = 8, column = 0, sticky = E)
        self.sections_entry               = Spinbox(root, from_ = 1, to= 2000)
        self.sections_entry.delete(0, END)
        self.sections_entry.insert(0, "10")        
        self.sections_entry.grid(         row = 8, column = 1, sticky = W) 
###############################################################################        
###############          3D design          ###################################
        self.status                 = IntVar()        
###############################################################################        
###############          drop down list          ##############################
# airfoil
        self.airfoil                = StringVar()
        self.label_airfoil          = Label(master, text = 'Airfoil:')
        self.label_airfoil.grid(    row = 4, column = 0, sticky = E)
        self.airfoils = ('NACA_2408', 'NACA_2412','NACA_2415','NACA_2418','NACA_4412','NACA_4415','NACA_4418', 'S822','S833', 'S1210', 'S2091', 'S4061', 'S4180', 'S4320', 'SD7034', 'SG6040', 'SG6041', 'SG6042', 'SG6043')
        self.airfoil.set('NACA_2418')
        self.selection_entry = OptionMenu(master, self.airfoil, *self.airfoils)
        self.selection_entry.grid(        row = 4, column = 1, sticky = W)        
###############################################################################  
        self.Method                = StringVar()
        self.label_Method          = Label(master, text = 'Correction method:')
        self.label_Method.grid(    row = 11, column = 0, sticky = E)
        self.Methods = ('Wilson-Walker', 'straight line', 'Glauert', 'Buhl')
        self.Method.set('Wilson-Walker')
        self.Methods_entry = OptionMenu(master, self.Method, *self.Methods)
        self.Methods_entry.grid(        row = 11, column = 1, sticky = W)        
###############################################################################         
###############          buttons          #####################################
        self.button_calc            = Button(root, text="Calculate", command=self.calculation)
        self.button_calc.config(    height = 1, width = 10)
        self.button_calc.grid(      row = 9, column = 1, sticky = W)
###############################################################################
        self.button_close           = Button(master, text="Close", command=master.destroy)
        self.button_close.config(   height = 1, width = 10)
        self.button_close.grid(     row = 10, column = 1, sticky = W)
###############################################################################
        self.button_fix             = Checkbutton(self.master,text='Fixed Number', variable=self.fixed_blades, command=self.fix_blades)
        self.button_fix.grid(       row = 5, column = 3, sticky = W)       
##############################################################################
        self.button_plot           = Checkbutton(self.master, text= 'Plot', variable=self.plot, command=self.plotting)
        self.button_plot.grid(     row = 9, column = 3, sticky = W)  
##############################################################################
        self.button_save           = Checkbutton(self.master, text= 'Save', variable=self.save, command=self.saving)
        self.button_save.grid(     row = 10, column = 3, sticky = W)          
###############################################################################
        self.button_3D             = Checkbutton(self.master,text='Show Blade', variable=self.status, command=self.shaping_3D)
        self.button_3D.grid(        row = 11, column = 3, sticky = W)                
###############################################################################      
###############          functions           ##################################
        self.aerodynamics()
        self.fix_blades()
###############################################################################
    def aerodynamics(self):        
        self.Profile            = self.airfoil.get()
        self.Re                 = 60000
#        link = 'Airfoils/polars/'+ str(self.Profile) +'/' + self.Profile + '_'+ str(int(self.Re/1000)) + 'k.txt'
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
###############################################################################
    def fix_blades(self):
        if self.fixed_blades.get() == 1:          
            self.blades_entry.config(state = NORMAL)
                
        elif self.fixed_blades.get() == 0:       
            self.blades_entry.config(state = DISABLED)
###############################################################################
    def calculation(self):
        '''
        Input:  rotor diameter, hub diameter, rotations, design wind velocity, Airfoil, amount of sections, efficiency, if desired: amount of blades
        Requirements: txt. data of the aerodynamic values of the used airfoils in the same folder
        Output: BEM code calculated parameters: local radii, chord-/ twist-distribution, thrust, torque, power (-coefficients)   
        '''
        self.aerodynamics()
 
        self.rotor_diameter     = float(self.rotor_diameter_entry.get())
        self.hub_diameter       = float(self.hub_diameter_entry.get())
        self.rotations          = float(self.rotations_entry.get())
        self.v_1                = float(self.v_1_entry.get())
        self.profile            = self.airfoil.get()

        self.method             = self.Method.get()

        self.sections           = int(self.sections_entry.get())+1
        self.efficiency         = float(self.efficiency_entry.get())/100
        
        self.frequency          = self.rotations /60
        self.Omega              = 2*np.pi*self.frequency
        
        self.nu                 = 15.11*10**-6
        self.rho                = 1.225

        self.Radius             = self.rotor_diameter/2
        self.hub_radius         = self.hub_diameter/2
        self.dr                = (self.Radius- self.hub_radius)/(self.sections-1)
        
        self.TSR                = self.Omega * self.Radius / self.v_1
        self.v_tip              = self.Omega*self.Radius 

        self.A                  = np.pi * (self.Radius**2 -  self.hub_radius**2)
        self.P_w                = 0.5* self.rho * self.v_1**3 * self.A
        self.P_t                = self.P_w * 16/27

        if self.fixed_blades.get() == 0:
            if self.TSR < 1:
                self.B = 5 
            elif self.TSR >= 1 and self.TSR < 4:
                self.B = 4
            elif self.TSR >= 4 and self.TSR < 8:
                self.B = 3
            elif self.TSR >= 8 :
                self.B = 2
        else:
            self.B                  = int(self.blades_entry.get())
#################################################################################
        print('{:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}| {:^8}'.format("a","a'",'Section', "Radius","chord", 'c/R', 'Phi','Twist', 'Re', 'solidity','c_P', 'dT', 'dQ', 'F', 'AoA')) 
        
        self.a_ax_list          = []
        self.a_rot_list         = []
        self.chord_list         = []
        self.c_P_list           = []
        self.c_T_list           = []
        self.F_list             = []
        
        self.M_list             = []
        self.P_list             = []
        self.T_list             = []
        
        self.radius_list        = []
        self.Re_list            = []
        self.Re2_list           = []        

        self.T_list             = []
        self.T2_list            = []
        self.twist_list         = []
        self.v_app_list         = []
        self.convergence_error = 10e-10
        for i in range(self.sections-1):                                                                       
            self.radius          = self.hub_radius + (i*self.dr+(i+1)*self.dr)/2

            self.a_ax_init       = 0.0 # can be any value below 0.3, must be smaller than a_ax
            self.a_rot_init         = 0
            self.delta_a         = 1.0
            self.delta_ap       = 1
            self.a_ax            = 0.3
            self.a_rot           = 0.001
            self.solidity        = 0.01
            
            self.chord           = 0
            self.count           = 0
            self.new             = 0
            #while abs(self.delta_a)> 10e-8:#self.a_ax_init < self.a_ax:#abs(a_ax - a_ax_init)>0:#a_ax > a_ax_init:#
            while self.new != 1:   
                self.v_rot          = self.radius*2*np.pi*self.frequency
                self.v_app          = np.sqrt((self.v_1*(1-self.a_ax))**2+(self.v_rot*(1+self.a_rot))**2) 
                
                self.TSR_local      = self.TSR * self.radius / self.Radius
        
                self.phi_rad        = abs(np.arctan(((1-self.a_ax)*self.v_1)/((1+self.a_rot)*(self.Omega*self.radius))))
                self.phi_deg        = np.rad2deg(self.phi_rad)
                
                self.twist          = self.phi_deg-self.alpha
###############        Prandtl Correction             ##########################                 
                self.f_tip              = 2/np.pi*(np.arccos(np.exp(-(self.B/2 *((self.Radius-self.radius)/(self.radius * np.sin(self.phi_rad)))))))
                self.f_hub              = 2/np.pi*(np.arccos(np.exp(-(self.B/2 *((self.radius-self.hub_radius)/(self.radius * np.sin(self.phi_rad)))))))
                self.F                  = self.f_tip*self.f_hub
                #chord       = 16*np.pi/(self.B) * radius/self.cl * (np.sin(1/3 * np.arctan(self.Radius/radius/self.TSR)))**2
                #chord       = 8 * F*radius*np.pi*np.sin(phi_rad)/(self.B*self.cl)*((1-self.TSR_local*np.tan(phi_rad))/(np.tan(phi_rad+self.TSR_local)))
                #self.chord       = 8*np.pi*self.radius/(self.B*self.cl)*(1-np.cos(self.phi_rad))
                self.chord          = abs(8*np.pi*self.radius*np.sin(self.phi_rad)/(3*self.B*self.cl*self.TSR_local))
                self.d_L            = 0.5*self.rho*self.chord*self.v_app**2*self.cl*self.dr#*radius
                self.d_D            = 0.5*self.rho*self.chord*self.v_app**2*self.cd*self.dr#*radius
                
                self.cn             = self.cl*np.cos(self.phi_rad) + self.cd*np.sin(self.phi_rad)
                self.ct             = self.cl*np.sin(self.phi_rad) - self.cd*np.cos(self.phi_rad)        
                
                self.a_rot          = 1/((4*self.F*(np.sin(self.phi_rad)*np.cos(self.phi_rad)))/(self.ct*self.solidity)-1) 
                self.a_ax           = 1/((4*self.F*(np.sin(self.phi_rad)*np.sin(self.phi_rad)))/(self.cn*self.solidity)+1)
                
                
                self.solidity       = self.B*self.chord/(2*np.pi*self.radius)  

###############        Correction Methods            ########################## 
                if self.method == 'Glauert':                    
                    if self.a_ax <= 0.4:
                        self.c_T = 4*self.a_ax*(1-self.a_ax)
                        print('Glauert <0.4')
                    else:
                        print('Glauert > 0.4')
                        self.a_ax   = 1/self.F*(0.143 + np.sqrt(0.0203-0.6427*(0.899-self.c_T)))  
                        self.c_T    = self.solidity*(1-self.a_ax)**2*(self.cn/(np.sin(self.phi_rad)*np.sin(self.phi_rad)))

                if self.method == 'straight line':
                    self.c_T1   = 1.816
                    if self.a_ax >= 0.4:
                        self.c_T    = self.c_T1 -4 *(np.sqrt(self.c_T1)-1)*(1-self.a_ax)
                    else:
                        self.c_T = 4*self.F*self.a_ax*(1-self.a_ax)
                
                if self.method == 'Wilson-Walker':
                    self.a_c       = 0.2
                    if self.a_ax <=0.2:                     
                        self.c_T    = 4*self.F*self.a_ax*(1-self.a_ax)                        
                    else:
                        self.K      = 4*self.F*np.sin(self.phi_rad)*np.sin(self.phi_rad)/(self.solidity*self.cn)
                        self.a_ax   = 1/2*(2+self.K*(1-2*self.a_c)-np.sqrt((self.K*(1-2*self.a_c)+2)**2+4*(self.K*self.a_c**2-1)))
                        self.c_T    = self.solidity*(1-self.a_ax)**2*self.cn/(np.sin(self.phi_rad)*np.sin(self.phi_rad))

                if self.method == 'Buhl':
                    self.c_T1   = 1.816
                    if self.a_ax >= 0.4:
                        self.c_T    = 8/9+(4*self.F-20*self.c_T1/9)*self.a_ax +(25*self.c_T1/9-4*self.F)* self.a_ax**2
                        self.a_ax = (18*self.F-20-3*np.sqrt(self.c_T*(50-36*self.F)+12*self.F*(3*self.F-4)))/(36*self.F-50)
                    else:
                        self.c_T = 4*self.F*self.a_ax*(1-self.a_ax)

                self.Re = self.v_app*self.chord/self.nu
                self.base = 5000
                self.Re = int(self.base * round(float(self.Re)/self.base))
                
                self.epsilon = self.cd/self.cl
                self.eta = (1-self.epsilon*(1/np.tan(self.phi_rad)*(1/np.tan(self.phi_rad))/((1/np.tan(self.phi_rad)+self.epsilon))))
                
                self.c_P = 4*self.a_ax*(1-self.a_ax)**2
                self.d_c_P = 2*(1-self.a_ax)*self.a_rot *self.TSR**2*(self.radius/self.Radius)**4                

                self.delta_a        = abs(self.a_ax_init - self.a_ax)
                self.delta_ap       = abs(self.a_rot_init- self.a_rot)
                
                self.count          += 1
                self.a_ax_init = (self.a_ax_init+self.a_ax)/2
                self.a_rot_init = (self.a_rot_init+self.a_rot)/2
                if  self.delta_a < self.convergence_error:
                    
                    self.dT             = 4*self.F*self.a_ax*(1-self.a_ax)*self.rho*self.v_1**2*np.pi*self.radius*self.dr
                    self.dQ             = 4*self.F*self.a_rot*(1-self.a_ax)*self.v_1*self.Omega*self.radius**2*np.pi*self.radius*self.dr*self.rho
                    
                    self.dT_MT          = 4*self.F*self.a_ax*(1-self.a_ax)*self.rho*self.v_1**2*np.pi*self.radius*self.dr
                    self.dQ_MT          = 4*self.F*self.a_rot*(1-self.a_ax)*self.v_1*self.Omega*self.radius**2*np.pi*self.radius*self.dr*self.rho
            
                    self.dT_BET         = self.B*self.cn *0.5*self.rho*self.chord*(self.v_app**2)*self.dr#self.v_1**2*(1-self.a_ax)**2/(np.sin(self.phi_rad)*np.sin(self.phi_rad))
                    self.dT_BET = self.solidity*np.pi*self.rho*self.Omega**2*self.radius**2*(1+self.a_rot)**2*self.cn*self.radius*self.dr/(np.cos(self.phi_rad)*np.cos(self.phi_rad))
                    self.dQ_BET         = self.B*self.ct *0.5*self.rho*self.chord*(self.v_app**2)*self.radius*self.dr# self.v_app = np.sqrt(self.v_1*(1-self.a_ax)/np.sin(self.phi_rad) *self.Omega*self.radius*(1+self.a_rot)/(np.cos(self.phi_rad)))
                    
                    self.new = 1
                elif self.count > 1000:
                    break
                
            self.a_ax_list.append(self.a_ax),           self.a_rot_list.append(self.a_rot)
            self.chord_list.append(self.chord),         self.twist_list.append(self.twist);             self.radius_list.append(self.radius)
            self.c_T_list.append(self.c_T);             self.c_P_list.append(self.c_P)
            self.F_list.append(self.F)
            self.T_list.append(self.dT);                self.M_list.append(self.dQ)
#            self.P_list.append(self.P)

            self.Re_list.append(self.Re)#,          

            self.v_app_list.append(self.v_app)
            #print(self.dT_BET, self.dT_MT, self.dT_BET-self.dT_MT)            
            print('{:>8.3f}| {:>8.3f}| {:>8.0f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8,.0f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}| {:>8.3f}'.format(self.a_ax, self.a_rot, i+1, self.radius, self.chord, self.chord/self.Radius, self.phi_deg, self.twist, self.Re, self.solidity, self.c_P, self.dT, self.dQ, self.F, self.alpha))
        self.M_tot = sum(self.M_list)
        self.P_tot = sum(self.P_list)*self.efficiency
        self.T_tot = sum(self.T_list)
        self.F_av = sum(self.F_list)/len(self.F_list)
        #print(self.F_av)
        #print(sum(self.P_list)*self.F_av)

    #############################################################
        self.turb_name = '%skrpm %sm/s %s m %s %sB'%(self.rotations/1000, self.v_1, self.rotor_diameter, self.profile, self.B)
        self.plotting()
###############################################################################
###############          print           ######################################
        print(50 * '--')
        print('{:<15s}: {:>10.3f} m | {:<15s}: {:>10.0f}   | {:<15s}: {:>10.4f} m'.format      ('Rot Diameter', self.rotor_diameter,'Sections', self.sections -1, 'delta r', self.dr))
        print('{:<15s}: {:>10s}   | {:<15s}: {:>10.3f} m/s'.format          ('Profile', self.profile, 'des. wind speed', self.v_1))
        print('{:<15s}: {:>10.0f}   | {:<15s}: {:>10.3f} rpm'.format          ('Blades', self.B, 'Rotations', self.rotations))
        if self.v_tip > 80:
            print('{:<15s}: {:>10.3f} m/s -> above 80 m/s'.format   ('Tip velocity', self.v_tip, ))  
        else:
            print('{:<15s}: {:>10.3f} m/s -> below 80 m/s'.format   ('Tip velocity', self.v_tip))
        print('{:<15s}: {:>10.3f}   '.format       ('TSR', self.TSR, ))
        print('{:<15s}: {:>10.3f} W | {:<15s}: {:>10.3f} W | {:<15s}: {:>10.3f} W | {:<15s}: {:>10.3f} %'.format('Power wind', self.P_w, 'Power Betz', self.P_t, 'Power WT', self.P_tot, 'Efficiency η', 100*self.P_tot/self.P_w))
#        print('{:<15s}: {:>10.3f}   | {:<15s}: {:>10.3f}   | {:<15s}: {:>10.3f} ° | {:<15s}: {:>10.3f} '.format('lift coeff.', self.cl, 'drag coeff.', self.cd, 'AoA', self.alpha, 'L2D', self.cl/self.cd))
        print('{:<15s}: {:>10.3f} Nm| {:<15s}: {:>10.3f} N'.format('Torque WT', self.M_tot, 'Thrust WT', self.T_tot))  
        print(50*'--')
###############################################################################
###############          plot           #######################################
    def plotting(self):   
        repetition = 1
        if self.plot.get() == 1:    
            repetition += 1
            
            plt.figure(int(repetition), figsize=(10,5))
            plt.title(self.profile)
            ax1 = plt.subplot(321)
            plt.xlabel('Radius [m]')
            plt.ylabel('Twist [°]')
            plt.xlim(0,self.Radius*1.1)
            plt.ylim(-6.0, 36.0)
            plt.yticks([0.00,6.00,12.00,18.00,24.00,30.00, 36.00],['0.00','6.00','12.00','18.00','24.00','30.00','36.00'])
            plt.minorticks_on()
            plt.grid(True,which='major')
            plt.grid(True,which='minor', linestyle='--')
            ax1.plot(self.radius_list,self.twist_list, color = 'k', marker = '.')
            
            ax2 = plt.subplot(322)
            plt.xlabel('Radius [m]')
            plt.ylabel('Chord [m]')
            plt.xlim(0,self.Radius*1.1)
            plt.ylim(0.0, 0.12)
            plt.minorticks_on()
            plt.grid(True,which='major')
            plt.grid(True,which='minor', linestyle='--')
            ax2.plot(self.radius_list,self.chord_list, color = 'k', marker = '.')
            
            ax3 = plt.subplot(323)
            plt.xlabel('Radius [m]')
            plt.ylabel('Re [-]')
            plt.xlim(0,self.Radius*1.1)
            plt.ylim(0,100000)
            plt.yticks([0,25000,50000,75000,100000,125000, 150000], [0,'25,000', '50,000', '75,000', '100,000', '125,000', '150,000'])
            plt.minorticks_on()
            plt.grid(True,which='major')
            plt.grid(True,which='minor', linestyle='--')
            ax3.plot(self.radius_list,self.Re_list, color = 'k', marker = '.')            
            
            ax4 = plt.subplot(324)
            plt.xlabel('Radius [m]')
            plt.ylabel('apparent wind velocity [m/s]')
            plt.xlim(0,self.Radius*1.1)
            plt.ylim(0,80)
            plt.minorticks_on()
            plt.grid(True,which='major')
            plt.grid(True,which='minor', linestyle='--')
            ax4.plot(self.radius_list,self.v_app_list, color = 'k', marker = '.')            
            
            ax5 = plt.subplot(325)
            plt.xlabel('Radius [m]')
            plt.ylabel("a")
            plt.xlim(0,self.Radius*1.1)
            plt.ylim(0,0.5)
            plt.minorticks_on()
            plt.grid(True,which='major')
            plt.grid(True,which='minor', linestyle='--')
            ax5.plot(self.radius_list,self.a_ax_list, color = 'k', marker = '.')              
            
            ax6 = plt.subplot(326)
            plt.xlabel('Radius [m]')
            plt.ylabel("a'")
            plt.xlim(0,self.Radius*1.1)
            plt.ylim(0,0.5)
            plt.minorticks_on()
            plt.grid(True,which='major')
            plt.grid(True,which='minor', linestyle='--')
            ax6.plot(self.radius_list,self.a_rot_list, color = 'k', marker = '.')               
            plt.show()
            #self.saving()            
#        if self.plot.get() == 0:
#            plt.close("plots")        

###############################################################################
###############          save           #######################################
    def saving(self):
        if self.save.get() == 1:
            
            
            self.turb_save = '%skrpm_%smps_%scm_%s_%sB'%(int(self.rotations/1000), int(self.v_1), int(self.rotor_diameter*100), self.profile, self.B)
            plt.savefig('plots/'+ self.turb_save + '.png', dpi=1000)
            file_name = '%sB_%skrpm_%sm-s_%s_cm_%s'%(self.B, self.rotations/1000, self.v_1, self.rotor_diameter, self.profile)
            
            with open('AWT/' + file_name + '.txt','w') as f:
                print('Blade Design', file = f)

                t1 = 'Radial Position [m]'
                t2 = 'Chord Length [m]'
                t3 = 'Twist [deg]'
                t4 = 'Pitch Axis Offset [m]'
                t5 = 'Thread Axis in [% chord]'
                t6 = 'Airfoil Name'
                t7 = '360 Polar Name'
                print('{:>26s}{:>26s}{:>26s}{:>26s}{:>26s}{:>26s}{:>26s}'.format(t1,t2,t3,t4,t5,t6,t7), file = f)
                print(182*'-', file = f)
                for i in range(self.sections-1):
                    self.rey2 = str('T1_Re' + str(self.Re) +'_M0.00_N9.0 360 M')
                    self.Re2_list.append(self.rey2)
                    print('{:26.5e}{:26.5e}{:26.5e}{:26.5e}{:26.5e}{:>26s}{:>28s}'.format(self.radius_list[i], self.chord_list[i], self.twist_list[i], 0, 0.25, self.profile, self.Re2_list[i]), file = f)

###############################################################################
###############          3D           #########################################
    def shaping_3D(self):
        if self.status.get() == 1:
            if self.profile[:4] =='NACA':
                self.airfoil_number = self.airfoil.get()[5:]
                shell.enable_matplotlib(gui='qt') 
                fig                     = plt.figure('3D', figsize=(10,5))
                ax                      = fig.add_subplot(111, projection = '3d')
       
                self.m = int(self.airfoil_number[0])/100
                self.p = int(self.airfoil_number[1])/10
                self.t = int(self.airfoil_number[2])/10 + int(self.airfoil_number[3])/100
                self.c = 1.0
                self.x = np.linspace(0,1,200)
    
                self.term1 =  0.2969 * (np.sqrt(self.x/self.c))
                self.term2 = -0.1260 * (self.x/self.c)
                self.term3 = -0.3516 * (self.x/self.c)**2
                self.term4 =  0.2843 * (self.x/self.c)**3
                self.term5 = -0.1015 * (self.x/self.c)**4
                
                self.y_t= (5 * self.t * self.c * (self.term1 + self.term2 + self.term3 + self.term4 + self.term5))
        
                self.x_u_list = []                                                     # upper x-values
                self.x_l_list = []                                                     # lower x-values
                self.y_u_list = []                                                     # upper y-values
                self.y_l_list = []                                                     # lower y-values
        
                for i in range(len(self.x)):
                    if 0 < self.x[i] and self.x[i]<(self.c*self.p):
                        self.y_c = (self.m/self.p**2)*(2*self.p*(self.x[i]/self.c)-(self.x[i]/self.c)**2)
                        self.d_xy = (2* self.m/(self.p**2))*(self.p-(self.x[i]/self.c))
                    else:
                        self.y_c = (self.m/(1-self.p)**2)*((1-2*self.p)+2*self.p*(self.x[i]/self.c)-(self.x[i]/self.c)**2)
                        self.d_xy = (2* self.m/((1-self.p)**2))*(self.p-(self.x[i]/self.c))
        
                    self.angle = (np.arctan(self.d_xy))
                    self.x_u = self.x[i] - self.y_t[i] * np.sin(self.angle)
                    self.x_l = self.x[i] + self.y_t[i] * np.sin(self.angle)
                    
                    self.x_u_list.append(self.x_u) 
                    self.x_l_list.append(self.x_l) 
                    
                    self.y_u = self.y_c + self.y_t[i] * np.cos(self.angle)
                    self.y_l = self.y_c - self.y_t[i] * np.cos(self.angle)
                    
                    self.y_u_list.append(self.y_u) 
                    self.y_l_list.append(self.y_l) 
                    
                self.x_u_list = self.x_u_list[::-1] 
                self.y_u_list = self.y_u_list[::-1] 
                
                self.x_list = self.x_u_list+self.x_l_list        
                self.y_list = self.y_u_list+self.y_l_list    
##############################

            else:
                filepath= 'Airfoils/geometry/'+self.profile+'/'+self.profile+'.txt'
                if not os.path.exists(filepath):
                    print('No geometric data for this self.profile found. Please include the geometry in the geometry folder')
                else:
                    with open(filepath) as csvfile:
                        self.x_list = []
                        self.y_list = []    
                        self.z_list = []
                        i = 0
                    
                        for line in csvfile:
                            x,y = line.split()
            
                            self.x_list.append(eval(x))
                            self.y_list.append(eval(y))

############################
            for i in range(self.sections):
                self.scale = self.chord_list[i]
                self.distance = self.radius_list[i]                

                self.x_list_scaled = [i * self.scale for i in self.x_list]
                self.y_list_scaled = [i * self.scale for i in self.y_list]
                self.z_list = [self.distance for j in range(len(self.x_list))]                
    
                self.angle2 = self.twist_list[i]#self.twist_list[i]

                self.x_center = 1/4*max(self.x_list_scaled)
                self.y_center = 0
                self.center = np.array([self.x_center,self.y_center])
    
                self.mat = np.array([self.x_list_scaled,self.y_list_scaled])        
        
                self.rot_angle = -self.angle2*np.pi/180
                self.rot = np.array([[np.cos(self.rot_angle), -np.sin(self.rot_angle)], [np.sin(self.rot_angle), np.cos(self.rot_angle)]])
                
                for i in range(len(self.mat)):
                    for j in range(len(self.mat.T)):
                        self.mat[i,j] = self.mat[i,j]-self.center[i]
                
                self.rot_mat = self.rot @ self.mat
                for i in range(len(self.rot_mat)):
                    for j in range(len(self.rot_mat.T)):
                        self.rot_mat[i,j] = self.rot_mat[i,j]+self.center[i]
                self.x_rotated= self.rot_mat[0,:].flatten().tolist()
                self.y_rotated= self.rot_mat[1,:].flatten().tolist() 
                
                ax.plot(self.x_rotated, self.y_rotated, self.z_list, label = self.profile)
                ax.legend()
                ax.set_xlabel('x axis')
                ax.set_ylabel('y axis')
                ax.set_zlabel('z axis')
                fig.canvas.set_window_title(self.turb_name)
                plt.xlim(-(self.Radius/2), self.Radius)
                plt.ylim(-(self.Radius/2), self.Radius)
                ax.set_zlim(0,self.Radius)  

        if self.status.get() == 0:
            plt.close('3D')  
###############################################################################
root = Tk()
GUI = Bladedesign(root)
root.mainloop()