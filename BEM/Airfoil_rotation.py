import matplotlib.pyplot as plt
import numpy as np
import csv
import os

class Rotation(object):
    def rotate(self):
        self.profile = 'NACA_2418'
        self.sections= 1
        filepath= 'Airfoils/geometry/'+self.profile+'/'+self.profile+'.dat'
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
    
            self.chord_list = self.sections*[1]#[0.2,0.25,0.3]
            self.twist_list = self.sections*[0]
###########################            
            for i in range(self.sections):
                self.scale      = self.chord_list[i]
                self.distance   = self.sections*[0]

                self.x_list_scaled = [i * self.scale for i in self.x_list] 
                self.y_list_scaled = [i * self.scale for i in self.y_list] 
                self.z_list = [self.distance for j in range(len(self.x_list))]

                self.angle2 = self.twist_list[i]      

                self.x_center = 1/4*max(self.x_list)
                self.y_center = 1/4*max(self.y_list)
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
                
                plt.plot(self.x_rotated,self.y_rotated, label ='Angle of Attack: %d °' %self.angle2)
                plt.plot(self.x_center, self.y_center, marker ='o')
                plt.xlim(-2,2)
                plt.ylim(-0.5,0.5)
                plt.minorticks_on()
                plt.grid(True,which='major')
                plt.grid(True,which='minor', linestyle='--')         
                plt.show()
        print(self.chord_list)
if __name__ == "__main__":
    d = Rotation()
    d.rotate()