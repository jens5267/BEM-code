# -*- coding: utf-8 -*-
"""
Created on Fri Nov  2 13:15:41 2018

@author: Jens Reichwein
"""

import matplotlib.pyplot as plt
import numpy as np
import itertools

def airfoil_polars():
    family  = 'NACA'
    family2 = 'S'
    family3 = 'SG'
    
    series = [2408,2412,2415,2418]
    series1 = [4412,4415,4418]
    series2 = [822, 1210, 2091, 4061, 4180, 4320]
    series3 = [6040,6041,6042,6043]
    
    airfoils = []
    count = 0
    for number in series3:
        airfoils.append(family3 + str(number))
        
   # for number2 in series2:
   #     airfoils.append(family2 + str(number2))
        
   # airfoils.append('SD7034')
     
   # for number3 in series3:
   #     airfoils.append(family3 + str(number3))
        

    for airfoil in airfoils:

        with open(airfoil +'/' + airfoil + '.txt', 'r') as fp:

            data = fp.read()
            x   = []
            y   = []

            for line in data.splitlines():
                temp    = line.split()
             
                X   = eval(temp[0])
                Y   = eval(temp[1])
                x.append(X)
                y.append(Y)
                
        plt.plot(x,y, label = airfoil, linewidth = 1)
        plt.legend(loc = 'upper right')
        plt.xlim(0,1)
        plt.ylim(-0.5, 0.5)
        plt.xticks([])
        plt.yticks([])   
        plt.savefig(fname= 'First.pdf', dpi = 1000, format = 'pdf', facecolor  = 'w')
airfoil_polars()