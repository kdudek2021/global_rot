import numpy as np
import time
import math

import matplotlib.pylab as plt

from matplotlib import collections  as mc



class Angular_momentum():
    def __init__(self, points, velocities, top_layer_indices, bottom_layer_indices, mass_list):
        self.points = points
        self.velocities = velocities
        self.top_layer_indices = top_layer_indices
        self.bottom_layer_indices = bottom_layer_indices
        self.mass_list = mass_list
        
        self.calculation()
        
        
        
    def calculation(self):
        
        total_L_vec = np.array([0.0, 0.0, 0.0])
        L_up_vec = np.array([0.0, 0.0, 0.0])
        L_down_vec = np.array([0.0, 0.0, 0.0])
        
        for i in range(len(self.points)):
            mass_point = self.mass_list[i]
            r_point_vec = np.array(self.points[i])
            vel_vec = self.velocities[i]
            
            momentum_vec = mass_point * vel_vec
            
            L_i = np.cross(r_point_vec, momentum_vec)
            
            total_L_vec += L_i
            
            if self.top_layer_indices.count(i) != 0:
                L_up_vec += L_i
            elif self.bottom_layer_indices.count(i) != 0:
                L_down_vec += L_i
            else:
                print("||| Angular momentum error |||")
        
        self.ang_mom_total = total_L_vec
        self.ang_mom_up = L_up_vec
        self.ang_mom_down = L_down_vec
        
        
        
    def results(self):
        return self.ang_mom_total, self.ang_mom_up, self.ang_mom_down
            
            
            
            
            