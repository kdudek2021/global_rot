import numpy as np
import time
import math

import matplotlib.pylab as plt

from matplotlib import collections  as mc


class Eq_lengths_hinging_angles():
    def __init__(self, points, double_bonds_indices, tripple_bonds_indices, k_constant_list, k_vertical, h_eq, a):
        self.points = points
        self.double_bonds_indices = double_bonds_indices
        self.tripple_bonds_indices = tripple_bonds_indices
        self.k_constant_list = k_constant_list
        self.k_vertical = k_vertical
        self.h_eq = h_eq
        self.a = a
        
        self.equilibrium_parameters()
        
        
    
    def equilibrium_parameters(self):
        
        self.eq_lengths = self.eq_lengths_calculation()
        self.eq_hinging_angles = self.eq_tripple_bond_angles()        
        
        
        
    def eq_lengths_calculation(self):
        
        eq_lenghts = []
        
        for i in range(len(self.double_bonds_indices)):
            
            if self.k_constant_list[i] != self.k_vertical:
                id1 = self.double_bonds_indices[i][0]
                id2 = self.double_bonds_indices[i][1]
                
                p1 = self.points[id1]
                p2 = self.points[id2]
                
                dist = ( (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2 )**0.5
                
                eq_lenghts.append(dist)
                
            else:
                eq_lenghts.append(self.h_eq)
                
        
        return eq_lenghts
    
    
    
    def eq_tripple_bond_angles(self):
        #We assume that all angles are smaller than 180 degrees !!!!!!!!!!!!!!!!!!
        eq_angles = []
        
        for i in range(len(self.tripple_bonds_indices)):
            
            id1 = self.tripple_bonds_indices[i][0]
            id2 = self.tripple_bonds_indices[i][1]        
            id3 = self.tripple_bonds_indices[i][2]
            
            p1 = np.array(self.points[id1])
            p2 = np.array(self.points[id2])
            p3 = np.array(self.points[id3])
            
            v1 = p1 - p2
            v2 = p3 - p2
            
            v1_mag = (v1[0]**2 + v1[1]**2 + v1[2]**2)**0.5
            v2_mag = (v2[0]**2 + v2[1]**2 + v2[2]**2)**0.5   
            
            angle = np.arccos(np.dot(v1, v2) / (v1_mag * v2_mag))
            
            eq_angles.append(angle)
        
        return eq_angles
            
            
            
    
    def results(self):
        return self.eq_lengths, self.eq_hinging_angles