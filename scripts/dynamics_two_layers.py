import numpy as np
import time
import math
import copy

import matplotlib.pylab as plt

from matplotlib import collections  as mc

from two_body_bonded_acceleration import Double_bond_force
from three_body_bonded_acceleration import Tripple_bond_force
from four_body_bonded_acceleration import Four_body_bond_force


class Dynamics():
    def __init__(self, points, velocities, mass_list, double_bonds_indices, k_constant_list, eq_lengths, tripple_bonds_indices, Kh_list, eq_hinging_angles, four_body_bonds_indices, Kh_twist_list, eq_twist_angles, dt):
        self.points = points
        self.velocities = velocities
        self.mass_list = mass_list
        
        self.double_bonds_indices = double_bonds_indices
        self.k_constant_list = k_constant_list
        self.eq_lengths = eq_lengths
        
        self.tripple_bonds_indices = tripple_bonds_indices
        self.Kh_list = Kh_list
        self.eq_hinging_angles = eq_hinging_angles
        
        self.four_body_bonds_indices = four_body_bonds_indices
        self.Kh_twist_list = Kh_twist_list
        self.eq_twist_angles = eq_twist_angles
        
        self.dt = dt
        
        self.algorithm_of_motion()
        #self.yoshida_integration()
        
    def algorithm_of_motion(self):
        
        pos = np.array(self.points) * 1.0
        vel = np.array(self.velocities) * 1.0

        k1 = self.dt * self.acceleration(pos, vel)
        l1 = self.dt * vel
        
        k2 = self.dt * self.acceleration(pos + l1 / 2.0, vel + k1 / 2.0)
        l2 = self.dt * (vel + k1 / 2.0)
        
        k3 = self.dt * self.acceleration(pos + l2 / 2.0, vel + k2 / 2.0)
        l3 = self.dt * (vel + k2 / 2.0)
        
        k4 = self.dt * self.acceleration(pos + l3, vel + k3)
        l4 = self.dt * (vel + k3)
        
        pos = pos + (1.0 / 6.0) * (l1 + 2.0 * l2 + 2.0 * l3 + l4)
        vel = vel + (1.0 / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        
        self.pos = pos
        self.vel = vel
    
    
    
    def yoshida_integration(self):
        pos = np.array(self.points) * 1.0
        vel = np.array(self.velocities) * 1.0
        
        w0 = - ((2.0)**(1.0/3.0)) / (2.0 - (2.0)**(1.0/3.0))
        w1 = 1.0 / (2.0 - (2.0)**(1.0/3.0))
        c1 = w1 / 2.0
        c4 = w1 / 2.0
        c2 = (w0 + w1) / 2.0
        c3 = (w0 + w1) / 2.0
        d1 = w1
        d3 = w1
        d2 = w0
        
        pos_1 = pos + c1 * vel * self.dt
        vel_1 = vel + d1 * self.acceleration(pos_1, vel) * self.dt
        
        pos_2 = pos_1 + c2 * vel_1 * self.dt
        vel_2 = vel_1 + d2 * self.acceleration(pos_2, vel) * self.dt
        
        pos_3 = pos_2 + c3 * vel_2 * self.dt
        vel_3 = vel_2 + d3 * self.acceleration(pos_3, vel) * self.dt
        
        pos_new = pos_3 + c4 * vel_3 * self.dt
        vel_new = vel_3
        
        self.pos = pos_new
        self.vel = vel_new
        
    
    
    def acceleration(self, pos, vel):
        
        acc_list = []
        for i in range(len(pos)):
            acc_list.append(np.array([0.0, 0.0, 0.0]))
        
        pos = pos.tolist()
        
        #acc_list = Double_bond_force(pos, self.double_bonds_indices, self.mass_list, acc_list, self.k_constant_list, self.eq_lengths).results()
        
        #acc_list = Tripple_bond_force(pos, self.tripple_bonds_indices, self.Kh_list, self.mass_list, acc_list, self.eq_hinging_angles).results()
        
        #acc_list = Four_body_bond_force(pos, self.four_body_bonds_indices, self.Kh_twist_list, self.mass_list, acc_list, self.eq_twist_angles).results()
        
        acc_list = Double_bond_force(np.array(pos, np.float64), np.array(self.double_bonds_indices, np.int64), np.array(self.mass_list, np.float64), np.array(acc_list, np.float64), np.array(self.k_constant_list, np.float64), np.array(self.eq_lengths, np.float64)).results()
        
        acc_list = Tripple_bond_force(np.array(pos, np.float64), np.array(self.tripple_bonds_indices, np.int64), np.array(self.Kh_list, np.float64), np.array(self.mass_list, np.float64), np.array(acc_list, np.float64), np.array(self.eq_hinging_angles, np.float64)).results()
        
        acc_list = Four_body_bond_force(pos, self.four_body_bonds_indices, self.Kh_twist_list, self.mass_list, acc_list, self.eq_twist_angles).results()
        
        #empty_vector = np.array([0.0, 0.0, 0.0], np.float64)
        #empty_vector2 = np.array([0.0, 0.0, 0.0], np.float64)
        #acc_list = Four_body_bond_force(np.array(pos, np.float64), np.array(self.four_body_bonds_indices, np.int64), np.array(self.Kh_twist_list, np.float64), np.array(self.mass_list, np.float64), np.array(acc_list, np.float64), np.array(self.eq_twist_angles, np.float64), empty_vector, empty_vector2).results()
        
        return np.array(acc_list)
    
    
    def results(self):
        return self.pos, self.vel
        