import numpy as np
import time
import math

import matplotlib.pylab as plt

from matplotlib import collections  as mc


class Energy():
    def __init__(self, mass_list, velocities, points, double_bonds_indices, k_constant_list, eq_lengths, tripple_bonds_indices, Kh_list, eq_hinging_angles, four_body_bonds_indices, Kh_twist_list, eq_twist_angles):
        self.mass_list = mass_list
        self.velocities = velocities
        self.points = points
        
        self.double_bonds_indices = double_bonds_indices
        self.k_constant_list = k_constant_list
        self.eq_lengths = eq_lengths
        
        self.tripple_bonds_indices = tripple_bonds_indices
        self.Kh_list = Kh_list
        self.eq_hinging_angles = eq_hinging_angles
        
        self.four_body_bonds_indices = four_body_bonds_indices
        self.Kh_twist_list = Kh_twist_list
        self.eq_twist_angles = eq_twist_angles
        
        self.energies()
        
        
    
    def energies(self):
        
        self.kinetic_energy = self.kinetic_energy_calc()
        
        self.two_body_pot_energy = self.two_body_pot_energy_calc()
        self.three_body_pot_energy = self.three_body_pot_energy_calc()
        self.four_body_pot_energy = self.four_body_pot_energy_calc()
        
        self.total_energy = self.kinetic_energy + self.two_body_pot_energy + self.four_body_pot_energy + self.three_body_pot_energy# + self.four_body_pot_energy
        
        
    def kinetic_energy_calc(self):
        
        Ek = 0.0
        
        for i in range(len(self.points)):
            
            v_i = (self.velocities[i][0]**2 + self.velocities[i][1]**2 + self.velocities[i][2]**2)**0.5
            Ek += 0.5* self.mass_list[i] * (v_i**2)
        
        return Ek
    
    
    
    def two_body_pot_energy_calc(self):
        
        Ep = 0.0
        
        for i in range(len(self.double_bonds_indices)):
            k = self.k_constant_list[i]
            eq_length = self.eq_lengths[i]
            
            id1 = self.double_bonds_indices[i][0]
            id2 = self.double_bonds_indices[i][1]
            
            dist = ( (self.points[id1][0] - self.points[id2][0])**2 + (self.points[id1][1] - self.points[id2][1])**2 + (self.points[id1][2] - self.points[id2][2])**2 )**0.5
            
            Ep += 1.0 * k * ((dist - eq_length)**2)
        
        return Ep
    
    
    def three_body_pot_energy_calc(self):
        
        Ep = 0.0
        
        for j in range(len(self.tripple_bonds_indices)):
            
            id1 = self.tripple_bonds_indices[j][0]
            id2 = self.tripple_bonds_indices[j][1]
            id3 = self.tripple_bonds_indices[j][2]
            
            vec1 = np.array([self.points[id1][0] - self.points[id2][0], self.points[id1][1] - self.points[id2][1], self.points[id1][2] - self.points[id2][2]])
            
            vec2 = np.array([self.points[id3][0] - self.points[id2][0], self.points[id3][1] - self.points[id2][1], self.points[id3][2] - self.points[id2][2]])
            
            vec1_mag = (vec1[0]**2 + vec1[1]**2 + vec1[2]**2)**0.5
            vec2_mag = (vec2[0]**2 + vec2[1]**2 + vec2[2]**2)**0.5
            
            angle = np.arccos(np.dot(vec1, vec2)/ (vec1_mag * vec2_mag) )
            angle_eq = self.eq_hinging_angles[j]
            
            Kh = self.Kh_list[j]
            
            Ep += 0.5 * Kh * ((angle - angle_eq)**2)
        
        return Ep
    
    
    
    def four_body_pot_energy_calc(self):
        
        Ep = 0.0
        
        for i in range(len(self.four_body_bonds_indices)):
            
            id1 = self.four_body_bonds_indices[i][0]
            id2 = self.four_body_bonds_indices[i][1]
            id3 = self.four_body_bonds_indices[i][2]
            id4 = self.four_body_bonds_indices[i][3]
            
            r_12_vec = np.array([self.points[id1][0] - self.points[id2][0], self.points[id1][1] - self.points[id2][1], self.points[id1][2] - self.points[id2][2]])
            
            r_23_vec = np.array([self.points[id2][0] - self.points[id3][0], self.points[id2][1] - self.points[id3][1], self.points[id2][2] - self.points[id3][2]])
            
            r_23_vec_mag = (r_23_vec[0]**2 + r_23_vec[1]**2 + r_23_vec[2]**2)**0.5
            
            r_43_vec = np.array([self.points[id4][0] - self.points[id3][0], self.points[id4][1] - self.points[id3][1], self.points[id4][2] - self.points[id3][2]])
            
            m_vec = np.cross(r_12_vec, r_23_vec)
            n_vec = np.cross(r_43_vec, r_23_vec)
            
            m_vec_mag = (m_vec[0]**2 + m_vec[1]**2 + m_vec[2]**2)**0.5
            n_vec_mag = (n_vec[0]**2 + n_vec[1]**2 + n_vec[2]**2)**0.5
                        
            cos_angle = np.dot(m_vec, n_vec) / (m_vec_mag * n_vec_mag)
            sin_angle = (np.dot(n_vec, r_12_vec) * r_23_vec_mag) / (m_vec_mag * n_vec_mag)
            
            angle = - np.arctan(sin_angle / cos_angle)
            angle_eq = self.eq_twist_angles[i]
            
            Kh_twist = self.Kh_twist_list[i]
            
            Ep += Kh_twist * (angle - self.eq_twist_angles[i])**2
            
            #print cos_angle, "cos_angle"
            #print sin_angle, "sin_angle"
            #print m_vec_mag, "m_vec_mag"
            #print n_vec_mag, "n_vec_mag"
            #print Kh_twist, "Kh_twist"
            #print angle, "angle"
            #print self.eq_twist_angles[i], "self.eq_twist_angles[i]"
            #print Kh_twist * (angle - self.eq_twist_angles[i])**2, "Kh_twist * (angle - self.eq_twist_angles[i])**2"
            #time.sleep(1000)
            
        return Ep
    
    
    
    def results(self):
        return self.total_energy, self.kinetic_energy, self.two_body_pot_energy, self.three_body_pot_energy, self.four_body_pot_energy
