import numpy as np
import time
import math
#import sympy as sp

#import matplotlib.pylab as plt

#from matplotlib import collections  as mc

from numba import jitclass
from numba import float64 as f64
from numba import int64 as i64


#@jitclass([
    #('points', f64[:,:]),
    #('four_body_bonds_indices', i64[:,:]),
    #('Kh_twist_list', f64[:]),
    #('eq_twist_angles', f64[:]),
    #('mass_list', f64[:]),
    #('empty_vector', f64[:]),
    #('empty_vector2', f64[:]),
    #('acc_list', f64[:,:]),
#])
class Four_body_bond_force():
    def __init__(self, points, four_body_bonds_indices, Kh_twist_list, mass_list, acc_list, eq_twist_angles):
        self.points = points
        self.four_body_bonds_indices = four_body_bonds_indices
        self.Kh_twist_list = Kh_twist_list
        self.mass_list = mass_list
        self.eq_twist_angles = eq_twist_angles
        self.acc_list = acc_list
        #self.empty_vector = empty_vector
        #self.empty_vector2 = empty_vector2
        
        self.acc_update()
        
        
    def acc_update(self):
        
        for i in range(len(self.four_body_bonds_indices)):
            
            id1 = self.four_body_bonds_indices[i][0]
            id2 = self.four_body_bonds_indices[i][1]
            id3 = self.four_body_bonds_indices[i][2]
            id4 = self.four_body_bonds_indices[i][3]
            
            p1_vec = np.array(self.points[id1])
            p2_vec = np.array(self.points[id2])
            p3_vec = np.array(self.points[id3])
            p4_vec = np.array(self.points[id4])
            #p1_vec = self.points[id1]
            #p2_vec = self.points[id2]
            #p3_vec = self.points[id3]
            #p4_vec = self.points[id4]

            
            r_12_vec = p1_vec - p2_vec
            r_23_vec = p2_vec - p3_vec
            r_43_vec = p4_vec - p3_vec
            
            r_12_vec_mag = ( r_12_vec[0]**2 + r_12_vec[1]**2 + r_12_vec[2]**2 )**0.5
            r_23_vec_mag = ( r_23_vec[0]**2 + r_23_vec[1]**2 + r_23_vec[2]**2 )**0.5
            r_43_vec_mag = ( r_43_vec[0]**2 + r_43_vec[1]**2 + r_43_vec[2]**2 )**0.5
            
            m_vec = np.cross(r_12_vec, r_23_vec)
            n_vec = np.cross(r_43_vec, r_23_vec)
            
            m_vec_mag = (m_vec[0]**2 + m_vec[1]**2 + m_vec[2]**2)**0.5
            n_vec_mag = (n_vec[0]**2 + n_vec[1]**2 + n_vec[2]**2)**0.5
            
            cos_angle = np.dot(m_vec, n_vec) / (m_vec_mag * n_vec_mag)
            sin_angle = (np.dot(n_vec, r_12_vec) * r_23_vec_mag) / (m_vec_mag * n_vec_mag)
            
            #print sin_angle, "sin_angle"
            #print cos_angle, "cos_angle"
            
            angle = - np.arctan(sin_angle / cos_angle)
            angle_eq = self.eq_twist_angles[i]
            
            Kh_twist = self.Kh_twist_list[i]
            
            forces = self.forces_calculation(r_12_vec_mag, r_23_vec_mag, r_43_vec_mag, m_vec_mag, n_vec_mag, cos_angle, sin_angle, angle, angle_eq, Kh_twist, n_vec, m_vec, r_12_vec, r_23_vec, r_43_vec, p1_vec, p2_vec, p3_vec, p4_vec)
            
            acc1 = forces[0] / self.mass_list[id1]
            acc2 = forces[1] / self.mass_list[id2]           
            acc3 = forces[2] / self.mass_list[id3]
            acc4 = forces[3] / self.mass_list[id4]
            
            self.acc_list[id1] += acc1
            self.acc_list[id2] += acc2
            self.acc_list[id3] += acc3
            self.acc_list[id4] += acc4


        
    def forces_calculation(self, r_12_vec_mag, r_23_vec_mag, r_43_vec_mag, m_vec_mag, n_vec_mag, cos_angle, sin_angle, angle, angle_eq, Kh_twist, n_vec, m_vec, r_12_vec, r_23_vec, r_43_vec, p1_vec, p2_vec, p3_vec, p4_vec):
        
        term1 = - 2.0 * Kh_twist * (angle - angle_eq)
        
        aux_factor = - (1.0 / (1.0 + (sin_angle / cos_angle)**2))
        
        #term3a = (r_23_vec_mag / (m_vec_mag * n_vec_mag)) * np.array(n_vec)
        
        #term3bx = (- r_43_vec[0] * r_23_vec[2] + r_43_vec[2] * r_23_vec[0]) * (- r_23_vec[0]) + (r_43_vec[0] * r_23_vec[1] - r_43_vec[1] * r_23_vec[0]) * r_23_vec[1]
        
        #term3by = (r_43_vec[1] * r_23_vec[2] - r_43_vec[2] * r_23_vec[1]) * r_23_vec[2] + (r_43_vec[0] * r_23_vec[1] - r_43_vec[1] * r_23_vec[0]) * (-r_23_vec[0])
        
        #term3bz = (r_43_vec[1] * r_23_vec[2] - r_43_vec[2] * r_23_vec[1]) * (- r_23_vec[1]) + (- r_43_vec[0] * r_23_vec[2] + r_43_vec[2] * r_23_vec[0]) * r_23_vec[0]
        
        #term3b = (1.0 / (m_vec_mag * n_vec_mag)) * np.array([term3bx, term3by, term3bz])
        
        #factor1 = (term3a * cos_angle - term3b * sin_angle) / (cos_angle**2)
        
        
        #term2 = aux_factor * factor1
        
        #force1 = term1 * term2
        
        force1 = self.force1_calc(r_12_vec_mag, r_23_vec_mag, r_43_vec_mag, m_vec_mag, n_vec_mag, cos_angle, sin_angle, angle, angle_eq, Kh_twist, n_vec, m_vec, r_12_vec, r_23_vec, r_43_vec, p1_vec, p2_vec, p3_vec, p4_vec, term1, aux_factor)
        
        force2 = self.force2_calc(r_12_vec_mag, r_23_vec_mag, r_43_vec_mag, m_vec_mag, n_vec_mag, cos_angle, sin_angle, angle, angle_eq, Kh_twist, n_vec, m_vec, r_12_vec, r_23_vec, r_43_vec, p1_vec, p2_vec, p3_vec, p4_vec, term1, aux_factor)

        force3 = self.force3_calc(r_12_vec_mag, r_23_vec_mag, r_43_vec_mag, m_vec_mag, n_vec_mag, cos_angle, sin_angle, angle, angle_eq, Kh_twist, n_vec, m_vec, r_12_vec, r_23_vec, r_43_vec, p1_vec, p2_vec, p3_vec, p4_vec, term1, aux_factor)
        
        force4 = self.force4_calc(r_12_vec_mag, r_23_vec_mag, r_43_vec_mag, m_vec_mag, n_vec_mag, cos_angle, sin_angle, angle, angle_eq, Kh_twist, n_vec, m_vec, r_12_vec, r_23_vec, r_43_vec, p1_vec, p2_vec, p3_vec, p4_vec, term1, aux_factor)
        
        return force1, force2, force3, force4
            
    
    
    def force1_calc(self, r_12_vec_mag, r_23_vec_mag, r_43_vec_mag, m_vec_mag, n_vec_mag, cos_angle, sin_angle, angle, angle_eq, Kh_twist, n_vec, m_vec, r_12_vec, r_23_vec, r_43_vec, p1_vec, p2_vec, p3_vec, p4_vec, term1, aux_factor):
        
        r_1x = p1_vec[0]
        r_1y = p1_vec[1]
        r_1z = p1_vec[2]

        r_2x = p2_vec[0]
        r_2y = p2_vec[1]
        r_2z = p2_vec[2]

        r_3x = p3_vec[0]
        r_3y = p3_vec[1]
        r_3z = p3_vec[2]
        
        r_4x = p4_vec[0]
        r_4y = p4_vec[1]
        r_4z = p4_vec[2]
        
        term3a_x = -(r_2y - r_3y)*(-r_3z + r_4z) + (r_2z - r_3z)*(-r_3y + r_4y)
        term3a_y = (r_2x - r_3x)*(-r_3z + r_4z) - (r_2z - r_3z)*(-r_3x + r_4x)
        term3a_z = -(r_2x - r_3x)*(-r_3y + r_4y) + (r_2y - r_3y)*(-r_3x + r_4x)
        
        term3b_x = (r_2y - r_3y)*(-(r_2x - r_3x)*(-r_3y + r_4y) + (r_2y - r_3y)*(-r_3x + r_4x)) + (-r_2z + r_3z)*((r_2x - r_3x)*(-r_3z + r_4z) - (r_2z - r_3z)*(-r_3x + r_4x))
        term3b_y = (-r_2x + r_3x)*(-(r_2x - r_3x)*(-r_3y + r_4y) + (r_2y - r_3y)*(-r_3x + r_4x)) + (r_2z - r_3z)*(-(r_2y - r_3y)*(-r_3z + r_4z) + (r_2z - r_3z)*(-r_3y + r_4y))
        term3b_z = (r_2x - r_3x)*((r_2x - r_3x)*(-r_3z + r_4z) - (r_2z - r_3z)*(-r_3x + r_4x)) + (-r_2y + r_3y)*(-(r_2y - r_3y)*(-r_3z + r_4z) + (r_2z - r_3z)*(-r_3y + r_4y))

        term3a_vec = (r_23_vec_mag / (m_vec_mag * n_vec_mag)) * np.array([term3a_x, term3a_y, term3a_z])

        term3b_vec = (1.0 / (m_vec_mag * n_vec_mag)) * np.array([term3b_x, term3b_y, term3b_z])
        
        factor1 = (term3a_vec * cos_angle - term3b_vec * sin_angle) / (cos_angle**2)
        
        force1 = term1 * aux_factor * factor1
        
        return force1



    def force2_calc(self, r_12_vec_mag, r_23_vec_mag, r_43_vec_mag, m_vec_mag, n_vec_mag, cos_angle, sin_angle, angle, angle_eq, Kh_twist, n_vec, m_vec, r_12_vec, r_23_vec, r_43_vec, p1_vec, p2_vec, p3_vec, p4_vec, term1, aux_factor):
        
        r_1x = p1_vec[0]
        r_1y = p1_vec[1]
        r_1z = p1_vec[2]

        r_2x = p2_vec[0]
        r_2y = p2_vec[1]
        r_2z = p2_vec[2]

        r_3x = p3_vec[0]
        r_3y = p3_vec[1]
        r_3z = p3_vec[2]
        
        r_4x = p4_vec[0]
        r_4y = p4_vec[1]
        r_4z = p4_vec[2]
        
        term3a_x = (r_1y - r_2y)*(-r_3z + r_4z) + (r_1z - r_2z)*(r_3y - r_4y) + (r_2y - r_3y)*(-r_3z + r_4z) - (r_2z - r_3z)*(-r_3y + r_4y)
        term3a_y = (r_1x - r_2x)*(r_3z - r_4z) + (r_1z - r_2z)*(-r_3x + r_4x) - (r_2x - r_3x)*(-r_3z + r_4z) + (r_2z - r_3z)*(-r_3x + r_4x)
        term3a_z = (r_1x - r_2x)*(-r_3y + r_4y) + (r_1y - r_2y)*(r_3x - r_4x) + (r_2x - r_3x)*(-r_3y + r_4y) - (r_2y - r_3y)*(-r_3x + r_4x) 
        
        term3b_x = (-r_1y + r_3y)*(-(r_2x - r_3x)*(-r_3y + r_4y) + (r_2y - r_3y)*(-r_3x + r_4x)) + (r_1z - r_3z)*((r_2x - r_3x)*(-r_3z + r_4z) - (r_2z - r_3z)*(-r_3x + r_4x)) + (r_3y - r_4y)*((r_1x - r_2x)*(r_2y - r_3y) - (r_1y - r_2y)*(r_2x - r_3x)) + (-r_3z + r_4z)*(-(r_1x - r_2x)*(r_2z - r_3z) + (r_1z - r_2z)*(r_2x - r_3x))
        
        term3b_y = (r_1x - r_3x)*(-(r_2x - r_3x)*(-r_3y + r_4y) + (r_2y - r_3y)*(-r_3x + r_4x)) + (-r_1z + r_3z)*(-(r_2y - r_3y)*(-r_3z + r_4z) + (r_2z - r_3z)*(-r_3y + r_4y)) + (-r_3x + r_4x)*((r_1x - r_2x)*(r_2y - r_3y) - (r_1y - r_2y)*(r_2x - r_3x)) + (r_3z - r_4z)*((r_1y - r_2y)*(r_2z - r_3z) - (r_1z - r_2z)*(r_2y - r_3y))
        
        term3b_z = (-r_1x + r_3x)*((r_2x - r_3x)*(-r_3z + r_4z) - (r_2z - r_3z)*(-r_3x + r_4x)) + (r_1y - r_3y)*(-(r_2y - r_3y)*(-r_3z + r_4z) + (r_2z - r_3z)*(-r_3y + r_4y)) + (r_3x - r_4x)*(-(r_1x - r_2x)*(r_2z - r_3z) + (r_1z - r_2z)*(r_2x - r_3x)) + (-r_3y + r_4y)*((r_1y - r_2y)*(r_2z - r_3z) - (r_1z - r_2z)*(r_2y - r_3y))
        
        term3a_vec = (r_23_vec_mag / (m_vec_mag * n_vec_mag)) * np.array([term3a_x, term3a_y, term3a_z])
        term3b_vec = (1.0 / (m_vec_mag * n_vec_mag)) * np.array([term3b_x, term3b_y, term3b_z])
        
        factor1 = (term3a_vec * cos_angle - term3b_vec * sin_angle) / (cos_angle**2)
        
        force2 = term1 * aux_factor * factor1
        
        return force2
        


    def force3_calc(self, r_12_vec_mag, r_23_vec_mag, r_43_vec_mag, m_vec_mag, n_vec_mag, cos_angle, sin_angle, angle, angle_eq, Kh_twist, n_vec, m_vec, r_12_vec, r_23_vec, r_43_vec, p1_vec, p2_vec, p3_vec, p4_vec, term1, aux_factor):
        
        r_1x = p1_vec[0]
        r_1y = p1_vec[1]
        r_1z = p1_vec[2]

        r_2x = p2_vec[0]
        r_2y = p2_vec[1]
        r_2z = p2_vec[2]

        r_3x = p3_vec[0]
        r_3y = p3_vec[1]
        r_3z = p3_vec[2]
        
        r_4x = p4_vec[0]
        r_4y = p4_vec[1]
        r_4z = p4_vec[2]
        
        term3a_x = (r_1y - r_2y)*(r_2z - r_4z) + (r_1z - r_2z)*(-r_2y + r_4y)
        term3a_y = (r_1x - r_2x)*(-r_2z + r_4z) + (r_1z - r_2z)*(r_2x - r_4x)
        term3a_z = (r_1x - r_2x)*(r_2y - r_4y) + (r_1y - r_2y)*(-r_2x + r_4x)
        
        term3b_x = (r_1y - r_2y)*(-(r_2x - r_3x)*(-r_3y + r_4y) + (r_2y - r_3y)*(-r_3x + r_4x)) + (-r_1z + r_2z)*((r_2x - r_3x)*(-r_3z + r_4z) - (r_2z - r_3z)*(-r_3x + r_4x)) + (-r_2y + r_4y)*((r_1x - r_2x)*(r_2y - r_3y) - (r_1y - r_2y)*(r_2x - r_3x)) + (r_2z - r_4z)*(-(r_1x - r_2x)*(r_2z - r_3z) + (r_1z - r_2z)*(r_2x - r_3x))
        
        term3b_y = (-r_1x + r_2x)*(-(r_2x - r_3x)*(-r_3y + r_4y) + (r_2y - r_3y)*(-r_3x + r_4x)) + (r_1z - r_2z)*(-(r_2y - r_3y)*(-r_3z + r_4z) + (r_2z - r_3z)*(-r_3y + r_4y)) + (r_2x - r_4x)*((r_1x - r_2x)*(r_2y - r_3y) - (r_1y - r_2y)*(r_2x - r_3x)) + (-r_2z + r_4z)*((r_1y - r_2y)*(r_2z - r_3z) - (r_1z - r_2z)*(r_2y - r_3y))
        
        term3b_z = (r_1x - r_2x)*((r_2x - r_3x)*(-r_3z + r_4z) - (r_2z - r_3z)*(-r_3x + r_4x)) + (-r_1y + r_2y)*(-(r_2y - r_3y)*(-r_3z + r_4z) + (r_2z - r_3z)*(-r_3y + r_4y)) + (-r_2x + r_4x)*(-(r_1x - r_2x)*(r_2z - r_3z) + (r_1z - r_2z)*(r_2x - r_3x)) + (r_2y - r_4y)*((r_1y - r_2y)*(r_2z - r_3z) - (r_1z - r_2z)*(r_2y - r_3y))
        
        term3a_vec = (r_23_vec_mag / (m_vec_mag * n_vec_mag)) * np.array([term3a_x, term3a_y, term3a_z])
        term3b_vec = (1.0 / (m_vec_mag * n_vec_mag)) * np.array([term3b_x, term3b_y, term3b_z])
        
        factor1 = (term3a_vec * cos_angle - term3b_vec * sin_angle) / (cos_angle**2)
        
        force3 = term1 * aux_factor * factor1
        
        return force3




    def force4_calc(self, r_12_vec_mag, r_23_vec_mag, r_43_vec_mag, m_vec_mag, n_vec_mag, cos_angle, sin_angle, angle, angle_eq, Kh_twist, n_vec, m_vec, r_12_vec, r_23_vec, r_43_vec, p1_vec, p2_vec, p3_vec, p4_vec, term1, aux_factor):
        
        r_1x = p1_vec[0]
        r_1y = p1_vec[1]
        r_1z = p1_vec[2]

        r_2x = p2_vec[0]
        r_2y = p2_vec[1]
        r_2z = p2_vec[2]

        r_3x = p3_vec[0]
        r_3y = p3_vec[1]
        r_3z = p3_vec[2]
        
        r_4x = p4_vec[0]
        r_4y = p4_vec[1]
        r_4z = p4_vec[2]
        
        term3a_x = (r_1y - r_2y)*(-r_2z + r_3z) + (r_1z - r_2z)*(r_2y - r_3y)
        term3a_y = (r_1x - r_2x)*(r_2z - r_3z) + (r_1z - r_2z)*(-r_2x + r_3x)
        term3a_z = (r_1x - r_2x)*(-r_2y + r_3y) + (r_1y - r_2y)*(r_2x - r_3x)
        
        term3b_x = (r_2y - r_3y)*((r_1x - r_2x)*(r_2y - r_3y) - (r_1y - r_2y)*(r_2x - r_3x)) + (-r_2z + r_3z)*(-(r_1x - r_2x)*(r_2z - r_3z) + (r_1z - r_2z)*(r_2x - r_3x))
        
        term3b_y = (-r_2x + r_3x)*((r_1x - r_2x)*(r_2y - r_3y) - (r_1y - r_2y)*(r_2x - r_3x)) + (r_2z - r_3z)*((r_1y - r_2y)*(r_2z - r_3z) - (r_1z - r_2z)*(r_2y - r_3y))
        
        term3b_z = (r_2x - r_3x)*(-(r_1x - r_2x)*(r_2z - r_3z) + (r_1z - r_2z)*(r_2x - r_3x)) + (-r_2y + r_3y)*((r_1y - r_2y)*(r_2z - r_3z) - (r_1z - r_2z)*(r_2y - r_3y))
        
        #print r_23_vec_mag, "r_23_vec_mag"
        #print m_vec_mag, "m_vec_mag"
        #print n_vec_mag, "n_vec_mag"
        #print (m_vec_mag * n_vec_mag), "(m_vec_mag * n_vec_mag)"
        
        term3a_vec = (r_23_vec_mag / (m_vec_mag * n_vec_mag)) * np.array([term3a_x, term3a_y, term3a_z])
        term3b_vec = (1.0 / (m_vec_mag * n_vec_mag)) * np.array([term3b_x, term3b_y, term3b_z])
        
        factor1 = (term3a_vec * cos_angle - term3b_vec * sin_angle) / (cos_angle**2)
        
        force4 = term1 * aux_factor * factor1
        
        return force4
    
    
    
    def results(self):
        return self.acc_list
            
            
            
            
            
            
            
            
            
            
            
            