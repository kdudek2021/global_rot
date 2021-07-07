import numpy as np
import time
import math

import matplotlib.pylab as plt

from matplotlib import collections  as mc


class Geometry():
    def __init__(self, a, h0, k_L, k_S, Kh_L, Kh_S, Kh_planes_bending, m_L, m_S, Kh_twist, k_vertical, m1, m2, m3, m4):
        self.a = a
        self.h0 = h0
        self.k_L = k_L
        self.k_S = k_S
        self.Kh_L = Kh_L
        self.Kh_S = Kh_S
        self.m_L = m_L
        self.m_S = m_S
        self.Kh_twist = Kh_twist
        self.k_vertical = k_vertical
        self.Kh_planes_bending = Kh_planes_bending
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.m4 = m4
        
        self.geometry_def()
        
        
        
    def geometry_def(self):
        self.points = []
        self.top_layer_indices = []
        self.bottom_layer_indices = []
        h_triangle = 0.5 * self.a * (2.0**0.5)
        
        h = 0.0
        P0 = np.array([0.0, 0.0, h])
        P1 = np.array([-self.a / 2.0, -h_triangle, h])
        P2 = np.array([-self.a, 0.0, h])
        P3 = np.array([-self.a / 2.0, h_triangle, h])
        P4 = np.array([self.a / 2.0, h_triangle, h])
        P5 = np.array([self.a, 0.0, h])
        P6 = np.array([self.a / 2.0, -h_triangle, h])
        
        self.points.extend([P0, P1, P2, P3, P4, P5, P6])
        self.top_layer_indices.extend([0, 1, 2, 3, 4, 5, 6])
        
        h = self.h0
        P7 = np.array([0.0, 0.0, h])
        P8 = np.array([-self.a / 2.0, -h_triangle, h])
        P9 = np.array([-self.a, 0.0, h])
        P10 = np.array([-self.a / 2.0, h_triangle, h])
        P11 = np.array([self.a / 2.0, h_triangle, h])
        P12 = np.array([self.a, 0.0, h])
        P13 = np.array([self.a / 2.0, -h_triangle, h])
        
        self.points.extend([P7, P8, P9, P10, P11, P12, P13])
        self.bottom_layer_indices.extend([7, 8, 9, 10, 11, 12, 13])
        
        #P14 = np.array([(P0[0] + P7[0])/2.0 - 0.4 * self.a, (P0[1] + P7[1])/2.0, (P0[2] + P7[2])/2.0])
        #P15 = np.array([(P0[0] + P7[0])/2.0 + 0.4 * self.a, (P0[1] + P7[1])/2.0, (P0[2] + P7[2])/2.0])
        #P15 = np.array([self.a * (1.5), 0.0, h])
        #self.points.extend([P14, P15])
        #self.bottom_layer_indices.extend([14,15])
        #self.points.extend([P14, P15])
        #self.top_layer_indices.extend([14,15])
        
        #h = 0.0
        
        #P16 = np.array([-self.a * (1.5), 0.0, h])
        #P17 = np.array([self.a * (1.5), 0.0, h])
        
        #self.points.extend([P16, P17])
        #self.top_layer_indices.extend([16,17])
        
        #/////////////////////////////////////////////////////////////
        #x_transl = 2.2 * self.a
        
        #h = 0.0
        #P14 = np.array([0.0 + x_transl, 0.0, h])
        #P15 = np.array([-self.a / 2.0 + x_transl, -h_triangle, h])
        #P16 = np.array([-self.a + x_transl, 0.0, h])
        #P17 = np.array([-self.a / 2.0 + x_transl, h_triangle, h])
        #P18 = np.array([self.a / 2.0 + x_transl, h_triangle, h])
        #P19 = np.array([self.a + x_transl, 0.0, h])
        #P20 = np.array([self.a / 2.0 + x_transl, -h_triangle, h])
        
        #self.points.extend([P14, P15, P16, P17, P18, P19, P20])
        #self.bottom_layer_indices.extend([14, 15, 16, 17, 18, 19, 20])
        ##self.top_layer_indices.extend([14, 15, 16, 17, 18, 19, 20])
        
        #h = self.h0
        #P21 = np.array([0.0 + x_transl, 0.0, h])
        #P22 = np.array([-self.a / 2.0 + x_transl, -h_triangle, h])
        #P23 = np.array([-self.a + x_transl, 0.0, h])
        #P24 = np.array([-self.a / 2.0 + x_transl, h_triangle, h])
        #P25 = np.array([self.a / 2.0 + x_transl, h_triangle, h])
        #P26 = np.array([self.a + x_transl, 0.0, h])
        #P27 = np.array([self.a / 2.0 + x_transl, -h_triangle, h])
        
        #self.points.extend([P21, P22, P23, P24, P25, P26, P27])
        #self.top_layer_indices.extend([21, 22, 23, 24, 25, 26, 27])
        ##self.bottom_layer_indices.extend([21, 22, 23, 24, 25, 26, 27])
        
        #P28 = np.array([(P7[0] + P21[0])/2.0, (P7[1] + P21[1])/2.0, P7[2]])
        #P29 = np.array([(P0[0] + P14[0])/2.0, (P0[1] + P14[1])/2.0, P0[2]])

        #self.points.extend([P28, P29])
        #self.top_layer_indices.extend([28, 29])
        
        
        self.double_and_tripple_bonds()
        
        self.mass_list = Mass_assignment(self.m_L, self.m_S, self.bottom_layer_indices, self.top_layer_indices, self.points, self.m1, self.m2, self.m3, self.m4).results()
        
        
        
    def double_and_tripple_bonds(self):
        
        double_bonds_results = Double_bonds(self.k_L, self.k_S, self.k_vertical).results()
        self.double_bonds_indices = double_bonds_results[0]
        self.k_constant_list = double_bonds_results[1]
        
        tripple_bonds_results = Tripple_bonds(self.Kh_L, self.Kh_S, self.Kh_planes_bending).results()
        self.tripple_bonds_indices = tripple_bonds_results[0]
        self.Kh_list = tripple_bonds_results[1]
        
        four_body_bonds_results = Four_body_bonds(self.Kh_twist, self.points).results()
        self.four_body_bonds_indices = four_body_bonds_results[0]
        self.Kh_twist_list = four_body_bonds_results[1]
        self.eq_twist_angles = four_body_bonds_results[2]
        
        
    
    def results(self):
        return self.points, self.top_layer_indices, self.bottom_layer_indices, self.double_bonds_indices, self.k_constant_list, self.tripple_bonds_indices, self.Kh_list, self.mass_list, self.four_body_bonds_indices, self.Kh_twist_list, self.eq_twist_angles





class Double_bonds():
    def __init__(self, k_L, k_S, k_vertical):
        self.k_L = k_L
        self.k_S = k_S
        self.k_vertical = k_vertical
        
        self.assignment_of_bonds()
    
    
    def assignment_of_bonds(self):
        
        self.double_bonds_indices = []
        self.k_constant_list = []
        
        #top layer
        self.double_bonds_indices.extend([[0, 1], [0, 2], [0,3], [0,4], [0,5], [0,6]])
        self.k_constant_list.extend([self.k_L, self.k_L, self.k_L, self.k_L, self.k_L, self.k_L])

        self.double_bonds_indices.extend([[1, 2], [2, 3], [3,4], [4,5], [5,6], [6,1]])
        self.k_constant_list.extend([self.k_L, self.k_L, self.k_L, self.k_L, self.k_L, self.k_L])
        
        self.double_bonds_indices.extend([[8, 11], [9, 12], [10, 13]])
        self.k_constant_list.extend([self.k_L, self.k_L, self.k_L])
        
        #bottom layer
        self.double_bonds_indices.extend([[7, 8], [7, 9], [7,10], [7,11], [7,12], [7,13]])
        self.k_constant_list.extend([self.k_L, self.k_L, self.k_L, self.k_L, self.k_L, self.k_L])

        self.double_bonds_indices.extend([[8, 9], [9, 10], [10,11], [11,12], [12,13], [13,8]])
        self.k_constant_list.extend([self.k_L, self.k_L, self.k_L, self.k_L, self.k_L, self.k_L])  
        
        self.double_bonds_indices.extend([[3, 6], [4, 1], [2, 5]])
        self.k_constant_list.extend([self.k_L, self.k_L, self.k_L])
        
        #ligaments connecting two layers
        self.double_bonds_indices.extend([[1, 9], [2, 10], [3, 11], [4, 12], [5, 13], [6, 8]])
        self.k_constant_list.extend([ self.k_S, self.k_S, self.k_S, self.k_S, self.k_S, self.k_S])
        
        #vertical ligaments inducing the deformation
        self.double_bonds_indices.extend([[0,7]])
        self.k_constant_list.extend([self.k_vertical])
        
        #ligaments to show global rotation
        #self.double_bonds_indices.extend([[0,14], [0,15], [7,14], [7,15]])
        #self.k_constant_list.extend([ self.k_S, self.k_S, self.k_S, self.k_S])
        
        





        
    def results(self):
        return self.double_bonds_indices, self.k_constant_list





class Tripple_bonds():
    def __init__(self, Kh_L, Kh_S, Kh_planes_bending):
        
        self.assignment_of_tripple_bonds(Kh_L, Kh_S, Kh_planes_bending)
        
        
    
    def assignment_of_tripple_bonds(self, Kh_L, Kh_S, Kh_planes_bending):
        
        self.tripple_bonds_indices = []
        self.Kh_list = []
        
        #bottom layer
        self.tripple_bonds_indices.extend([[1,0,2], [2,0,3], [3,0,4], [4,0,5], [5,0,6], [6,0,1]])
        self.Kh_list.extend([Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L])
        
        self.tripple_bonds_indices.extend([[0,2,1], [2,1,0], [3,2,0], [0,3,2], [4,3,0], [0,4,3], [0,5,4], [5,4,0], [0,6,5], [6,5,0], [0,1,6], [1,6,0]])
        self.Kh_list.extend([Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L])
        
        self.tripple_bonds_indices.extend([[1,0,7], [2,0,7], [3,0,7], [4,0,7], [5,0,7], [6,0,7]])
        self.Kh_list.extend([Kh_planes_bending, Kh_planes_bending, Kh_planes_bending, Kh_planes_bending, Kh_planes_bending, Kh_planes_bending])
        
        #top layer
        self.tripple_bonds_indices.extend([[8,7,9], [9,7,10], [10,7,11], [11,7,12], [12,7,13], [13,7,8]])
        self.Kh_list.extend([Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L])
        
        self.tripple_bonds_indices.extend([[9,8,7], [7,9,8], [10,9,7], [7,10,9], [11,10,7], [7,11,10], [7,12,11], [12,11,7], [7,13,12], [13,12,7], [7,8,13], [8,13,7]])
        self.Kh_list.extend([Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L, Kh_L])
        
        self.tripple_bonds_indices.extend([[9,7,0], [10,7,0], [11,7,0], [12,7,0], [13,7,0], [8,7,0]])
        self.Kh_list.extend([Kh_planes_bending, Kh_planes_bending, Kh_planes_bending, Kh_planes_bending, Kh_planes_bending, Kh_planes_bending])


        
        ##ligaments connecting two layers
        self.tripple_bonds_indices.extend([[1,6,8], [13,8,6], [8,9,1], [2,1,9], [10,2,3], [2,10,9], [3,11,10], [11,3,4], [4,12,11], [12,4,5], [12,13,5], [6,5,13]])
        self.Kh_list.extend([Kh_S, Kh_S, Kh_S, Kh_S, Kh_S, Kh_S, Kh_S, Kh_S, Kh_S, Kh_S, Kh_S, Kh_S])


        
    
    def results(self):
        return self.tripple_bonds_indices, self.Kh_list





class Mass_assignment():
    def __init__(self, m_L, m_S, bottom_layer_indices, top_layer_indices, points, m1, m2, m3, m4):
        self.m_L = m_L
        self.m_S = m_S
        self.bottom_layer_indices = bottom_layer_indices
        self.top_layer_indices = top_layer_indices
        self.points = points
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.m4 = m4
        
        self.mass_selection()
        
        
    def mass_selection(self):
        
        #self.mass_list = []
        
        self.mass_list = [self.m3, self.m4, self.m4, self.m4, self.m4, self.m4, self.m4, self.m1, self.m2, self.m2, self.m2, self.m2, self.m2, self.m2]
        
        self.mass_list.extend([self.m1, self.m1])
        #for i in range(len(self.points)):
            #if self.top_layer_indices.count(i) == 0:
                #self.mass_list.append(self.m_S)
            #else:
                #self.mass_list.append(self.m_L)
        
    
    def results(self):
        return self.mass_list







class Four_body_bonds():
    def __init__(self, Kh_twist, points):
        self.Kh_twist = Kh_twist
        self.points = points
        
        self.indices_assignment()
                
        
    def indices_assignment(self):
        
        self.four_body_bonds_indices = [[1,6,8,13], [2,1,9,8], [3,2,10,9], [4,3,11,10], [5,4,12,11], [6,5,13,12]]
        self.Kh_twist_list = [self.Kh_twist, self.Kh_twist, self.Kh_twist, self.Kh_twist, self.Kh_twist, self.Kh_twist]
        
        #self.four_body_bonds_indices.extend([[17,14,0,6], [4,0,14,15], [3,0,14,20], [18,14,0,1]])
        #self.Kh_twist_list.extend([(10**5)*self.Kh_twist, (10**5)*self.Kh_twist, (10**5)*self.Kh_twist, (10**5)*self.Kh_twist])
        
        #self.four_body_bonds_indices.extend([[11,7,21,22], [13,7,21,24], [8,7,21,25], [10,7,21,27]])
        #self.Kh_twist_list.extend([(10**5)*self.Kh_twist, (10**5)*self.Kh_twist, (10**5)*self.Kh_twist, (10**5)*self.Kh_twist])
        
        self.eq_twist_angles = self.eq_twist_angles_calc()
        
        
    def eq_twist_angles_calc(self):
        
        self.eq_twist_angles = []
        
        for i in range(len(self.four_body_bonds_indices)):
            
            id1 = self.four_body_bonds_indices[i][0]
            id2 = self.four_body_bonds_indices[i][1]
            id3 = self.four_body_bonds_indices[i][2]
            id4 = self.four_body_bonds_indices[i][3]

            vec_p1 = np.array([self.points[id1][0], self.points[id1][1], self.points[id1][2]])
            vec_p2 = np.array([self.points[id2][0], self.points[id2][1], self.points[id2][2]])
            vec_p3 = np.array([self.points[id3][0], self.points[id3][1], self.points[id3][2]])
            vec_p4 = np.array([self.points[id4][0], self.points[id4][1], self.points[id4][2]])
            
            r_12 = vec_p1 - vec_p2
            r_23 = vec_p2 - vec_p3
            r_43 = vec_p4 - vec_p3
            
            r_12_mag = (r_12[0]**2 + r_12[1]**2 + r_12[2]**2)**0.5
            r_23_mag = (r_23[0]**2 + r_23[1]**2 + r_23[2]**2)**0.5
            r_43_mag = (r_43[0]**2 + r_43[1]**2 + r_43[2]**2)**0.5
            
            m_vec = np.cross(r_12, r_23)
            m_vec_mag = (m_vec[0]**2 + m_vec[1]**2 + m_vec[2]**2)**0.5
            
            n_vec = np.cross(r_43, r_23)
            n_vec_mag = (n_vec[0]**2 + n_vec[1]**2 + n_vec[2]**2)**0.5
            
            cos_angle = np.dot(m_vec, n_vec) / (m_vec_mag * n_vec_mag)
            
            sin_angle = (np.dot(n_vec, r_12) * r_23_mag) / (m_vec_mag * n_vec_mag)
            
            angle = - np.arctan(sin_angle / cos_angle)
            
            self.eq_twist_angles.append(angle)
        
        return self.eq_twist_angles
            
            
            
            
    def results(self):
        return self.four_body_bonds_indices, self.Kh_twist_list, self.eq_twist_angles
        
