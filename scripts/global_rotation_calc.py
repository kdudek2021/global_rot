import numpy as np
import time
import math

import matplotlib.pylab as plt

from matplotlib import collections  as mc



class Global_rotation():
    def __init__(self, points_all_times):
        self.points_all_times = points_all_times
        
        self.calculation()
        
        
        
    def calculation(self):
        
        indices_opposite = [[6,8], [3,11]]
        self.list_of_angles_rot = []
        
        for i in range(len(self.points_all_times)):
            points = self.points_all_times[i]

            first_pair_id = indices_opposite[0]
            second_pair_id = indices_opposite[1]
            
            pair_1_id_1 = first_pair_id[0]
            pair_1_id_2 = first_pair_id[1]        
        
            pair_2_id_1 = second_pair_id[0]
            pair_2_id_2 = second_pair_id[1]

            point1_middle = np.array([0.5 * (points[pair_1_id_1][0] + points[pair_1_id_2][0]), 0.5 * (points[pair_1_id_1][1] + points[pair_1_id_2][1]), 0.5 * (points[pair_1_id_1][2] + points[pair_1_id_2][2])])

            point2_middle = np.array([0.5 * (points[pair_2_id_1][0] + points[pair_2_id_2][0]), 0.5 * (points[pair_2_id_1][1] + points[pair_2_id_2][1]), 0.5 * (points[pair_2_id_1][2] + points[pair_2_id_2][2])])
            
            if i == 0:
                vector_initial = point2_middle - point1_middle
                vector_initial_projection = np.array([vector_initial[0], vector_initial[1], 0.0])
                vector_initial_mag_proj = (vector_initial_projection[0]**2 + vector_initial_projection[1]**2 + vector_initial_projection[2]**2)**0.5
            else:
                vector_connect_ligaments = point2_middle - point1_middle
                vector_connect_ligaments_proj = np.array([vector_connect_ligaments[0], vector_connect_ligaments[1], 0.0])
                vector_connect_ligaments_mag_proj = (vector_connect_ligaments_proj[0]**2 + vector_connect_ligaments_proj[1]**2 + vector_connect_ligaments_proj[2]**2)**0.5
                
                angle_rot = np.arccos(np.dot(vector_initial_projection, vector_connect_ligaments_proj) / (vector_initial_mag_proj * vector_connect_ligaments_mag_proj))
                
                self.list_of_angles_rot.append(angle_rot)
    
    
    
    def results(self):
        return self.list_of_angles_rot
                
                
        
