import numpy as np
import time
import math

import matplotlib.pylab as plt

from matplotlib import collections  as mc





class XY_Visualisation_mag_moments_entire_structure():
    def __init__(self, all_pos_lists, two_body_bonds_indices, i, dt, top_layer_indices, bottom_layer_indices, points_all_times, mass_list, m1, m2, m3, m4):
        
        self.points = all_pos_lists
        self.top_layer_indices = top_layer_indices
        self.bottom_layer_indices = bottom_layer_indices
        self.points_all_times = points_all_times
        self.mass_list = mass_list
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.m4 = m4
        
        self.i_time = i * dt
        
        self.drawing(all_pos_lists, two_body_bonds_indices)#, mag_vectors, index_of_mag_mom_to_visualise)
        
        
        
    def drawing(self, all_pos_lists, two_body_bonds_indices):#, mag_vectors_to_visualise, index_of_mag_mom_to_visualise):
        
        vector_graphical_translation = [0.0, 0.0, 0.0]#[-0.025, -0.025]
        
        points_coordinates = self.points_to_plot()  
        points_top_x = points_coordinates[0]
        points_top_y = points_coordinates[1]
        points_bottom_x = points_coordinates[2]
        points_bottom_y = points_coordinates[3]
        
        different_masses_points = self.different_masses()
        points_yellow_x = different_masses_points[0]
        points_yellow_y = different_masses_points[1]
        points_red_x = different_masses_points[2]
        points_red_y= different_masses_points[3]
        points_blue_x= different_masses_points[4]
        points_blue_y= different_masses_points[5]
        points_green_x= different_masses_points[6]
        points_green_y= different_masses_points[7]

        
        
        lines = []
        
        for j in range(len(two_body_bonds_indices)):
            point1_index = two_body_bonds_indices[j][0]
            point2_index = two_body_bonds_indices[j][1]
            
            line_segment = [ (np.array([self.points[point1_index][0], self.points[point1_index][1]]) ).tolist(), (np.array([self.points[point2_index][0], self.points[point2_index][1]]) ).tolist()]
            lines.append(line_segment)
        
        plt.clf()#///////////////////////////////////////////////////////////clear canvas
        
        rot_lines_res = self.additional_rot_lines()
        
        p1_initial = rot_lines_res[0]
        p2_initial = rot_lines_res[1]
        p1_current = rot_lines_res[2]
        p2_current = rot_lines_res[3]
        
        lines_initial_rot_reference = []
        lines_current_rot_reference = []
        
        line_initial_segment = [[p1_initial[0], p1_initial[1]], [p2_initial[0], p2_initial[1]] ]
        line_current_segment = [[p1_current[0], p1_current[1]], [p2_current[0], p2_current[1]] ]
        
        lines_initial_rot = [line_initial_segment]
        lines_current_rot = [line_current_segment]
        
        
        #plt.xlim([-0.08, 0.08])
        #plt.ylim([-0.08, 0.08])
        plt.xlim(-0.4, 0.4) 
        plt.ylim(-0.2, 0.8)
        
        lc = mc.LineCollection(lines, linewidths = 1, color = 'black')
        
        lc_initial_rot = mc.LineCollection(lines_initial_rot, linewidths = 1, color = 'black')
        lc_current_rot = mc.LineCollection(lines_current_rot, linewidths = 1, color = 'red')        
        #lc_vec = mc.LineCollection(lines_mag_vectors, linewidths = 2, color = 'red')
        
        #plt.axis('equal')
        #plt.axis([-0.05, 0.15, -0.05, 0.1])
        #plt.axes([-0.5, -0.5, 1.0, 1.0])
        #ax = plt.axes([0., 0., 1., 1.])
        
        ax = plt.gca()		#I move both xticks abd yticks away from the graph
        ax.set_aspect('equal')
        ax.set_aspect(aspect=1.0)
        
        ax.plot(points_yellow_x, points_yellow_y, 'yo', markersize = 8)
        ax.plot(points_red_x, points_red_y, 'ro', markersize = 8)
        ax.plot(points_blue_x, points_blue_y, 'bo', markersize = 8)
        ax.plot(points_green_x, points_green_y, 'go', markersize = 8)
        
        #ax.tick_params(direction='out', pad=10)
        plt.draw()
        
        ax.add_collection(lc)
        ax.add_collection(lc_initial_rot)
        ax.add_collection(lc_current_rot)
        #ax.add_collection(lc_vec)
        #plt.show()
        plt.savefig("structure_XY_at_t = "+str(self.i_time)+".png") 
        #plt.savefig("structure_XY_at_t = "+str(self.i_time)+".svg") 



    def additional_rot_lines(self):
        
        initial_points = self.points_all_times[0]   
        current_points = self.points_all_times[-1]
        
        points_both_times_moments = [initial_points, current_points]
        
        indices_opposite = [[6,8], [3,11]]
        
        x_initial_vec = []
        y_initial_vec = []
        z_initial_vec = []

        x_current_vec = []
        y_current_vec = []
        z_current_vec = []
        
        for i in range(2):
            points = points_both_times_moments[i]

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
                vector_initial_mag = (vector_initial[0]**2 + vector_initial[1]**2 + vector_initial[2]**2)**0.5
                #point1_middle_extended = point1_middle - vector_initial * 1.0
                #point2_middle_extended = point1_middle + vector_initial * 2.0
                
                point1_middle_extended = np.array([0.0, 0.0, point2_middle[2]])
                point2_middle_extended = np.array([0.0, 0.0, point2_middle[2]]) + (vector_initial / vector_initial_mag) * 0.35
                
                p1_initial = point1_middle_extended
                p2_initial = point2_middle_extended

            else:
                vector_current = point2_middle - point1_middle
                vector_current_mag = (vector_current[0]**2 + vector_current[1]**2 + vector_current[2]**2)**0.5
                #point1_middle_extended = point1_middle - vector_current * 1.0
                #point2_middle_extended = point1_middle + vector_current * 2.0

                point1_middle_extended = np.array([0.0, 0.0, point2_middle[2]])
                point2_middle_extended = np.array([0.0, 0.0, point2_middle[2]]) + (vector_current / vector_current_mag) * 0.35
                
                p1_current = point1_middle_extended
                p2_current = point2_middle_extended
        
        return p1_initial, p2_initial, p1_current, p2_current
    
    
    def points_to_plot(self):
        
        points_top_y = []
        points_top_z = []
        points_bottom_y = []
        points_bottom_z = []
        
        for i in self.top_layer_indices:
            #print i, "i"
            points_top_y.append(self.points[i][0])
            points_top_z.append(self.points[i][1])
        
        for j in self.bottom_layer_indices:
            #print j, "j"
            points_bottom_y.append(self.points[j][0])
            points_bottom_z.append(self.points[j][1])
        
        points_top_y = np.array(points_top_y)
        points_top_z = np.array(points_top_z)       
        points_bottom_y = np.array(points_bottom_y)
        points_bottom_z = np.array(points_bottom_z)
        
        
        return points_top_y, points_top_z, points_bottom_y, points_bottom_z
        
        
        
    def different_masses(self):
        points_yellow_x = []
        points_yellow_y = []
        
        points_red_x = []
        points_red_y = []
        
        points_blue_x = []
        points_blue_y = []
        
        points_green_x = []
        points_green_y = []
        
        for i in range(len(self.points)):
            
            if self.mass_list[i] == self.m1:
                points_yellow_x.append(self.points[i][0])
                points_yellow_y.append(self.points[i][1])
            elif self.mass_list[i] == self.m2:
                points_red_x.append(self.points[i][0])
                points_red_y.append(self.points[i][1])    
            elif self.mass_list[i] == self.m4:
                points_blue_x.append(self.points[i][0])
                points_blue_y.append(self.points[i][1])    
            elif self.mass_list[i] == self.m3:
                points_green_x.append(self.points[i][0])
                points_green_y.append(self.points[i][1])
            else:
                None
        
        return points_yellow_x, points_yellow_y, points_red_x, points_red_y, points_blue_x, points_blue_y, points_green_x, points_green_y
