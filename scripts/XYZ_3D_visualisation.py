import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt




class Entire_structure():
    def __init__(self, all_pos_lists, two_body_bonds_indices, i, dt, top_layer_indices, bottom_layer_indices, k_constant_list, k_L, k_S, points_all_times, mass_list, m1, m2, m3, m4):
        
        self.points = all_pos_lists
        self.top_layer_indices = top_layer_indices
        self.bottom_layer_indices = bottom_layer_indices
        self.two_body_bonds_indices = two_body_bonds_indices
        
        self.mass_list = mass_list
        self.m1 = m1
        self.m2 = m2
        self.m3 = m3
        self.m4 = m4
        
        self.k_constant_list = k_constant_list
        self.k_L = k_L
        self.k_S = k_S
        
        self.points_all_times = points_all_times
        
        self.i_time = i * dt
        
        self.drawing(all_pos_lists)#, mag_vectors, index_of_mag_mom_to_visualise)
        
        
        
    def drawing(self, all_pos_lists):#, mag_vectors_to_visualise, index_of_mag_mom_to_visualise):
        
        mpl.rcParams['legend.fontsize'] = 10

        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax._axis3don = False
        
        vertices = self.points_to_plot()
        points_top_x = vertices[0]
        points_top_y = vertices[1]
        points_top_z = vertices[2]
        points_bottom_x = vertices[3]
        points_bottom_y = vertices[4]
        points_bottom_z = vertices[5]
        
        different_masses_points = self.different_masses()
        points_yellow_x = different_masses_points[0]
        points_yellow_y = different_masses_points[1]
        points_yellow_z = different_masses_points[2]
        points_red_x = different_masses_points[3]
        points_red_y = different_masses_points[4]
        points_red_z = different_masses_points[5]
        points_blue_x = different_masses_points[6]
        points_blue_y = different_masses_points[7]
        points_blue_z = different_masses_points[8]
        points_green_x = different_masses_points[9]
        points_green_y = different_masses_points[10]
        points_green_z = different_masses_points[11]
        
        
        plt.axis('equal')
        ax.set_xlim3d(-0.5, 0.5)
        ax.set_ylim3d(-0.5,0.5)
        ax.set_zlim3d(-0.05,0.5)
        
        for i in range(len(self.two_body_bonds_indices)):
            id1 = self.two_body_bonds_indices[i][0]
            id2 = self.two_body_bonds_indices[i][1]
            
            x_coor_two_points = np.array([self.points[id1][0], self.points[id2][0]]) 
            y_coor_two_points = np.array([self.points[id1][1], self.points[id2][1]])
            z_coor_two_points = np.array([self.points[id1][2], self.points[id2][2]])
            
            if self.k_constant_list[i] == self.k_L:
                ax.plot(x_coor_two_points, y_coor_two_points, z_coor_two_points, c = 'black',label='parametric curve')
            else:
                ax.plot(x_coor_two_points, y_coor_two_points, z_coor_two_points, c = 'green',label='parametric curve')
        
        
        rot_lines_res = self.additional_rot_lines()
        
        x_initial_vec = rot_lines_res[0]
        y_initial_vec = rot_lines_res[1]
        z_initial_vec = rot_lines_res[2]
        
        #print x_initial_vec, "x_initial_vec"
        #print y_initial_vec, "y_initial_vec"
        #print z_initial_vec, "z_initial_vec"
                
        x_current_vec = rot_lines_res[3]
        y_current_vec = rot_lines_res[4]
        z_current_vec = rot_lines_res[5]
        
        #print x_current_vec, "x_current_vec"
        #print y_current_vec, "y_current_vec"
        #print z_current_vec, "z_current_vec"
        
        ax.plot(x_initial_vec, y_initial_vec, z_initial_vec, c = 'black',label='parametric curve')
        ax.plot(x_current_vec, y_current_vec, z_current_vec, c = 'red',label='parametric curve')
        
        #for j in range(len(points_top_x)):
            ##ax.scatter3D(points_top_x, points_top_y, points_top_z, c='blue', s = 50)
            #ax.scatter3D(points_top_x, points_top_y, points_top_z, c='blue')
        #for j in range(len(points_bottom_x)):
            ##ax.scatter3D(points_bottom_x, points_bottom_y, points_bottom_z, c='red', s = 50)
            #ax.scatter3D(points_bottom_x, points_bottom_y, points_bottom_z, c='red')
        
        for j in range(len(points_yellow_x)):
            #ax.scatter3D(points_top_x, points_top_y, points_top_z, c='blue', s = 50)
            ax.scatter3D(points_yellow_x, points_yellow_y, points_yellow_z, c='yellow', s=50)
        for j in range(len(points_red_x)):
            #ax.scatter3D(points_bottom_x, points_bottom_y, points_bottom_z, c='red', s = 50)
            ax.scatter3D(points_red_x, points_red_y, points_red_z, c='red', s=50)
        for j in range(len(points_blue_x)):
            #ax.scatter3D(points_bottom_x, points_bottom_y, points_bottom_z, c='red', s = 50)
            ax.scatter3D(points_blue_x, points_blue_y, points_blue_z, c='blue', s=50)
        for j in range(len(points_green_x)):
            #ax.scatter3D(points_bottom_x, points_bottom_y, points_bottom_z, c='red', s = 50)
            ax.scatter3D(points_green_x, points_green_y, points_green_z, c='green', s=50)

        
        #plt.axis('equal')
        #ax.legend()
        #plt.show()
        plt.savefig("structure_XYZ_at_t="+str("%.4f" % self.i_time)+".png") 
        #plt.savefig("structure_XYZ_at_t="+str("%.4f" % self.i_time)+".svg")         
        
        
        
    
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
                
                #x_initial_vec.extend([point1_middle_extended[0], point2_middle_extended[0]])
                #y_initial_vec.extend([point1_middle_extended[1], point2_middle_extended[1]])
                #z_initial_vec.extend([point1_middle_extended[2], point2_middle_extended[2]])

                point2_middle_extended =  np.array([0.0, 0.0, point2_middle[2]]) + (vector_initial / vector_initial_mag) * 0.35
                
                x_initial_vec.extend([0.0, -point2_middle_extended[0]*1.5])
                y_initial_vec.extend([0.0, -point2_middle_extended[1]*1.5])
                z_initial_vec.extend([point2_middle[2], point2_middle_extended[2]])
            
            else:
                vector_current = point2_middle - point1_middle
                vector_current_mag = (vector_current[0]**2 + vector_current[1]**2 + vector_current[2]**2)**0.5
                #point1_middle_extended = point1_middle - vector_current * 1.0
                #point2_middle_extended = point1_middle + vector_current * 2.0
                
                #x_current_vec.extend([point1_middle_extended[0], point2_middle_extended[0]])
                #y_current_vec.extend([point1_middle_extended[1], point2_middle_extended[1]])
                #z_current_vec.extend([point1_middle_extended[2], point2_middle_extended[2]])
                
                point2_middle_extended = np.array([0.0, 0.0, point2_middle[2]]) + (vector_current / vector_current_mag) * 0.35
                
                x_current_vec.extend([0.0, -point2_middle_extended[0]*1.5])
                y_current_vec.extend([0.0, -point2_middle_extended[1]*1.5])
                z_current_vec.extend([point2_middle[2], point2_middle_extended[2]])
                
        return x_initial_vec, y_initial_vec, z_initial_vec, x_current_vec, y_current_vec, z_current_vec



        
    def points_to_plot(self):
        
        points_top_x = []
        points_top_y = []
        points_top_z = []

        points_bottom_x = []
        points_bottom_y = []
        points_bottom_z = []
        
        for i in self.top_layer_indices:
            #print i, "i"
            points_top_x.append(self.points[i][0])
            points_top_y.append(self.points[i][1])
            points_top_z.append(self.points[i][2])
        
        for j in self.bottom_layer_indices:
            #print j, "j"
            points_bottom_x.append(self.points[j][0])
            points_bottom_y.append(self.points[j][1])
            points_bottom_z.append(self.points[j][2])
        
        points_top_x = np.array(points_top_x)
        points_top_y = np.array(points_top_y)       
        points_top_z = np.array(points_top_z)
        points_bottom_x = np.array(points_bottom_x)
        points_bottom_y = np.array(points_bottom_y)
        points_bottom_z = np.array(points_bottom_z)        
        
        return points_top_x, points_top_y, points_top_z, points_bottom_x, points_bottom_y, points_bottom_z
    
    
    
    def different_masses(self):
        
        points_yellow_x = []
        points_yellow_y = []
        points_yellow_z = []
        
        points_red_x = []
        points_red_y = []
        points_red_z = []
        
        points_blue_x = []
        points_blue_y = []
        points_blue_z = []
        
        points_green_x = []
        points_green_y = []
        points_green_z = []

        for i in range(len(self.points)):
            
            if self.mass_list[i] == self.m1:
                points_yellow_x.append(self.points[i][0])
                points_yellow_y.append(self.points[i][1])
                points_yellow_z.append(self.points[i][2])
            elif self.mass_list[i] == self.m2:
                points_red_x.append(self.points[i][0])
                points_red_y.append(self.points[i][1])    
                points_red_z.append(self.points[i][2])
            elif self.mass_list[i] == self.m4:
                points_blue_x.append(self.points[i][0])
                points_blue_y.append(self.points[i][1])
                points_blue_z.append(self.points[i][2])
            elif self.mass_list[i] == self.m3:
                points_green_x.append(self.points[i][0])
                points_green_y.append(self.points[i][1])
                points_green_z.append(self.points[i][2])
            else:
                None
        
        return points_yellow_x, points_yellow_y, points_yellow_z, points_red_x, points_red_y, points_red_z, points_blue_x, points_blue_y, points_blue_z, points_green_x, points_green_y, points_green_z
