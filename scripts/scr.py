import numpy as np
import time
import math

import matplotlib.pylab as plt

from matplotlib import collections  as mc

#////////////////
from geometry_and_bonds import Geometry
from XY_visualisation_structure import XY_Visualisation_mag_moments_entire_structure
from YZ_visualisation_structure import YZ_Visualisation_mag_moments_entire_structure
from XYZ_3D_visualisation import Entire_structure
from eq_lengths_tripple_bond_angles import Eq_lengths_hinging_angles
from interaction1 import Double_bond_force
from interaction2 import Tripple_bond_force
from four_body_bonded_inter_expressions import Expressions_for_four_points
from interaction3 import Four_body_bond_force
from velocities_assign import Velocities_assignment
from dynamics_two_layers import Dynamics
from three_body_bonded_inter_expressions import Expressions_for_three_points
from energy_conservation import Energy
from angular_momentum_calc import Angular_momentum
from global_rotation_calc import Global_rotation


if __name__ == '__main__':
    a = 0.2
    h0 = 0.3
    
    
    m_total = 0.16
    
    n_multiplier = 1#200 
    m_S = 0.08
    m_L = 0.08 * n_multiplier
    
    m_top_layer = 0.5# [kg] 0.08 * 6.0
    m1 = m_top_layer / 7.0
    m2 = (m_top_layer - m1) / 6.0
    m4 = m2 / 5.0#m2 / 5.0
    m3 = m_top_layer - (6.0 * m4)

    k_L = 10**8
    k_S = 10**9
    k_vertical = 120.0
    
    Kh_L = 1.0
    Kh_S = 1.0
    Kh_planes_bending = 1000.0
    
    Kh_twist = 0.25
    
    h_eq = 0.18
        
    results_geometry = Geometry(a, h0, k_L, k_S, Kh_L, Kh_S, Kh_planes_bending, m_L, m_S, Kh_twist, k_vertical, m1, m2, m3, m4).results()
    
    points = results_geometry[0]
    top_layer_indices = results_geometry[1]
    bottom_layer_indices = results_geometry[2]
    double_bonds_indices = results_geometry[3]
    k_constant_list = results_geometry[4]
    tripple_bonds_indices = results_geometry[5]
    Kh_list = results_geometry[6]
    mass_list = results_geometry[7]
    four_body_bonds_indices = results_geometry[8]
    Kh_twist_list = results_geometry[9]
    eq_twist_angles = results_geometry[10]
    
    results_eq_lengths_tripple_bond_angles = Eq_lengths_hinging_angles(points, double_bonds_indices, tripple_bonds_indices, k_constant_list, k_vertical, h_eq, a).results()
    
    eq_lengths = results_eq_lengths_tripple_bond_angles[0]
    eq_hinging_angles = results_eq_lengths_tripple_bond_angles[1]
    
    velocities = Velocities_assignment(points).results()
    
    dt =  10**(-6)#10**(-6)
    t = 1.0 * 10**(-1)#1.0 * 10**(-1)#1.2 * 10**(-1)#2.0 * 10**(-1)#5.0 * 10**(-1)
    
    with open("initial_parameters_ver6_dt="+ str(dt) +"_k_L=" + str(k_L) + "_k_S=" + str(k_S) + ".dat", "w") as w:
        w.write("a = " + str(a))
        w.write("\n")
        w.write("m_top_layer = " + str(m_top_layer))
        w.write("\n")
        w.write("h0 = " + str(h0))
        w.write("\n")
        w.write("m_S = " + str(m_S))
        w.write("\n")
        w.write("m_L = " + str(m_L))
        w.write("\n")
        w.write("k_L = " + str(k_L))
        w.write("\n")
        w.write("k_S = " + str(k_S))
        w.write("\n")
        w.write("k_vertical = " + str(k_vertical))
        w.write("\n")
        w.write("Kh_L = " + str(Kh_L))
        w.write("\n")
        w.write("Kh_S = " + str(Kh_S))
        w.write("\n")
        w.write("Kh_twist = " + str(Kh_twist))
        w.write("\n")
        w.write("h_eq = " + str(h_eq))
        w.write("\n")
        w.write("dt = " + str(dt))
        w.write("\n")
        w.write("t = " + str(t))
        w.write("\n")
    
    i = 0
    
    #Expressions_for_four_points() ////////////////////////////////////////////////////////////////////////////////// IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!
    #print " "  
    
    #Expressions_for_three_points()  #////////////////////////////////////////////////////////////////////////////////// IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!
    points_all_times = []
    
    with open("energy_conservation_outer_masses_ratio_5.dat", "w") as g:
        for j in np.arange(0.0, t, dt):
            
            dynamics_results = Dynamics(points, velocities, mass_list, double_bonds_indices, k_constant_list, eq_lengths, tripple_bonds_indices, Kh_list, eq_hinging_angles, four_body_bonds_indices, Kh_twist_list, eq_twist_angles, dt).results()
            
            points = dynamics_results[0]
            velocities = dynamics_results[1]
            
            #print points[0][0], points[0][1]
            
            points_all_times.append(points)
            
            ang_momentum_results = Angular_momentum(points, velocities, top_layer_indices, bottom_layer_indices, mass_list).results()
            ang_mom_total = ang_momentum_results[0]
            ang_mom_up = ang_momentum_results[1]
            ang_mom_down = ang_momentum_results[2]
            
            ang_mom_total_mag = (ang_mom_total[0]**2 + ang_mom_total[1]**2 + ang_mom_total[2]**2)**0.5
            ang_mom_up_mag = (ang_mom_up[0]**2 + ang_mom_up[1]**2 + ang_mom_up[2]**2)**0.5
            ang_mom_down_mag = (ang_mom_down[0]**2 + ang_mom_down[1]**2 + ang_mom_down[2]**2)**0.5
            
            energy_results = Energy(mass_list, velocities, points, double_bonds_indices, k_constant_list, eq_lengths, tripple_bonds_indices, Kh_list, eq_hinging_angles, four_body_bonds_indices, Kh_twist_list, eq_twist_angles).results()
            
            total_energy = energy_results[0]
            kinetic_energy = energy_results[1]
            two_body_pot_energy = energy_results[2]
            three_body_pot_energy = energy_results[3]
            four_body_pot_energy = energy_results[4]
            
            if i == 0 or i % 100 == 0:
                g.write(str(i * dt))                        #1
                g.write(" ")
                g.write(str(total_energy))                  #2
                g.write(" ")
                g.write(str(kinetic_energy))                #3
                g.write(" ")
                g.write(str(two_body_pot_energy))           #4
                g.write(" ")
                g.write(str(three_body_pot_energy))         #5
                g.write(" ")
                g.write(str(four_body_pot_energy))          #6
                g.write(" ")
                g.write(str(ang_mom_total_mag))             #7
                g.write(" ")
                g.write(str(ang_mom_up_mag))                #8
                g.write(" ")
                g.write(str(ang_mom_down_mag))              #9
                g.write(" ")
                g.write(str(ang_mom_up_mag - ang_mom_down_mag)) #10
                g.write('\n')

            
            if i == 0 or i % 10000 == 0:
                None
            
            i += 1
    
    
    list_of_angle_global_rot = Global_rotation(points_all_times).results()
    
    with open("data_aux.dat", "w") as m:
        for i in range(len(list_of_angle_global_rot)):
            global_rot_extent = list_of_angle_global_rot[i] * (180.0 / np.pi)#[deg]
            
            points_at_a_given_time_step = points_all_times[i]
            points_initial = points_all_times[0]
            
            point_centre_top_z_initial = points_initial[7][2]
            point_centre_top_z = points_at_a_given_time_step[7][2]
            point_centre_bottom_z_initial = points_initial[0][2]
            point_centre_bottom_z = points_at_a_given_time_step[0][2]
            
            change_in_height_top_z = abs(point_centre_top_z_initial - point_centre_top_z)
            change_in_height_bottom_z = abs(point_centre_bottom_z_initial - point_centre_bottom_z)
            
            m.write(str(i * dt))
            m.write(" ")
            #m.write(str(global_rot_extent))#2
            #m.write(" ")
            #m.write(str(point_centre_top_z))#3
            #m.write(" ")
            #m.write(str(point_centre_bottom_z))#4
            #m.write(" ")
            #m.write(str(change_in_height_top_z))#5
            #m.write(" ")
            m.write(str(change_in_height_bottom_z))#6
            m.write('\n')
        
    

