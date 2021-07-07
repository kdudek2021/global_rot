#import time
import math

import numpy as np
from numba import jitclass
from numba import float64 as f64
from numba import int64 as i64


@jitclass([
    ('points', f64[:,:]),
    ('tripple_bonds_indices', i64[:,:]),
    ('Kh_list', f64[:]),
    ('mass_list', f64[:]),
    ('acc_list', f64[:,:]),
    ('eq_hinging_angles', f64[:]),
])
class Tripple_bond_force(object):
    def __init__(self, points, tripple_bonds_indices, Kh_list, mass_list, acc_list, eq_hinging_angles):
        self.points = points
        self.tripple_bonds_indices = tripple_bonds_indices
        self.Kh_list = Kh_list
        self.mass_list = mass_list
        self.acc_list = acc_list
        self.eq_hinging_angles = eq_hinging_angles
        
        self.acc_update()
        
        
    
    def acc_update(self):
        
        for i in range(len(self.tripple_bonds_indices)):
            
            id1 = self.tripple_bonds_indices[i][0]
            id2 = self.tripple_bonds_indices[i][1]
            id3 = self.tripple_bonds_indices[i][2]
            
            p1_vec = self.points[id1]
            p2_vec = self.points[id2]
            p3_vec = self.points[id3]
            
            r_12_vec = p1_vec - p2_vec
            r_32_vec = p3_vec - p2_vec
            
            #print r_12_vec, "r_12_vec"
            #print r_32_vec, "r_32_vec"
            #time.sleep(1000)
            
            r_12_vec_mag = ( r_12_vec[0]**2 + r_12_vec[1]**2 + r_12_vec[2]**2 )**0.5
            r_32_vec_mag = ( r_32_vec[0]**2 + r_32_vec[1]**2 + r_32_vec[2]**2 )**0.5
            
            Kh = self.Kh_list[i]
            
            angle_value = self.angle_123_calc(r_12_vec, r_32_vec, r_12_vec_mag, r_32_vec_mag)
            angle_eq = self.eq_hinging_angles[i]
            
            if angle_value != 0.0 and angle_value != np.pi:
            
                factor1 = - Kh * (angle_value - angle_eq)
                
                factor2 = self.force_second_factor(p1_vec, p2_vec, p3_vec, r_12_vec, r_32_vec, r_12_vec_mag, r_32_vec_mag)
                
                factor2_for_point1 = factor2[0]
                factor2_for_point2 = factor2[1]
                factor2_for_point3 = factor2[2]

                acc_change_p1 = factor1 * factor2_for_point1 / self.mass_list[id1]
                acc_change_p2 = factor1 * factor2_for_point2 / self.mass_list[id2]   
                acc_change_p3 = factor1 * factor2_for_point3 / self.mass_list[id3]
                
                self.acc_list[id1] += acc_change_p1
                self.acc_list[id2] += acc_change_p2            
                self.acc_list[id3] += acc_change_p3
            
            
            
    
    def angle_123_calc(self, r_12_vec, r_32_vec, r_12_vec_mag, r_32_vec_mag):
        
        angle = np.arccos( np.dot(r_12_vec, r_32_vec)/ (r_12_vec_mag * r_32_vec_mag) )
        
        #checking_if_everything_works = np.isnan(angle)
        #if checking_if_everything_works == 0:
            #None
        #else:
            #print r_12_vec, "r_12_vec"
            #print r_32_vec, "r_32_vec"
            #print r_12_vec_mag, "r_12_vec_mag"
            #print r_32_vec_mag, "r_32_vec_mag"
            #print np.dot(r_12_vec, r_32_vec), "np.dot(r_12_vec, r_32_vec)"
            #print np.dot(r_12_vec, r_32_vec)/ (r_12_vec_mag * r_32_vec_mag), "np.dot(r_12_vec, r_32_vec)/ (r_12_vec_mag * r_32_vec_mag)"
            #print np.arccos(np.around( np.dot(r_12_vec, r_32_vec)/ (r_12_vec_mag * r_32_vec_mag), 7 ) ), "np.arccos(np.dot(r_12_vec, r_32_vec)/ (r_12_vec_mag * r_32_vec_mag))"
            #time.sleep(1000)

        
        return angle
    
    
    
    def force_second_factor(self, p1_vec, p2_vec, p3_vec, r_12_vec, r_32_vec, r_12_vec_mag, r_32_vec_mag):
        
        r_1x = p1_vec[0] * 1.0
        r_1y = p1_vec[1] * 1.0
        r_1z = p1_vec[2] * 1.0

        r_2x = p2_vec[0] * 1.0
        r_2y = p2_vec[1] * 1.0
        r_2z = p2_vec[2] * 1.0
        
        r_3x = p3_vec[0] * 1.0
        r_3y = p3_vec[1] * 1.0
        r_3z = p3_vec[2] * 1.0
        
        r_12_mag = r_12_vec_mag
        r_32_mag = r_32_vec_mag
        #We try to calculate a complex derivatives d\theta / d r_12_vec etc.



        
        d_theta_d_r_1x = -((-1.0*r_1x + 1.0*r_2x)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5) + (-r_2x + r_3x)*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5))/np.sqrt(-((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))**2*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.0)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.0) + 1)
        
        d_theta_d_r_1y = -((-1.0*r_1y + 1.0*r_2y)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5) + (-r_2y + r_3y)*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5))/np.sqrt(-((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))**2*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.0)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.0) + 1)
        
        d_theta_d_r_1z = -((-1.0*r_1z + 1.0*r_2z)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5) + (-r_2z + r_3z)*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5))/np.sqrt(-((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))**2*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.0)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.0) + 1)

        derivative_d_theta_d_r_1 = np.array([d_theta_d_r_1x, d_theta_d_r_1y, d_theta_d_r_1z])

        d_theta_d_r_2x = -((1.0*r_1x - 1.0*r_2x)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5) + (-1.0*r_2x + 1.0*r_3x)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.5) + (-r_1x + 2*r_2x - r_3x)*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5))/np.sqrt(-((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))**2*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.0)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.0) + 1)
        
        d_theta_d_r_2y = -((1.0*r_1y - 1.0*r_2y)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5) + (-1.0*r_2y + 1.0*r_3y)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.5) + (-r_1y + 2*r_2y - r_3y)*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5))/np.sqrt(-((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))**2*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.0)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.0) + 1)
        
        d_theta_d_r_2z = -((1.0*r_1z - 1.0*r_2z)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5) + (-1.0*r_2z + 1.0*r_3z)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.5) + (-r_1z + 2*r_2z - r_3z)*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5))/np.sqrt(-((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))**2*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.0)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.0) + 1)
        
        derivative_d_theta_d_r_2 = np.array([d_theta_d_r_2x, d_theta_d_r_2y, d_theta_d_r_2z])
        
        d_theta_d_r_3x = -((r_1x - r_2x)*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5) + (1.0*r_2x - 1.0*r_3x)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.5))/np.sqrt(-((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))**2*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.0)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.0) + 1)
        
        d_theta_d_r_3y = -((r_1y - r_2y)*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5) + (1.0*r_2y - 1.0*r_3y)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.5))/np.sqrt(-((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))**2*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.0)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.0) + 1)
        
        d_theta_d_r_3z = -((r_1z - r_2z)*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-0.5) + (1.0*r_2z - 1.0*r_3z)*((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-0.5)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.5))/np.sqrt(-((r_1x - r_2x)*(-r_2x + r_3x) + (r_1y - r_2y)*(-r_2y + r_3y) + (r_1z - r_2z)*(-r_2z + r_3z))**2*((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**(-1.0)*((-r_2x + r_3x)**2 + (-r_2y + r_3y)**2 + (-r_2z + r_3z)**2)**(-1.0) + 1)
        
        derivative_d_theta_d_r_3 = np.array([d_theta_d_r_3x, d_theta_d_r_3y, d_theta_d_r_3z])
        
        return derivative_d_theta_d_r_1, derivative_d_theta_d_r_2, derivative_d_theta_d_r_3
    


    def results(self):
        return self.acc_list
        



# Testing:
if __name__ == '__main__':
    points = np.array([[0, 0, 0], [0.1, 1, 0], [-0.1, 1, 0]], np.float64)
    tripple_bonds_indices = np.array([[1, 0, 2]], np.int64)
    Kh_list = np.array([1.0], np.float64)
    mass_list = np.array([1, 1, 1], np.float64) * 20
    acc_list = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], np.float64)
    eq_hinging_angles = np.array([np.pi*0.5], np.float64)


    tripple_bond_force = Tripple_bond_force(points, tripple_bonds_indices, Kh_list, mass_list, acc_list, eq_hinging_angles)

    tripple_bond_force.acc_update()

    print(tripple_bond_force.acc_list)

    # visualize test case:
    VISUALIZE = 1
    if VISUALIZE:
        import pylab as plt
        for t in tripple_bonds_indices:
            plt.plot([points[t[i],0] for i in range(3)]
                    ,[points[t[i],1] for i in range(3)], "k-o")
        for i in range(len(points)):
            VEC = np.array([points[i], points[i] + acc_list[i]])
            plt.plot(VEC[:,0], VEC[:,1], "r-")
        plt.gca().set_aspect('equal')
        plt.show()
