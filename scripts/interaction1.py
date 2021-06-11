import numpy as np

from numba import jitclass
from numba import float64 as f64
from numba import int64 as i64


@jitclass([
    ('points', f64[:,:]),
    ('double_bonds_indices', i64[:,:]),
    ('k_constant_list', f64[:]),
    ('eq_lengths', f64[:]),
    ('mass_list', f64[:]),
    ('acc_list', f64[:,:]),
])
class Double_bond_force(object):
    def __init__(self, points, double_bonds_indices, mass_list, acc_list, k_constant_list, eq_lengths):
        self.points = points
        self.double_bonds_indices = double_bonds_indices
        self.k_constant_list = k_constant_list
        self.eq_lengths = eq_lengths
        self.mass_list = mass_list
        self.acc_list = acc_list
        
        self.acceleration_update()
        
        
    def acceleration_update(self):
        
        for i in range(len(self.k_constant_list)):
            k_stiff = self.k_constant_list[i]
            
            id1 = self.double_bonds_indices[i][0]#index of the first point
            id2 = self.double_bonds_indices[i][1]
            
            r_vec_12 = self.points[id1] - self.points[id2]
            r_vec_21 = - r_vec_12
            
            r_21_dist = ( (self.points[id1][0] - self.points[id2][0])**2 + (self.points[id1][1] - self.points[id2][1])**2 + (self.points[id1][2] - self.points[id2][2])**2 )**0.5
            r_21_eq = self.eq_lengths[i]
            
            force_1 = - 2.0 * k_stiff * (r_21_dist - r_21_eq) * (r_vec_12 / r_21_dist)
            force_2 = - force_1
            
            a1_change = force_1 / self.mass_list[id1]
            a2_change = force_2 / self.mass_list[id2]
            
            self.acc_list[id1] += a1_change
            self.acc_list[id2] += a2_change
    
    
    def results(self):
        return self.acc_list


if __name__ == '__main__':
    points = np.array([[0, 0, 0], [0.1, 1, 0], [-0.1, 1, 0]], np.float64)
    double_bonds_indices = np.array([[1, 0], [0, 2]], np.int64)
    k_constant_list = np.array([1.0, 1.0], np.float64)
    eq_lengths = np.array([1, 1], np.float64) * 1.01
    mass_list = np.array([1, 1, 1], np.float64) * 0.1
    acc_list = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]], np.float64)


    double_bond_force = Double_bond_force(
            points, double_bonds_indices, mass_list, acc_list, k_constant_list,
            eq_lengths)

    double_bond_force.acceleration_update()

    print(double_bond_force.acc_list)

    # visualize test case:
    VISUALIZE = 1
    if VISUALIZE:
        import pylab as plt
        for t in double_bonds_indices:
            plt.plot([points[t[i],0] for i in range(2)]
                    ,[points[t[i],1] for i in range(2)], "k-o")
        for i in range(len(points)):
            VEC = np.array([points[i], points[i] + acc_list[i]])
            plt.plot(VEC[:,0], VEC[:,1], "r-")
        plt.gca().set_aspect('equal')
        plt.show()
