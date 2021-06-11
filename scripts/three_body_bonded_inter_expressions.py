import numpy as np
import time
import math
import sympy as sp

import matplotlib.pylab as plt

from matplotlib import collections  as mc



class Expressions_for_three_points():
    def __init__(self):
        
        self.determine_expressions()
        
        
    def determine_expressions(self):
        
        r_1x = sp.Symbol("r_1x")
        r_1y = sp.Symbol("r_1y")
        r_1z = sp.Symbol("r_1z")
        
        r_2x = sp.Symbol("r_2x")
        r_2y = sp.Symbol("r_2y")
        r_2z = sp.Symbol("r_2z")
        
        r_3x = sp.Symbol("r_3x")
        r_3y = sp.Symbol("r_3y")
        r_3z = sp.Symbol("r_3z")
        
        r_12_mag = ((r_1x - r_2x)**2 + (r_1y - r_2y)**2 + (r_1z - r_2z)**2)**0.5        #r_12_mag = sp.Symbol("r_12_mag")
        r_32_mag = ((r_3x - r_2x)**2 + (r_3y - r_2y)**2 + (r_3z - r_2z)**2)**0.5        #r_32_mag = sp.Symbol("r_32_mag")
        
        r_12x = r_1x - r_2x
        r_12y = r_1y - r_2y
        r_12z = r_1z - r_2z

        r_32x = r_3x - r_2x
        r_32y = r_3y - r_2y
        r_32z = r_3z - r_2z
        
        theta = sp.acos((r_12x * r_32x + r_12y * r_32y + r_12z * r_32z) / (r_12_mag * r_32_mag))
        
        force1_expr_x = sp.diff(theta, r_1x)
        force1_expr_y = sp.diff(theta, r_1y)
        force1_expr_z = sp.diff(theta, r_1z)
        
        print(force1_expr_x, "  force1_expr_x  ")
        print(force1_expr_y, "  force1_expr_y  ")
        print(force1_expr_z, "  force1_expr_z  ")
        print(" |||||||||||||||||||||||||||||||||||||")
    
        force2_expr_x = sp.diff(theta, r_2x)
        force2_expr_y = sp.diff(theta, r_2y)
        force2_expr_z = sp.diff(theta, r_2z)

        print(force2_expr_x, "  force2_expr_x  ")
        print(force2_expr_y, "  force2_expr_y  ")
        print(force2_expr_z, "  force2_expr_z  ")
        print(" |||||||||||||||||||||||||||||||||||||")
        
        force3_expr_x = sp.diff(theta, r_3x)
        force3_expr_y = sp.diff(theta, r_3y)
        force3_expr_z = sp.diff(theta, r_3z)
        
        print(force3_expr_x, "  force3_expr_x  ")
        print(force3_expr_y, "  force3_expr_y  ")
        print(force3_expr_z, "  force3_expr_z  ")
        print(" |||||||||||||||||||||||||||||||||||||")