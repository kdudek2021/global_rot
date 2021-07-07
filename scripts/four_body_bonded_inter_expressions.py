import numpy as np
import time
import math
import sympy as sp

import matplotlib.pylab as plt

from matplotlib import collections  as mc



class Expressions_for_four_points():
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
        
        r_4x = sp.Symbol("r_4x")
        r_4y = sp.Symbol("r_4y")
        r_4z = sp.Symbol("r_4z")
        
        
        
        r_12x = r_1x - r_2x
        r_12y = r_1y - r_2y
        r_12z = r_1z - r_2z

        r_23x = r_2x - r_3x
        r_23y = r_2y - r_3y
        r_23z = r_2z - r_3z

        r_43x = r_4x - r_3x
        r_43y = r_4y - r_3y
        r_43z = r_4z - r_3z
        
        n_x = r_43y * r_23z - r_43z * r_23y
        n_y = r_43z * r_23x - r_43x * r_23z
        n_z = r_43x * r_23y - r_43y * r_23x
        
        self.expression_for_point1 = self.first_point_expression(n_x, n_y, n_z, r_12x, r_12y, r_12z, r_23x, r_23y, r_23z, r_43x, r_43y, r_43z, r_1x, r_1y, r_1z, r_2x, r_2y, r_2z, r_3x, r_3y, r_3z, r_4x, r_4y, r_4z)
        self.expression_for_point2 = self.second_point_expression(n_x, n_y, n_z, r_12x, r_12y, r_12z, r_23x, r_23y, r_23z, r_43x, r_43y, r_43z, r_1x, r_1y, r_1z, r_2x, r_2y, r_2z, r_3x, r_3y, r_3z, r_4x, r_4y, r_4z)        
        self.expression_for_point3 = self.third_point_expression(n_x, n_y, n_z, r_12x, r_12y, r_12z, r_23x, r_23y, r_23z, r_43x, r_43y, r_43z, r_1x, r_1y, r_1z, r_2x, r_2y, r_2z, r_3x, r_3y, r_3z, r_4x, r_4y, r_4z)
        self.expression_for_point4 = self.fourth_point_expression(n_x, n_y, n_z, r_12x, r_12y, r_12z, r_23x, r_23y, r_23z, r_43x, r_43y, r_43z, r_1x, r_1y, r_1z, r_2x, r_2y, r_2z, r_3x, r_3y, r_3z, r_4x, r_4y, r_4z)
        
        print(self.expression_for_point4[0], "  ||| term 3a_x |||")
        print(self.expression_for_point4[1], "  ||| term 3a_y |||")
        print(self.expression_for_point4[2], "  ||| term 3a_z |||")
        print(self.expression_for_point4[3], "  ||| term 3b_x |||")        
        print(self.expression_for_point4[4], "  ||| term 3b_y |||")
        print(self.expression_for_point4[5], "  ||| term 3b_z |||")



    def first_point_expression(self, n_x, n_y, n_z, r_12x, r_12y, r_12z, r_23x, r_23y, r_23z, r_43x, r_43y, r_43z, r_1x, r_1y, r_1z, r_2x, r_2y, r_2z, r_3x, r_3y, r_3z, r_4x, r_4y, r_4z):
        
        #term3a = ?
        
        expr_aux_3a = n_x * r_12x + n_y * r_12y + n_z * r_12z
        
        term3a_x = sp.diff(expr_aux_3a, r_1x)
        term3a_y = sp.diff(expr_aux_3a, r_1y)
        term3a_z = sp.diff(expr_aux_3a, r_1z)
        

        
        #term3b = ?
        
        term3b1_x = r_12y * r_23z - r_12z * r_23y
        term3b1_y = r_12z * r_23x - r_12x * r_23z
        term3b1_z = r_12x * r_23y - r_23x * r_12y
        
        term3b2_x = r_43y * r_23z - r_43z * r_23y
        term3b2_y = r_43z * r_23x - r_43x * r_23z
        term3b2_z = r_43x * r_23y - r_43y * r_23x
        
        term3b_aux = term3b1_x * term3b2_x + term3b1_y * term3b2_y + term3b1_z * term3b2_z
        
        term3b_x = sp.diff(term3b_aux, r_1x)
        term3b_y = sp.diff(term3b_aux, r_1y)
        term3b_z = sp.diff(term3b_aux, r_1z)
        
        return term3a_x, term3a_y, term3a_z, term3b_x, term3b_y, term3b_z
    
    
    
    def second_point_expression(self, n_x, n_y, n_z, r_12x, r_12y, r_12z, r_23x, r_23y, r_23z, r_43x, r_43y, r_43z, r_1x, r_1y, r_1z, r_2x, r_2y, r_2z, r_3x, r_3y, r_3z, r_4x, r_4y, r_4z):
        
        #term3a = ?
        
        expr_aux_3a = n_x * r_12x + n_y * r_12y + n_z * r_12z
        
        term3a_x = sp.diff(expr_aux_3a, r_2x)
        term3a_y = sp.diff(expr_aux_3a, r_2y)
        term3a_z = sp.diff(expr_aux_3a, r_2z)
        
        #term3b = ?
        
        term3b1_x = r_12y * r_23z - r_12z * r_23y
        term3b1_y = r_12z * r_23x - r_12x * r_23z
        term3b1_z = r_12x * r_23y - r_23x * r_12y
        
        term3b2_x = r_43y * r_23z - r_43z * r_23y
        term3b2_y = r_43z * r_23x - r_43x * r_23z
        term3b2_z = r_43x * r_23y - r_43y * r_23x
        
        term3b_aux = term3b1_x * term3b2_x + term3b1_y * term3b2_y + term3b1_z * term3b2_z
        
        term3b_x = sp.diff(term3b_aux, r_2x)
        term3b_y = sp.diff(term3b_aux, r_2y)
        term3b_z = sp.diff(term3b_aux, r_2z)

        return term3a_x, term3a_y, term3a_z, term3b_x, term3b_y, term3b_z
    
    
    def third_point_expression(self, n_x, n_y, n_z, r_12x, r_12y, r_12z, r_23x, r_23y, r_23z, r_43x, r_43y, r_43z, r_1x, r_1y, r_1z, r_2x, r_2y, r_2z, r_3x, r_3y, r_3z, r_4x, r_4y, r_4z):
        
        #term3a = ?
        
        expr_aux_3a = n_x * r_12x + n_y * r_12y + n_z * r_12z
        
        term3a_x = sp.diff(expr_aux_3a, r_3x)
        term3a_y = sp.diff(expr_aux_3a, r_3y)
        term3a_z = sp.diff(expr_aux_3a, r_3z)
        
        #term3b = ?
        
        term3b1_x = r_12y * r_23z - r_12z * r_23y
        term3b1_y = r_12z * r_23x - r_12x * r_23z
        term3b1_z = r_12x * r_23y - r_23x * r_12y
        
        term3b2_x = r_43y * r_23z - r_43z * r_23y
        term3b2_y = r_43z * r_23x - r_43x * r_23z
        term3b2_z = r_43x * r_23y - r_43y * r_23x
        
        term3b_aux = term3b1_x * term3b2_x + term3b1_y * term3b2_y + term3b1_z * term3b2_z
        
        term3b_x = sp.diff(term3b_aux, r_3x)
        term3b_y = sp.diff(term3b_aux, r_3y)
        term3b_z = sp.diff(term3b_aux, r_3z)

        return term3a_x, term3a_y, term3a_z, term3b_x, term3b_y, term3b_z
    
    
    def fourth_point_expression(self, n_x, n_y, n_z, r_12x, r_12y, r_12z, r_23x, r_23y, r_23z, r_43x, r_43y, r_43z, r_1x, r_1y, r_1z, r_2x, r_2y, r_2z, r_3x, r_3y, r_3z, r_4x, r_4y, r_4z):
        
        #term3a = ?
        
        expr_aux_3a = n_x * r_12x + n_y * r_12y + n_z * r_12z
        
        term3a_x = sp.diff(expr_aux_3a, r_4x)
        term3a_y = sp.diff(expr_aux_3a, r_4y)
        term3a_z = sp.diff(expr_aux_3a, r_4z)
        
        #term3b = ?
        
        term3b1_x = r_12y * r_23z - r_12z * r_23y
        term3b1_y = r_12z * r_23x - r_12x * r_23z
        term3b1_z = r_12x * r_23y - r_23x * r_12y
        
        term3b2_x = r_43y * r_23z - r_43z * r_23y
        term3b2_y = r_43z * r_23x - r_43x * r_23z
        term3b2_z = r_43x * r_23y - r_43y * r_23x
        
        term3b_aux = term3b1_x * term3b2_x + term3b1_y * term3b2_y + term3b1_z * term3b2_z
        
        term3b_x = sp.diff(term3b_aux, r_4x)
        term3b_y = sp.diff(term3b_aux, r_4y)
        term3b_z = sp.diff(term3b_aux, r_4z)
        
        return term3a_x, term3a_y, term3a_z, term3b_x, term3b_y, term3b_z
    
    
