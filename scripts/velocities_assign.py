import numpy as np
import time
import math

import matplotlib.pylab as plt

from matplotlib import collections  as mc



class Velocities_assignment():
    def __init__(self, points):
        self.points = points
        
        self.assignment()
        
        
    def assignment(self):
        
        self.velocities = []
        
        for i in range(len(self.points)):
            self.velocities.append(np.array([0.0, 0.0, 0.0]))
    
    
    def results(self):
        return self.velocities
