# -*- coding: utf-8 -*-
"""
Surface object and functions for miniTopSim

Class functions:
    normal(): calculates the normal vector and returns it
    plot(tend, dtime): plots the start and given values
    write(file, time): writes output file
"""


import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

from init_surface import init_surface
import parameters as par

class Surface:
    
    def __init__(self, time = None, xValues = None, yValues = None):
        """Init the object. Copy start value for printing."""
        
        if time is None or xValues is None or yValues is None:
            #create a new object based on parameter database
            num_points = (par.XMAX - par.XMIN) // par.DELTA_X + 1
            self.x = np.linspace(int(par.XMIN), int(par.XMAX), int(num_points))
            self.y = np.array(init_surface(self.x))
            self.startx = deepcopy(self.x)
            self.starty = deepcopy(self.y)
        else:
            #initial values given use them
            self.x = deepcopy(xValues)
            self.y = deepcopy(yValues)
            self.startx = deepcopy(self.x)
            self.starty = deepcopy(self.y)
        
    def normal(self):
        """Calc normal vector and return it."""        
        x = np.concatenate(([self.x[0]], self.x, [self.x[-1]]))
        y = np.concatenate(([self.y[0]], self.y, [self.y[-1]]))
        dx = self._calc_vector(x)
        dy = self._calc_vector(y)
        length = np.linalg.norm([dx,dy],axis=0)
        return dy/length, -dx/length
    
    def _calc_vector(self, value):
        """Subtract end coordinate with start coordinate.
        
        :param value: coordinate array
        """
        delta = value[2:] - value[:-2]
        return delta
               
    def write(self, file, time, mode):  
        """Write output file.
        
        :param file: filename
        :param time: current time
        """
        with open(file, mode) as fp:
            fp.write('surface: {} {} x-positions y-positions\n'.format(time, len(self.x)))
            for x, y in zip(self.x, self.y):
                fp.write("{} {}\n".format(x, y))
    
    def distance(self, other):
        """Calculates the distance of two surfaces.
        The distance is the mean of the distance of one surface to the points of the other
        and vice versa.
        """
        return (self.points_distance(list(zip(other.x, other.y))) +
                other.points_distance(list(zip(self.x, self.y)))) / 2
    
    def points_distance(self, points):
        """Calculates the distance of a given list of points from this surface"""
        sum = 0
        for point in points:
            sum+= self.point_distance(point)
        
        return sum / len(points)
    
    def point_distance(self, point):
        """Calculates the distance of a point to a surface"""
        minimum_distance = None
        cur_distance = None
        x_distance = None
        y_distance = None
        
        own_points = np.array(list(zip(self.x, self.y)))
        
        """
        Calculating distance of segment with point
        assume endpoints are A = (a,b) and B = (c,d)
        assume point is P = (x,y)
        
        We calculate:
        surface Vector is AB = (a-c, b-d)
        normal Vector is n = (d-b, a-c)
        Vector from endpoint to Point is AP = (a-x, b-y)
        dot product of AP and n is our wanted distance
        -> AP * n = (a-x)*(d-b) + (b-y)*(a-c) =
                  = a*d - b*c + x*(b-d) + y*(c-a)
        """
        s_points = own_points[:-1]
        e_points = own_points[1:]
        
        cur_distance = s_points[:, 0]*e_points[:, 1] - e_points[:, 0]*s_points[:, 1] + \
                point[0]*(s_points[:, 1] - e_points[:, 1]) + \
                point[1]*(e_points[:, 0] - s_points[:, 0])
        
        minimum_distance = (np.absolute(cur_distance)).min()
        
        #now do it with start and endpoint
        x_distance = own_points[0][0] - point[0]
        y_distance = own_points[0][1] - point[1]
        cur_distance = np.sqrt(x_distance*x_distance + y_distance*y_distance)
        if cur_distance < minimum_distance:
            minimum_distance = cur_distance
        
        x_distance = own_points[-1][0] - point[0]
        y_distance = own_points[-1][1] - point[1]
        cur_distance = np.sqrt(x_distance*x_distance + y_distance*y_distance)
        if cur_distance < minimum_distance:
            minimum_distance = cur_distance
        
        return minimum_distance
    
    def deloop(self):
        """Remove all loops from the surface"""
        
        x=np.array((self.x[:-1], self.x[1:]-self.x[:-1]))   #(x,dx)
        y=np.array((self.y[:-1], self.y[1:]-self.y[:-1]))   #(y,dy)
        
        for i in range(len(self.x)-1):
            #start from i+2, to skip neighbour segment
            for j in np.arange(i+2,len(self.x)-1):
                
                #coefficient matrix
                A = np.array(((x[1,i], -x[1,j]),(y[1,i], -y[1,j])))
                #inhomogeneity
                b = np.array((x[0,j]-x[0,i],y[0,j]-y[0,i]))
                              
                try:
                    #the linear equasion has exactly one solution
                    st = np.linalg.solve(A,b)   
                                   
                    #0 <= s,t < 1
                    if np.all(np.logical_and(st>=0, st<1)):
                        
                        #create a mask which deletes the survace points
                        #from the loop
                        mask_cut = np.fromfunction(lambda k: 
                            np.logical_not(np.logical_and(k>(i+1),k<=j)),
                            (len(self.x),))
                        
                        self.x[i+1] = x[0,i]+ x[1,i]*st[0]
                        self.y[i+1] = y[0,i]+ y[1,i]*st[0]
                        self.x =self.x[mask_cut]
                        self.y =self.y[mask_cut]
                        x=np.array((self.x[:-1], self.x[1:]-self.x[:-1]))   
                        y=np.array((self.y[:-1], self.y[1:]-self.y[:-1]))   
                        
                        break
                    
                except np.linalg.linalg.LinAlgError: 
                    #the linear equasion has more then one or no solution
                    None

def load(file, wanted_time = None):
    """
    Loads the data from a .srt file.
    Loads the nearest time to that specified. If None given it will use last possible time
    
    :param file filename to load the data from
    :param time time to load from file
    
    :return return a surface object with best matching time inside file
    """
    time = -1
    numel = None
    xValues = []
    yValues = []
    
    #get relevant data from save file
    with open(file) as data_file:
        for line in data_file:
            # Check where new surface begins and get time and number of elements
            if 'surface:' in line:
                splittedLine = line.split(' ')
                cur_time = float(splittedLine[1])
                
                #we assume the time to be sorted so we can simply wait till the diff gets bigger
                if wanted_time is not None and abs(cur_time - wanted_time) > abs(time - wanted_time):
                        #best found simply stop loop
                        break
                
                time = cur_time
                numel = float(splittedLine[2])
                
                # Reset lists to be empty
                xValues.clear()
                yValues.clear()
            else:
                # Extract xy pair from line
                xyValuePair = line.rstrip('\n').split(' ')
                xValue = xyValuePair[0]
                yValue = xyValuePair[1]
                xValues.append(float(xValue))
                yValues.append(float(yValue))
    
    #transform it into a surface object
    return Surface(time, xValues, yValues)
