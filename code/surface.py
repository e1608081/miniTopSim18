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
        
        #define start & enpoint of our surface parts
        self.x = np.array(self.x)
        self.y = np.array(self.y)
        s_x = self.x[:-1]
        s_y = self.y[:-1]
        e_x = self.x[1:]
        e_y = self.y[1:]
        
        #define a shorter variable name for point
        px = point[0]
        py = point[1]
        
        #calculate vector and normal vector of segment
        ab_x = e_x - s_x
        ab_y = e_y - s_y
        n_x = ab_y
        n_y = -ab_x
        
        #calculate the normal distance
        ap_x = px - s_x
        ap_y = py - s_y
        
        normal_distance = (n_x*ap_x + n_y*ap_y)/(n_x*n_x + n_y*n_y)
        
        #now calculate where we would be inside the segment
        #0 stands for start point
        #1 stands for end point
        segment_part = (ab_x*ap_x + ab_y*ap_y)/(ab_x*ab_x + ab_y*ab_y)
        
        #use this information to create a mask that shows validity of normal distance
        norm_d_valid = np.logical_and(segment_part >= 0, segment_part <= 1)
        
        #calculate minimum distance based on that
        minimum_distance = None
        if np.any(norm_d_valid):
            minimum_distance = (np.absolute(normal_distance[np.where(norm_d_valid)])).min()
        
        #get distance from start and enpoint and compare
        x_distance = s_x - px
        y_distance = s_y - py
        cur_distance = np.sqrt(x_distance*x_distance + y_distance*y_distance).min()
        if minimum_distance is None or cur_distance < minimum_distance:
            minimum_distance = cur_distance
        
        x_distance = e_x - px
        y_distance = e_y - py
        cur_distance = np.sqrt(x_distance*x_distance + y_distance*y_distance).min()
        if minimum_distance is None or cur_distance < minimum_distance:
            minimum_distance = cur_distance
        
        return minimum_distance
    
    def deloop(self):
        """Remove all loops from the surface"""
    
        n_segments = len(self.y)-1      #amount of segments
        
        #make 2D arrays containing all possible combinations for x,y and dx,dy
        x_i, x_j = np.meshgrid(self.x[:-1],self.x[:-1],indexing='ij')
        y_i, y_j = np.meshgrid(self.y[:-1],self.y[:-1],indexing='ij')
        dx=np.array(self.x[1:]-self.x[:-1])
        dy=np.array(self.y[1:]-self.y[:-1])
        dx_i, dx_j = np.meshgrid(dx, dx, indexing='ij')
        dy_i, dy_j = np.meshgrid(dy, dy, indexing='ij')
        
        
        #creating the arrays for the linear equations A*x=b
        #use np.stack to bring them in the right shape:
        #shape(A)=(M,M,2,2)
        #shape(b)=(M,M,2)
        b_temp = np.array((x_j-x_i,y_j-y_i))
        A_temp = np.array((np.stack((dx_i,-dx_j),axis=2),
                           np.stack((dy_i,-dy_j),axis=2)))
        b = np.stack(b_temp, axis=2)
        A = np.stack(A_temp, axis=2)
        
        #calculate the rank of each matrix in A.
        #If the rank is one the matrix is singular
        rank = np.linalg.matrix_rank(A)
        
        #create a mask, because not all equations in A*x=b
        #can or must be solved:
        #mask equations without a solution
        mask_singular = np.not_equal(rank,1)
        #mask the triangele matrix of (M,M), so we don't solve the same 
        #equation twice: ij = ji
        mask_sing_tri = np.triu(mask_singular) 
        #mask equations where i=j and j= i+1
        mask_diagonal = np.logical_not(np.logical_or(np.eye(n_segments, k=1),
                                                     np.eye(n_segments, k=0)))
        #combine all masks
        mask = np.logical_and(mask_sing_tri,mask_diagonal)
        
        #initiate the array for the solutions
        st=np.full_like(b,99999, dtype=float)
       
        #solve the equations A*st=b
        st[mask] = np.linalg.solve(A[mask],b[mask])
        #just solutions where: 0<= s,t <1 are relevant for crossings
        crossings = np.all(np.logical_and(np.greater_equal(st,0),
                                          np.less(st,1)), axis=2)

        #get the indexes of the intersecting segments        
        indexes_i, indexes_j = np.meshgrid(np.arange(n_segments),
                                           np.arange(n_segments),indexing='ij')
        i = indexes_i[crossings]
        j = indexes_j[crossings]
        #indexes has the shape: [[i_1,j_1], [i_2, j_2]]
        indexes = np.stack((i,j), axis =1)
        
        #create a mask to remove the loop-segments and add new surface-points
        #for the intersections
        mask_remove = np.ones(n_segments+1, dtype=bool)
        #remove the segments form the loop
        for temp_i, temp_j in indexes:
            mask_remove_temp = np.fromfunction(lambda k:np.logical_not(
                    np.logical_and(k>(temp_i+1),k<=temp_j)),(len(self.x),))
            mask_remove = np.logical_and(mask_remove,mask_remove_temp)
            self.x[temp_i+1] = self.x[temp_i]+ dx[temp_i]*st[(temp_i,temp_j,0)]
            self.y[temp_i+1] = self.y[temp_i]+ dy[temp_i]*st[(temp_i,temp_j,0)]
               
        self.x = self.x[mask_remove]
        self.y = self.y[mask_remove]   

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
