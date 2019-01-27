# -*- coding: utf-8 -*-
"""
Surface object and functions for miniTopSim

Class functions:
    normal(): calculates the normal vector and returns it
    plot(tend, dtime): plots the start and given values
    write(file, time): writes output file
"""

import sys
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
from copy import deepcopy

from init_surface import init_surface
import parameters as par

class Surface:
    def __init__(self, time = None, xValues = None, yValues = None):
        """Init the object. Copy start value for printing."""

        if time is None or xValues is None or yValues is None:
            # create a new object based on parameter database

            # use the surface file if specified
            if par.INITIAL_SURFACE_FILE != '':
                surface = load(par.INITIAL_SURFACE_FILE)
                self.x = surface.x
                self.y = surface.y
                self.startx = deepcopy(self.x)
                self.starty = deepcopy(self.y)

            # if no surface file is specified, use the chosen surface type
            else:
                num_points = (par.XMAX - par.XMIN) // par.DELTA_X + 1
                self.x = np.linspace(int(par.XMIN), int(par.XMAX), int(num_points))
                self.y = np.array(init_surface(self.x))
                self.startx = deepcopy(self.x)
                self.starty = deepcopy(self.y)

        else:
            #initial values given use them
            self.x = np.array(xValues)
            self.y = np.array(yValues)
            self.startx = deepcopy(self.x)
            self.starty = deepcopy(self.y)

    def normal(self):
        """Calc normal vector and return it."""
        x = np.concatenate(([self.x[0]], self.x, [self.x[-1]]))
        y = np.concatenate(([self.y[0]], self.y, [self.y[-1]]))

        # Find duplicate nodes. Mask indices match node indices on [1:-1]
        duplicate_mask = np.logical_and(np.equal(x[:-1], x[1:]), np.equal(y[:-1], y[1:]))
        dm_left = np.append(duplicate_mask, False)
        dm_right = np.insert(duplicate_mask, 0, False)
        duplicate_mask = np.logical_or(dm_left, dm_right)

        duplicate_mask[-2] = False
        duplicate_mask[1] = False
        unique_mask = np.invert(duplicate_mask)
        duplicate_mask[0] = False
        duplicate_mask[-1] = False
        duplicate_index, = np.nonzero(duplicate_mask)
        unique_index, = np.nonzero(unique_mask)

        # Calculate normal vectors for unique nodes.
        dx = np.zeros_like(x)
        dy = np.zeros_like(y)
        dx[unique_index] = x[unique_index+1] - x[unique_index-1]
        dy[unique_index] = y[unique_index+1] - y[unique_index-1]

        length = np.linalg.norm([dx, dy], axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            norm_x = dy / length
            norm_y = -dx / length

        # Calculate angle of existing normal vectors.
        norm_angles = np.angle(norm_x + 1j*norm_y)
        angle_index = []

        # Traverse through duplicate node mask to get indices and amounts of duplicate nodes.
        dpm = np.array(duplicate_mask)
        for i in range(len(dpm)-1):
            index = deepcopy(i)
            while dpm[i]:
                dpm[i] = False
                i = i+1
            angle_index.append((index, i-index))

        # Calculate angles for duplicate nodes.
        for index, elem in angle_index:
            new_angles = np.linspace(norm_angles[index-1], norm_angles[index+elem],
                                     elem+2, endpoint=True)
            norm_angles[index:index + elem] = new_angles[1:-1]

        norm_x[duplicate_index] = np.cos(norm_angles[duplicate_index])
        norm_y[duplicate_index] = np.sin(norm_angles[duplicate_index])

        return norm_x[1:-1], norm_y[1:-1]
               
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

    def eliminate_overhangs(self):
        """Eliminate overhanging parts of the surface.
        """
        for i in range(self.x.size-1):
            self.x[i+1:] = np.where(
                    np.logical_and(self.x[i+1:] < self.x[i],
                                   self.y[i+1:] < self.y[i]),
                    self.x[i],
                    self.x[i+1:])
            self.x[:i] = np.where(
                    np.logical_and(self.x[:i] > self.x[i],
                                   self.y[:i] < self.y[i]),
                    self.x[i],
                    self.x[:i])

    def calc_viewfactor(self):
        """Calculates the view-factor from surface parameters

        :return: nxn matrix representing the view-factor
        """
        # create arrays for nodes i & j (j has equal rows, i has equal columns)
        x_j = np.ones(shape=(self.x.size, self.x.size)) * self.x
        y_j = np.ones(shape=(self.y.size, self.y.size)) * self.y

        x_i = x_j.transpose()
        y_i = y_j.transpose()

        # print("x_i: {}".format(x_i[0]))
        # print("x_j: {}".format(x_j[0]))

        # calculate distances between nodes i & j
        x_ij = x_i - x_j
        y_ij = y_i - y_j

        # create arrays for the normal vectors of nodes i & j (j has equal rows, i has equal columns)
        x_normal, y_normal = self.normal()
        x_i_normal = np.ones_like(x_i) * -x_normal[:, np.newaxis]
        y_i_normal = np.ones_like(y_i) * -y_normal[:, np.newaxis]
        x_j_normal = np.ones_like(x_j) * -x_normal
        y_j_normal = np.ones_like(y_j) * -y_normal

        # calculate cosines of angles (cos_a, cos_b) with [(a*b)/(|a|*|b|)]
        # deactivate warnings when division by zero -> is covered with np.nan_to_num
        with np.errstate(divide='ignore', invalid='ignore'):
            cos_a = np.nan_to_num((x_j_normal*x_ij + y_j_normal*y_ij) / (np.sqrt(x_j_normal**2+y_j_normal**2) * np.sqrt(x_ij**2+y_ij**2)))
            cos_b = np.nan_to_num((x_i_normal*-x_ij + y_i_normal*-y_ij) / (np.sqrt(x_i_normal**2+y_i_normal**2) * np.sqrt(x_ij**2+y_ij**2)))

        # calculate distances between nodes i & j (d_ij)
        d_ij = np.sqrt(x_ij ** 2 + y_ij ** 2)

        # calculate surface length of node (delta_l)
        x = np.concatenate(([2*self.x[0]-self.x[1]], self.x, [3*self.x[-1]-self.x[-2]]))
        y = np.concatenate(([self.y[0]], self.y, [self.y[-1]]))
        d_x = x[1:] - x[:-1]
        d_y = y[1:] - y[:-1]
        d_l = np.sqrt(d_x**2 + d_y**2)
        d_l_avg = (d_l[1:]+d_l[:-1]) / 2
        delta_l = np.ones_like(x_i) * d_l_avg

        # calculate view-factor and mask out all values where cos_a<0 and cos_b<0 and all elements on diagonal
        # deactivate warnings when division by zero ->is covered with np.nan_to_num
        with np.errstate(divide='ignore', invalid='ignore'):
            f_ij = np.nan_to_num((cos_a * cos_b * delta_l) / (2 * d_ij))
        mask = (cos_a > np.zeros_like(x_i)) * (cos_b > np.zeros_like(x_i)) * np.invert(np.eye(self.x.size, dtype=bool))
        return f_ij * mask

    def calc_viewfactor_derivative(self):
        """Calculates the derivative of the view-factor from surface parameters

        :return: nxn matrix representing the derivative of the view-factor
        """
        # create arrays for nodes i & j (j has equal rows, i has equal columns)
        x_j = np.ones(shape=(self.x.size, self.x.size)) * self.x
        y_j = np.ones(shape=(self.y.size, self.y.size)) * self.y

        x_i = x_j.transpose()
        y_i = y_j.transpose()

        # print("x_i: {}".format(x_i[0]))
        # print("x_j: {}".format(x_j[0]))

        # calculate distances between nodes i & j
        x_ij = x_i - x_j
        y_ij = y_i - y_j

        # create arrays for the normal vectors of nodes i & j (j has equal rows, i has equal columns)
        x_normal, y_normal = self.normal()
        x_i_normal = np.ones_like(x_i) * -x_normal[:, np.newaxis]
        y_i_normal = np.ones_like(y_i) * -y_normal[:, np.newaxis]
        x_j_normal = np.ones_like(x_j) * -x_normal
        y_j_normal = np.ones_like(y_j) * -y_normal

        # calculate cosines of angles (cos_a, cos_b) with [(a*b)/(|a|*|b|)]
        # deactivate warnings when division by zero -> is covered with np.nan_to_num
        with np.errstate(divide='ignore', invalid='ignore'):
            cos_a = np.nan_to_num((x_j_normal*x_ij + y_j_normal*y_ij) / (np.sqrt(x_j_normal**2+y_j_normal**2) * np.sqrt(x_ij**2+y_ij**2)))
            cos_b = np.nan_to_num((x_i_normal*-x_ij + y_i_normal*-y_ij) / (np.sqrt(x_i_normal**2+y_i_normal**2) * np.sqrt(x_ij**2+y_ij**2)))
            sin_b = ((x_i_normal * x_ij + y_i_normal * y_ij) / (
            np.sqrt(x_i_normal ** 2 + y_i_normal ** 2) * np.sqrt(x_ij ** 2 + y_ij ** 2)))
            sin_b = np.sqrt(1 - sin_b**2)
            sin_b = np.nan_to_num(sin_b)
            sin_b = np.copysign(sin_b, x_normal)

        # calculate distances between nodes i & j (d_ij)
        d_ij = np.sqrt(x_ij ** 2 + y_ij ** 2)

        # calculate surface length of node (delta_l)
        x = np.concatenate(([2*self.x[0]-self.x[1]], self.x, [3*self.x[-1]-self.x[-2]]))
        y = np.concatenate(([self.y[0]], self.y, [self.y[-1]]))
        d_x = x[1:] - x[:-1]
        d_y = y[1:] - y[:-1]
        d_l = np.sqrt(d_x**2 + d_y**2)
        d_l_avg = (d_l[1:]+d_l[:-1]) / 2
        delta_l = np.ones_like(x_i) * d_l_avg

        # calculate view-factor and mask out all values where cos_a<0 and cos_b<0 and all elements on diagonal
        # deactivate warnings when division by zero ->is covered with np.nan_to_num
        with np.errstate(divide='ignore', invalid='ignore'):
            f_ij = np.nan_to_num((cos_a * sin_b * delta_l) / (2 * d_ij))
        mask = (cos_a > np.zeros_like(x_i)) * (cos_b > np.zeros_like(x_i)) * np.invert(np.eye(self.x.size, dtype=bool))
        return f_ij * mask

    def adapt(self):
        """"Adds or removes nodes to keep surface length and normal vector angle near constant
        """
        self._delete_nodes()
        try:
            self._insert_by_distance()
        except ValueError as Err:
            raise Err
        self._insert_by_angle()

    def _insert_by_distance(self, interpol_type='quadratic'):
        """Inserts nodes if distance exceeds par.MAX_SEGLEN"""
        # Raise Error if points not in strict ascending order.
        if np.any(self.x[1:] <= self.x[:-1]):
            index = tuple(np.where(self.x[1:] <= self.x[:-1])[0])
            # self.x = np.delete(self.x, index)
            # self.y = np.delete(self.y, index)
            raise ValueError('Adaptive Grid: Nodes not in strictly ascending order at indices:', index)

        # Calculate distances of nodes
        x_dist = np.abs(self.x[1:] - self.x[:-1])
        y_dist = np.abs(self.y[1:] - self.y[:-1])
        distance = np.insert(np.sqrt(x_dist ** 2 + y_dist ** 2), 0, 0)

        # Find where and how many nodes to insert by distance criteria
        dist_index, = np.where(distance > par.MAX_SEGLEN)
        # insert_nodes = np.floor_divide(distance[dist_index], par.MAX_SEGLEN)
        insert_nodes = np.int32(np.ceil(distance[dist_index] / par.MAX_SEGLEN))
        new_points = np.split(self.x, dist_index)

        # Insert equally spaced points in x-axis
        for i in range(len(dist_index)):
            line = np.linspace(self.x[dist_index[i]-1],
                               self.x[dist_index[i]],
                               insert_nodes[i] + 1,
                               endpoint=True)
            line = line[1:-1]
            new_points[i] = np.append(new_points[i], line)

        interpolate = scipy.interpolate.interp1d(self.x, self.y, interpol_type)
        self.x = np.concatenate(new_points)
        self.y = interpolate(self.x)

    def _insert_by_angle(self):
        surface_x = self.x[1:] - self.x[:-1]
        surface_y = self.y[1:] - self.y[:-1]
        surface_norms = np.angle(surface_y - 1j * surface_x, deg=True)
        surface_angles = np.abs(surface_norms[1:] - surface_norms[:-1])

        # Find where to insert nodes by angle criteria.
        angle_index, = np.where(surface_angles > par.MAX_ANGLE)

        # Split before Nodes to duplicate.
        x = np.split(self.x, angle_index + 1)
        y = np.split(self.y, angle_index + 1)

        # Insert nodes.
        for i in range(len(angle_index)):
            amt = np.int32(np.floor(np.abs(surface_angles[angle_index[i]] / par.MAX_ANGLE)))-1

            x[i] = np.concatenate((x[i], np.full(amt, x[i+1][0])))
            y[i] = np.concatenate((y[i], np.full(amt, y[i+1][0])))

        self.x = np.concatenate(x)
        self.y = np.concatenate(y)

    def _delete_nodes(self):
        x = np.concatenate(([self.x[0]], self.x, [self.x[-1]]))
        y = np.concatenate(([self.y[0]], self.y, [self.y[-1]]))

        # Calculate where node removal meets angle criteria.
        surface_x = x[1:] - x[:-1]
        surface_y = y[1:] - y[:-1]
        surface_norms = np.angle(surface_y - 1j * surface_x, deg=True)
        surface_angles = np.abs(surface_norms[1:] - surface_norms[:-1])
        angle_mask = surface_angles < par.MAX_ANGLE

        # Calculate where node removal meets distance criteria.
        x_dist2 = x[2:] - x[:-2]
        y_dist2 = y[2:] - y[:-2]
        dist2 = np.sqrt(x_dist2**2 + y_dist2**2)

        distance_mask = dist2 < par.MAX_SEGLEN
        distance_mask[0] = False
        distance_mask[-1] = False
        # Get indices where distance criteria is met.
        index, = np.where(distance_mask)
        rem = np.vstack((index, dist2[index]))

        # Remove nodes with smallest new distance first, to prevent chain removal of adjacent nodes.
        while rem.size:
            minimum = np.argmin(rem, axis=1)
            index = int(minimum[1])
            # Recalculate distance criteria for preceding node.
            if rem[0][index] - 1 in rem[0]:
                j = int(rem[0][index])

                # Calculate index offset.
                off_l = 2
                while distance_mask[index-off_l] and index-off_l not in rem[0]:
                    off_l = off_l + 1
                off_r = 1
                while distance_mask[index + off_r] and index+off_r not in rem[0]:
                    off_r = off_r + 1

                # Update distance in preceding node.
                x_d = np.abs(x[j + off_r] - x[j - off_l])
                y_d = np.abs(y[j + off_r] - y[j - off_l])
                rem[1][index - 1] = np.sqrt(x_d ** 2 + y_d ** 2)

            # Recalculate distance criteria for following node.
            if rem[0][index] + 1 in rem[0]:
                j = int(rem[0][index])

                # Calculate index offset.
                off_l = 1
                while distance_mask[index-off_l] and index-off_l not in rem[0]:
                    off_l = off_l + 1
                off_r = 2
                while distance_mask[index + off_r] and index+off_r not in rem[0]:
                    off_r = off_r + 1

                # Update distance in preceding node.
                x_d = np.abs(x[j + off_r] - x[j - off_l])
                y_d = np.abs(y[j + off_r] - y[j - off_l])
                rem[1][index - 1] = np.sqrt(x_d ** 2 + y_d ** 2)

            # Mark current node for deletion.
            if rem[1][index] < par.MAX_SEGLEN:
                rem = np.delete(rem, index, 1)
            else:
                distance_mask[int(rem[0][index])] = False
                rem = np.delete(rem, index, 1)

        # Remove nodes.
        delete, = np.where(np.logical_and(angle_mask, distance_mask))

        self.x = np.delete(self.x, delete)
        self.y = np.delete(self.y, delete)

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
