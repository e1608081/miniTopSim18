# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 15:59:54 2018

@author: Jakob
"""

import os
import parameters as par
import numpy as np
from scipy import special, constants


class Beam:

    
    def __init__(self, type=None):
        """Initialize Beam from config file. Beam type can be overwritten.
        
        :param type: Beam type
        """
        global scannerinstance
        scannerinstance = self.scanning()
        if type is None:
            self.type = par.BEAM_TYPE
        else:
            self.type = type

        self.current = par.BEAM_CURRENT
        # wz and wx are input in nm
        self.wz = par.SCAN_WIDTH
        self.wx = par.ERF_BEAM_WIDTH
        self.xc = par.BEAM_CENTER
        self.sigma = self.__sigma_from_fwhm(par.FWHM)
        self.j = par.BEAM_CURRENT_DENSITY

    def __sigma_from_fwhm(self, fwhm):
        """Calculates standard deviation from fwhm.
        
        :param fwhm: fwhm
        
        :returns: standard deviation sigma
        """
        return fwhm / np.sqrt(8 * np.log(2))

    def __calculate_gauss(self, x, beam_center):
        """Calculates beam flux density based on gaussian distribution.

        :param x: position on x axis in nm
        
        :returns: beam flux density in Atoms/cm^2
        """
        exp_z = np.square(x - beam_center)
        exp_n = 2 * np.square(self.sigma)
        exp = np.exp(-(exp_z / exp_n))
        # calculations are done in nm, F_beam is expected in A/cm^2. nm^2 -> cm^2 conversion factor is 1e-14
        const = self.current / (np.sqrt(2 * np.pi) * self.sigma * self.wz * constants.e * 1e-14)
        f_beam = const * exp
        return f_beam

    def __calculate_erf(self, x, beam_center):
        """Calculates beam flux density based on error function.
        
        :param x: position on x axis in nm
        
        :returns: beam flux density in Atoms/cm^2
        """
        x1 = self.xc - (self.wx / 2)
        x2 = self.xc + (self.wx / 2)
        erf1 = (- (x - x2) / (np.sqrt(2) * self.sigma))
        erf2 = (- (x - x1) / (np.sqrt(2) * self.sigma))
        # calculations are done in nm, F_beam is expected in A/cm^2. nm^2 -> cm^2 conversion factor is 1e-14
        const = self.current / (2 * self.wx * self.wz * constants.e * 1e-14)
        erfcalc1 = special.erf(erf1)
        erfcalc2 = special.erf(erf2)
        f_beam = const * (erfcalc1 - erfcalc2)
        return f_beam

    def __call__(self, x):
        """Calculates beam flux density.
        
        :param x: position on x axis
        
        :returns: beam flux density
        """
        if self.type is 'Gaussian':
            try:
                calc = 0
                dwell_list = list()
                for beam_pos, dwell_time in scannerinstance:
                    
                    calc += self.__calculate_gauss(x,beam_pos)
                    dwell_list.append(dwell_time) ##new!!
                    #print("list", dwell_list)
                return calc
            except StopIteration:
                print("Check Parameters. Might be wrong!")
                
        elif self.type is 'error_function':
            try:
                calc = 0
                dwell_list = list()
                for beam_pos, dwell_time in scannerinstance:
                    calc += self.__calculate_gauss(x,beam_pos)
                    dwell_list.append(dwell_time)
                return calc
            except StopIteration:
                print("Check Parameters. Might be wrong!")
                
        else:
            # use constant "broad beam"
            # beam current density J = F_beam * e
            f_beam = self.j / constants.e
            return f_beam


    def scanning(self):
        """Switch between 3 possible scanning modes.
        
        returns: Pixels and dwell-time.
        """        
        scan_type = par.SCAN_TYPE
        dwell_time = par.DWELL_TIME
                
        if scan_type == 'None':
            """Beam stays on same spot.
            """
            while True:
                yield par.BEAM_CENTER, dwell_time       

        if scan_type == 'raster' or scan_type == 'serpentine':
            """Generates the Pixel position for raster and serpentine mode.
            
            In raster mode the beam is scanning across and gets back to the start.
            In serpentine mode the beam scans forward and back.
            """
            
            pixel_spacing = par.PIXEL_SPACING
            n_scans = par.N_SCANS
            start_position = par.BEAM_CENTER
            end_position = par.BEAM_CENTER + (par.N_PIXELS) * par.PIXEL_SPACING

            if scan_type == 'raster':
                
                for scans in range(int(n_scans)):
                    current_pixel = start_position
                    #jumps back after pass
                    
                    while(end_position - current_pixel) > 0.:
                        yield current_pixel, dwell_time 
                        
                        current_pixel += pixel_spacing 

            if scan_type == 'serpentine':
                
                current_pixel = start_position
                current_end = end_position
                current_start = start_position
                
                for scans in range(int(n_scans)):
                    while abs(current_end - current_pixel) > 0.:
                        yield current_pixel, dwell_time
                        
                        current_pixel += pixel_spacing
                    
#                    print("current_pixel",current_pixel)
#                    print("pixel_spacing", pixel_spacing)
#                    print("start_position",current_start)
#                    print("end_position", current_end)
#                    print("")
                    pixel_spacing *= (-1)
                    temp= current_pixel
                    current_pixel= current_end
                    current_end = current_start
                    current_start = temp

        if scan_type == 'stream file':
            """gets the position and dwell time from file
            
            File has to have more then one entry!
            """

            if os.path.isfile(par.STREAM_FILE):
                pixels, dwell_times = np.loadtxt(par.STREAM_FILE, unpack = True)
                index = 0
                
                for scans in range(int(par.N_SCANS)):
                    while index < pixels.size:
                        yield pixels[index], dwell_times[index]
                        index += 1
                    index = 0
                        
                    
            else: #no Stream file found
                print("Stream File not found!")
                