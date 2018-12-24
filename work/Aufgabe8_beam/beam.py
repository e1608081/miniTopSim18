# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 15:59:54 2018

@author: Jakob
"""

import parameters as par
import numpy as np
from scipy import special




class Beam:

    def __init__(self, type=None):
        """Initialize Beam from config file. Beam type can be overwritten.
        
        :param type: Beam type
        """
        self.current = par.BEAM_CURRENT
        self.wz = par.SCAN_WIDTH
        self.xc = par.BEAM_CENTER
        self.sigma = self._sigma_from_fwhm(par.FWHM)
        self.wx = par.ERF_BEAM_WIDTH
        if type is None:
            self.type = par.BEAM_TYPE
        else:
            self.type = type
            
            
    def _sigma_from_fwhm(self, fwhm):
        """Calculates standard deviation from fwhm.
        
        :param fwhm: fwhm
        
        :returns: standard deviation sigma
        """
        return fwhm / np.sqrt(8 * np.log(2))
    
    
    def _calculate_gauss(self, x):
        """Calculates beam flux density based on gaussian distribution.
        
        :param x: position on x axis
        
        :returns: beam flux density
        """
        exp_z = np.square(x - self.xc)
        exp_n = 2 * np.square(self.sigma)
        exp = ( - (exp_z / exp_n))
        const = self.current / (np.sqrt(2 * np.pi) * self.sigma * self.wz * np.e)
        f_beam = const * exp
        return f_beam
    
    
    def _calculate_erf(self, x):
        """Calculates beam flux density based on error function.
        
        :param x: position on x axis
        
        :returns: beam flux density
        """
        x1 = self.xc - (self.wx / 2)
        x2 = self.xc + (self.wx / 2)
        erf1 = np.exp( - (np.square(x - x2) / (2 * np.square(self.sigma))))
        erf2 = np.exp( - (np.square(x - x1) / (2 * np.square(self.sigma))))
        const = self.current / (2 * self.wx * self.wz)
        f_beam = const * (special.erf(erf2) - special.erf(erf1))
        return f_beam
       
    
    def __call__(self, x):
        """Calculates beam flux density.
        
        :param x: position on x axis
        
        :returns: beam flux density
        """
        if self.type is 'Gaussian':
            return self._calculate_gauss(x)    
        elif self.type is 'error function':
            return self._calculate_erf(x)
        else:
            #use constant "broad beam" 
            # beam current density J = F_beam * e
            j = par.BEAM_CURRENT_DENSITY
            f_beam = j / np.e
            return f_beam
        
    

