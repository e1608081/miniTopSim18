# -*- coding: utf-8 -*-
"""
USAGE: $ python3 process_parameters.py [configfile]

Test the configfile and write output if conformed

Configfiles:
    test1.cfg -> Replace float by Total_Time -> Correct
    test2.cfg -> Undefined Parameter -> Fail
    test3.cfg -> Datatypes don't match -> Fail
    test4.cfg -> No value -> Fail
    test5.cfg -> Invalid Condition -> Fail
    test6.cfg -> Several failurs -> Fail
    test7.cfg -> All Parameters new defined -> Correct
    
Created on Sun Nov 25 10:28:49 2018
@author: Christian Eder
"""

import sys
import os
import configparser

import parameters as par

filedir = os.path.dirname(__file__)
codedir = os.path.join(filedir, '..', '..', 'code')
sys.path.insert(0, codedir)

def _write_output(output_file):
    """Write the values to output file
    
    :param output_file: file for output data
    """
    cfg = configparser.ConfigParser()

    section = 'Setup'
    cfg.add_section(section)
    cfg.set(section, 'ETCHING', str(par.ETCHING))
    cfg.set(section, 'REDEP', str(par.REDEP))
    
    section = 'Initial Conditions'
    cfg.add_section(section)
    cfg.set(section, 'XMIN', str(par.XMIN))
    cfg.set(section, 'XMAX', str(par.XMAX))
    cfg.set(section, 'DELTA_X', str(par.DELTA_X))
    cfg.set(section, 'INITIAL_SURFACE_TYPE', str(par.INITIAL_SURFACE_TYPE))
    cfg.set(section, 'FUN_XMIN', str(par.FUN_XMIN))
    cfg.set(section, 'FUN_XMAX', str(par.FUN_XMAX))
    
    section = 'Numerics'
    cfg.add_section(section)
    cfg.set(section, 'TIME_STEP', str(par.TIME_STEP))
    
    section = 'Beam'
    cfg.add_section(section)
    cfg.set(section, 'TOTAL_TIME', str(par.TOTAL_TIME))
    
    section = 'Physics'
    cfg.add_section(section)
    cfg.set(section, 'ETCH_RATE', str(par.ETCH_RATE))
    
    section = 'Output'
    cfg.add_section(section)
    cfg.set(section, 'PLOT_SURFACE', str(par.PLOT_SURFACE))

    with open(output_file, 'w') as output:
        cfg.write(output)



try:
    check_flag = par.set_parameters(sys.argv[1])
except FileNotFoundError as err:
    print(err)
    sys.exit(-1)
except IndexError:
    print('No file specified as systemargument')
    sys.exit(-1)

_write_output(str(sys.argv[1]).replace('.cfg','.dat'))

        


