# -*- coding: utf-8 -*-
"""
Read parameters from database file and config file

Module functions:
    set_parameters(config_file): Set the parameters as global variables

Created on Sat Nov 17 11:54:44 2018
@author: Christian Eder
"""

import sys
import os
import configparser


def set_parameters(config_file):
    """set global variables and check conditions
    
    :param config_file: name of the config file
    """
    db_file = os.path.join(os.path.dirname(__file__), 'parameters.db')
    
    #read values from file
    db_params = _read_file(db_file)
    db_config = _read_file(config_file)
    
    #set global variables
    for key, value in db_params.items():
        globals()[key] = value[0]
        
    #parameter and type checks
    check_params = _check_parameters(db_params, db_config)
    check_types = _check_types(db_params, db_config)
        
    if check_params and check_types:
        for key, value in db_config.items(): 
            if type(globals()[key]) is float:
                globals()[key] = float(value)
            else:
                globals()[key] = value  
        
        #condition check
        if not _check_conditions(db_params):
            sys.exit(-1)
    else:
        sys.exit(-1)   
    

def _read_file(file):
    """Read values from file with ConfigParser module and return dictionary
    
    :param file: name of the file with ConfigParser shape
    
    :return cfg_dict: dictionary with values from file
    """
    if os.path.isfile(file) is False:
        raise FileNotFoundError(file + ' Not Found')
        
    cfg = configparser.ConfigParser()
    cfg.read(file)
    
    cfg_dict = {}
    
    for section in cfg.sections():
        for option in cfg.options(section):
            cfg_dict[option.upper()] = eval(cfg[section][option])
            
    return cfg_dict
     
       
def _check_parameters(db_params, db_config):
    """Check if parameters with type in database defined in config file
    
    :param db_params: dictionary with values from database file
    :param db_config: dictionary with values from config file
    
    :return check_flag: return True if all checks are conformed, else False
    """
    check_flag = True
    
    for key, value in db_params.items():
        if type(value[0]) is type:
           if key not in db_config:
                check_flag = False
                print('Parameter not given ' + key)
        
    return check_flag


def _check_types(db_params, db_config):
    """Check if all types match in database and config file 
    
    Exception: int should be converted in float
    
    :param db_params: dictionary with values from database file
    :param db_config: dictionary with values from config file
    
    :return check_flag: return True if all checks are conformed, else False
    """
    check_flag = True
    
    for key, value in db_config.items():
        if key in globals():      
            if type(globals()[key]) is not type(value):
                if type(globals()[key]) is float and type(value) is not int:
                    print('Datatypes not match - ' + key)
                    check_flag = False
        else:
            print('Undefined parameter - ' + key)
            check_flag = False
    
    return check_flag
  
          
def _check_conditions(db_params):
    """Check if all conditions conformed
    
    :param db_params: dictionary with values from database file
    
    :return check_flag: return True if all checks are conformed, else False
    """
    check_flag = True
    
    for key, value in db_params.items():
        if value[1] is not None:
            if not eval(value[1]):
                print('Invalid Condition - ' + value[1])
                check_flag = False
    
    return check_flag


