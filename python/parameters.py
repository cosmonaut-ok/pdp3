#!/usr/bin/env python

import os
from xml.dom import minidom
from numpy import *

class Parameters:
    def __init__(self, parameters_file, movie_file=None, clim_e_field_r=[0,1], clim_e_field_z=[0,1]):
        '''
        read parameters from properties.xml, recalculate and set parameters as object fields
        '''

        self.config_path = os.path.dirname(parameters_file)
        self.movie_file = movie_file if movie_file else os.path.join(self.config_path, 'field_movie.avi')

        ## read parameters_file
        dom_root = minidom.parse(parameters_file)

        # get geometry parameters for gird and ticks
        geometry = dom_root.getElementsByTagName('geometry')[0]
        self.r_grid_count = int(geometry.getElementsByTagName('n_grid_r')[0].firstChild.data)-1
        self.z_grid_count = int(geometry.getElementsByTagName('n_grid_z')[0].firstChild.data)-1
        self.r_size = float(geometry.getElementsByTagName('r_size')[0].firstChild.data)
        self.z_size = float(geometry.getElementsByTagName('z_size')[0].firstChild.data)

        # get file to save parameters
        file_save_parameters = dom_root.getElementsByTagName('file_save_parameters')[0]

        # get data_path
        local_data_path = file_save_parameters.getElementsByTagName('path_to_result')[0].firstChild.data

        # calculate data path
        if str.startswith(local_data_path, '/'):
            self.data_path = local_data_path
        else:
            self.data_path = os.path.join(self.config_path, local_data_path)

        # get system state path
        local_system_state_path = file_save_parameters.getElementsByTagName('path_to_save_state')[0].firstChild.data

        # getData frame number per data file
        self.frames_per_file = int(file_save_parameters.getElementsByTagName('frames_per_file')[0].firstChild.data)

        # calculate data path
        if str.startswith(local_system_state_path, '/'):
            self.system_state_path = local_system_state_path
        else:
            self.system_state_path = os.path.join(self.config_path, local_system_state_path)

        ## get normalization parameters
        bunch = dom_root.getElementsByTagName('particles_bunch');
        bunch_density = float(bunch.item(0).getElementsByTagName('density')[0].firstChild.data);

        self.bunch_density = bunch_density # particles density of electron/ion bunch

        self.clim_e_field_r = clim_e_field_r # r-component of E field
        self.clim_e_field_z = clim_e_field_z # z-component of E field
        self.clim_e_field_bunch = [-(bunch_density*1.6e-19), 0] # E field of electron/ion bunch
