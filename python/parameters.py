#!/usr/bin/env python

import os
# import argparse





from xml.dom import minidom
# from numpy import *
# from pylab import *

# from matplotlib import *
# import matplotlib.pyplot as plt
# import matplotlib.animation as ani

class Parameters:
    def __init__(self, parameters_file, clim_e_field_r=[0,1], clim_e_field_z=[0,1], movie_file=None):
        '''
        + r_size
        + z_size
        + r_grid_count
        + z_grid_count
        ===
        start_time
        end_time
        step_time # delta_t
        ===
        + config_path
        + data_path
        + system_state_path
        + movie_file
        ===
        + bunch_density
        ===
        + clim_e_field_r     # r-component of E field
        + clim_e_field_z     # z-component of E field
        + clim_e_field_bunch # E field of electron/ion bunch
        ===
        frames_per_file
        '''

        # self.clim_e_field_r = clim_e_field_r
        # self.clim_e_field_z = clim_e_field_z
        # self.clim_rho_beam = clim_rho_beam

        self.config_path = os.path.dirname(parameters_file)
        self.movie_file = movie_file if movie_file else os.path.join(self.config_path, 'field_movie_python.avi')

        ## read parameters_file
        dom_root = minidom.parse(parameters_file)

        # get geometry parameters for gird and ticks
        geometry = dom_root.getElementsByTagName('geometry')[0]
        self.r_grid_count = int(geometry.getElementsByTagName('n_grid_r')[0].firstChild.data)-1
        self.z_grid_count = int(geometry.getElementsByTagName('n_grid_z')[0].firstChild.data)-1
        self.r_size= float(geometry.getElementsByTagName('r_size')[0].firstChild.data)
        self.z_size = float(geometry.getElementsByTagName('z_size')[0].firstChild.data)

        # set number of marks with names in X and Y axes
        # x_ticks_number = 10;
        # y_ticks_number = 4;
        # %% Titles and Ticks
        # self.__x_axe_title = 'Z(m)';
        # self.__y_axe_title = 'R(m)';
        # self.__x_tick_range = linspace(0, x_tick_max, x_ticks_number) # we need 10 (or x_ticks_number) ticks
        # self.__x_tick_gird_size = linspace(0, self.__size_field_z, x_ticks_number) # from 0 to x_tick_max. it's required
        # self.__y_tick_range = linspace(0, y_tick_max, y_ticks_number) # to convert gird to real size (meters)
        # self.__y_tick_gird_size = linspace(0, self.__size_field_r, y_ticks_number) # Same for X and Y axes

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
        bunch = dom_root.getElementsByTagName('Particles_bunch');
        bunch_density = float(bunch.item(0).getElementsByTagName('density')[0].firstChild.data);

        self.bunch_density = bunch_density # particles density of electron/ion bunch

        self.clim_e_field_r = clim_e_field_r # r-component of E field
        self.clim_e_field_z = clim_e_field_z # z-component of E field
        self.clim_e_field_bunch = [-(bunch_density*1.6e-19), 0] # E field of electron/ion bunch


        # clim_rho_beam = [-(bunch_density*1.6e-19) 0];


        # self.__data_file_e_field_r = os.path.join(self.__data_path, 'e_field_r')
        # self.__data_file_e_field_z = os.path.join(self.__data_path, 'e_field_z')
        # self.__data_file_rho_beam = os.path.join(self.__data_path, 'rho_beam')
        # self.__figure = plt.figure()
