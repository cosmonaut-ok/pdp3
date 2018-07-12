import os
from xml.dom import minidom
from numpy import *
import json

class Parameters:
    def __init__(self, parameters_file, movie_file=None, clim_e_field_r=[0,1], clim_e_field_z=[0,1]):
        '''
        read parameters from properties.xml, recalculate and set parameters as object fields
        '''

        self.config_path = os.path.dirname(parameters_file)
        self.movie_file = movie_file if movie_file else os.path.join(self.config_path, 'field_movie.avi')

        ## read parameters_file
        self.dom_root = minidom.parse(parameters_file)

        # check, if debug used
        is_debug = json.loads(self.dom_root.getElementsByTagName('debug')[0].firstChild.data.lower())

        if is_debug:
            n_grid_r_name = 'debug_n_grid_r'
            n_grid_z_name = 'debug_n_grid_z'
        else:
            n_grid_r_name = 'n_grid_r'
            n_grid_z_name = 'n_grid_z'

        # get geometry parameters for gird and ticks
        geometry = self.dom_root.getElementsByTagName('geometry')[0]
        self.r_grid_count = int(geometry.getElementsByTagName(n_grid_r_name)[0].firstChild.data)-1
        self.z_grid_count = int(geometry.getElementsByTagName(n_grid_z_name)[0].firstChild.data)-1
        self.r_size = float(geometry.getElementsByTagName('r_size')[0].firstChild.data)
        self.z_size = float(geometry.getElementsByTagName('z_size')[0].firstChild.data)

        # get file to save parameters
        file_save_parameters = self.dom_root.getElementsByTagName('file_save_parameters')[0]

        # get data_path
        local_data_path = file_save_parameters.getElementsByTagName('result_path')[0].firstChild.data

        # calculate data path
        if str.startswith(local_data_path, '/'):
            self.data_path = local_data_path
        else:
            self.data_path = os.path.join(self.config_path, local_data_path)

        # get system state path
        local_system_state_path = file_save_parameters.getElementsByTagName('save_state_path')[0].firstChild.data

        # getData frame number per data file
        self.frames_per_file = int(file_save_parameters.getElementsByTagName('frames_per_file')[0].firstChild.data)

        # calculate data path
        if str.startswith(local_system_state_path, '/'):
            self.system_state_path = local_system_state_path
        else:
            self.system_state_path = os.path.join(self.config_path, local_system_state_path)

        ## get normalization parameters
        beam = self.dom_root.getElementsByTagName('particle_beam');
        self.bunch_density = float(beam.item(0).getElementsByTagName('bunch_density')[0].firstChild.data)
        self.bunch_radius = float(beam.item(0).getElementsByTagName('bunch_radius')[0].firstChild.data)
        self.bunch_length = float(beam.item(0).getElementsByTagName('bunch_length')[0].firstChild.data)
        self.number_bunches = float(beam.item(0).getElementsByTagName('number_bunches')[0].firstChild.data)
        self.bunches_distance = float(beam.item(0).getElementsByTagName('bunches_distance')[0].firstChild.data)
        self.beam_initial_velocity = float(beam.item(0).getElementsByTagName('initial_velocity')[0].firstChild.data)

        self.clim_e_field_r = clim_e_field_r # r-component of E field
        self.clim_e_field_z = clim_e_field_z # z-component of E field
        self.clim_e_field_beam = [-(self.bunch_density*1.6e-19), 0] # E field of electron/ion beam

        # get particles
        particles = self.dom_root.getElementsByTagName('particles')[0];
        self.particles = particles.getElementsByTagName('particle_specie')

        # get time
        time = self.dom_root.getElementsByTagName('time')[0];
        self.start_time = float(time.getElementsByTagName('start_time')[0].firstChild.data)
        self.end_time = float(time.getElementsByTagName('end_time')[0].firstChild.data)
        self.step_interval = float(time.getElementsByTagName('step_interval')[0].firstChild.data)
