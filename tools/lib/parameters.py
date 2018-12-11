import os
from xml.dom import minidom
from numpy import *
import json
import math
from os.path import normpath


class Probe:
    def __init__(self, probe_type, component, schedule, r_start=-1, z_start=-1, r_end=-1, z_end=-1):
        if (probe_type == 'frame' or probe_type == 'col' or probe_type == 'row' or probe_type == 'dot' or probe_type == 'mpframe'):
            self.type = probe_type

        self.component = component
        self.schedule = schedule
        self.r_start = r_start
        self.z_start = z_start
        self.r_end = r_end
        self.z_end = z_end


class ParticlesSpecie:
    def __init__(self, name, number_macro,
                 left_density, right_density, temperature,
                 charge=1, mass=1):
        self.name = name
        self.charge = charge
        self.mass = mass
        self.number_macro = number_macro
        self.left_density = left_density
        self.right_density = right_density
        self.temperature = temperature


class Parameters:
    def __init__(self, parameters_file): # , movie_file=None): # , clim_e_field_r=[0,1], clim_e_field_z=[0,1]):
        '''
        read parameters from properties.xml, recalculate and set parameters as object fields
        '''

        self.config_path = normpath(os.path.dirname(parameters_file))
        # self.movie_file = movie_file if movie_file else os.path.join(self.config_path, 'field_movie.avi')

        ## read parameters_file
        self.dom_root = minidom.parse(parameters_file)

        # check, if debug used
        is_debug = json.loads(self.dom_root.getElementsByTagName('debug')[0].firstChild.data.lower())

        # use hdf5 format
        self.use_hdf5 = json.loads(self.dom_root.getElementsByTagName('use_hdf5')[0].firstChild.data.lower())

        if is_debug:
            n_grid_r_name = 'debug_n_grid_r'
            n_grid_z_name = 'debug_n_grid_z'
        else:
            n_grid_r_name = 'n_grid_r'
            n_grid_z_name = 'n_grid_z'

        # get geometry parameters for gird and ticks
        geometry = self.dom_root.getElementsByTagName('geometry')[0]
        self.number_r_grid = int(geometry.getElementsByTagName(n_grid_r_name)[0].firstChild.data)
        self.number_z_grid = int(geometry.getElementsByTagName(n_grid_z_name)[0].firstChild.data)
        self.r_size = float(geometry.getElementsByTagName('r_size')[0].firstChild.data)
        self.z_size = float(geometry.getElementsByTagName('z_size')[0].firstChild.data)

        # get file to save parameters
        file_save_parameters = self.dom_root.getElementsByTagName('file_save_parameters')[0]

        probes = file_save_parameters.getElementsByTagName('probes')[0].getElementsByTagName('probe')
        self.probes = []

        for i in probes:
            p_type = i.getAttribute('type')
            p_compon = i.getAttribute('component')
            p_schedule = int(i.getAttribute('schedule'))

            r_start = -1
            r_end = -1
            z_start = -1
            z_end = -1

            if p_type == 'frame' or p_type == 'mpframe':
                r_start = int(i.getAttribute('r_start'))
                r_end = int(i.getAttribute('r_end'))
                z_start = int(i.getAttribute('z_start'))
                z_end = int(i.getAttribute('z_end'))
                p = Probe(p_type, p_compon, p_schedule,
                          r_start, z_start,
                          r_end, z_end)

            elif p_type == 'col':
                z_start = int(i.getAttribute('z'))
                p = Probe(p_type, p_compon, p_schedule,
                          r_start, z_start,
                          r_end, z_end)

            elif p_type == 'row':
                r_start = int(i.getAttribute('r'))
                p = Probe(p_type, p_compon, p_schedule,
                          r_start, z_start,
                          r_end, z_end)

            elif p_type == 'dot':
                r_start = int(i.getAttribute('r'))
                z_start = int(i.getAttribute('z'))
                p = Probe(p_type, p_compon, p_schedule,
                          r_start, z_start,
                          r_end, z_end)
            else:
                raise NameError('probe type can not be {}'.format(p_type))

            self.probes.append(p)

        # set data root path
        if is_debug:
            self.data_root = file_save_parameters.getElementsByTagName('debug_data_root')[0].firstChild.data
        else:
            self.data_root = file_save_parameters.getElementsByTagName('data_root')[0].firstChild.data

        # get data_path
        local_data_path = file_save_parameters.getElementsByTagName('result_path')[0].firstChild.data

        # calculate data path
        if str.startswith(local_data_path, '/'):
            self.data_path = local_data_path
        else:
            self.data_path = normpath(os.path.join(self.config_path, local_data_path))

        # getData frame number per data file
        self.frames_per_file = int(file_save_parameters.getElementsByTagName('frames_per_file')[0].firstChild.data)

        ## get normalization parameters
        beam = self.dom_root.getElementsByTagName('particle_beam');
        self.bunch_density = float(beam.item(0).getElementsByTagName('bunch_density')[0].firstChild.data)
        self.bunch_radius = float(beam.item(0).getElementsByTagName('bunch_radius')[0].firstChild.data)
        self.bunch_length = float(beam.item(0).getElementsByTagName('bunch_length')[0].firstChild.data)
        self.number_bunches = float(beam.item(0).getElementsByTagName('number_bunches')[0].firstChild.data)
        self.bunches_distance = float(beam.item(0).getElementsByTagName('bunches_distance')[0].firstChild.data)
        self.beam_initial_velocity = float(beam.item(0).getElementsByTagName('initial_velocity')[0].firstChild.data)

        # get particles
        particles = self.dom_root.getElementsByTagName('particles')[0].getElementsByTagName('particle_specie')
        self.particles = [] #  particles.getElementsByTagName('particle_specie')

        for i in particles:
            p_name = i.getAttribute('name')
            p_charge = int(i.getElementsByTagName('charge')[0].firstChild.data)
            p_mass = int(i.getElementsByTagName('mass')[0].firstChild.data)
            p_temperature = float(i.getElementsByTagName('temperature')[0].firstChild.data)
            p_number_macro = float(i.getElementsByTagName('number_macro')[0].firstChild.data)
            p_left_density = float(i.getElementsByTagName('left_density')[0].firstChild.data)
            p_right_density = float(i.getElementsByTagName('right_density')[0].firstChild.data)

            p = ParticlesSpecie(p_name, p_number_macro,
                                p_left_density, p_right_density, p_temperature,
                                p_charge, p_mass)

            self.particles.append(p)


        # get time
        time = self.dom_root.getElementsByTagName('time')[0];
        self.start_time = float(time.getElementsByTagName('start_time')[0].firstChild.data)
        self.end_time = float(time.getElementsByTagName('end_time')[0].firstChild.data)
        self.step_interval = float(time.getElementsByTagName('step_interval')[0].firstChild.data)

        # get plot parameters
        plot = self.dom_root.getElementsByTagName('plot')[0]

        # get plot figure parameters
        figure = plot.getElementsByTagName('figure')[0]
        # get video parameters
        self.figure_color = str(figure.getElementsByTagName('color')[0].firstChild.data)
        self.figure_width = float(figure.getElementsByTagName('width')[0].firstChild.data)
        self.figure_height = float(figure.getElementsByTagName('height')[0].firstChild.data)
        self.figure_dpi = int(figure.getElementsByTagName('dpi')[0].firstChild.data)

        font = figure.getElementsByTagName('font')[0]
        self.figure_font_name = str(font.getElementsByTagName('name')[0].firstChild.data)
        self.figure_font_family = str(font.getElementsByTagName('family')[0].firstChild.data)
        self.figure_font_size = int(font.getElementsByTagName('size')[0].firstChild.data)

        video = plot.getElementsByTagName('video')[0]
        self.video_codec = str(video.getElementsByTagName('codec')[0].firstChild.data)
        self.video_fps = int(video.getElementsByTagName('fps')[0].firstChild.data)
        self.video_bitrate = int(video.getElementsByTagName('bitrate')[0].firstChild.data)


    def get_frame_number_by_timestamp(self, timestamp, dump_interval):
        number_frames = round(timestamp / self.step_interval / dump_interval)
        return number_frames


    def get_timestamp_by_frame_number(self, frame_number, dump_interval):
        timestamp = frame_number * self.step_interval * dump_interval
        return timestamp


    def get_clim_estimation(self):
        p_factor_array = []
        factor = 1e-5 # very empirical coefitient :)
        counter = 0

        for i in self.particles:
            n_l = i.left_density
            n_r = i.right_density
            t = i.temperature

            n = n_l if n_l > n_r else n_r

            p_factor_array.append(factor * n * t * math.sqrt(self.step_interval))
            counter = counter+1

        return sum(p_factor_array)


    def get_row_by_radius(self, radius):
        return(int(round(radius * self.number_r_grid / self.r_size)))


    def get_col_by_longitude(self, longitude):
        return(int(round(longitude * self.number_z_grid / self.z_size)))
