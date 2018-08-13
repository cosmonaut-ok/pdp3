import os
from xml.dom import minidom
from numpy import *
import json
import math

class Parameters:
    def __init__(self, parameters_file): # , movie_file=None): # , clim_e_field_r=[0,1], clim_e_field_z=[0,1]):
        '''
        read parameters from properties.xml, recalculate and set parameters as object fields
        '''

        self.config_path = os.path.dirname(parameters_file)
        # self.movie_file = movie_file if movie_file else os.path.join(self.config_path, 'field_movie.avi')

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

        if is_debug:
            self.data_dump_interval = int(file_save_parameters.getElementsByTagName('debug_data_dump_interval')[0].firstChild.data)
        else:
            self.data_dump_interval = int(file_save_parameters.getElementsByTagName('data_dump_interval')[0].firstChild.data)

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

        # get particles
        particles = self.dom_root.getElementsByTagName('particles')[0];
        self.particles = particles.getElementsByTagName('particle_specie')

        # get time
        time = self.dom_root.getElementsByTagName('time')[0];
        self.start_time = float(time.getElementsByTagName('start_time')[0].firstChild.data)
        self.end_time = float(time.getElementsByTagName('end_time')[0].firstChild.data)
        self.step_interval = float(time.getElementsByTagName('step_interval')[0].firstChild.data)

        # get plot parameters
        plot = self.dom_root.getElementsByTagName('plot')[0]

        # get video parameters
        video = plot.getElementsByTagName('video')[0]
        self.video_codec = str(plot.getElementsByTagName('codec')[0].firstChild.data)
        self.video_fps = int(plot.getElementsByTagName('fps')[0].firstChild.data)
        self.video_dpi = int(plot.getElementsByTagName('dpi')[0].firstChild.data)
        self.video_bitrate = int(plot.getElementsByTagName('bitrate')[0].firstChild.data)
        # get clim estimation
        self.get_clim_estimation()

    def get_file_frame_number_by_timestamp(self, timestamp):
        number_frames = timestamp / self.step_interval / self.data_dump_interval
        fpf = self.frames_per_file
        file_number = int(number_frames // fpf)
        frame_number = int(number_frames % fpf)
        return [file_number, frame_number]

    def get_frame_from_data(self, image, frame_number):
        fpf = self.frames_per_file
        sr = self.r_grid_count
        sz = self.z_grid_count

        if frame_number > fpf:
            raise ValueError('Can not process frame %i. Out of range (%i)' % (frame_number, fpf))
        else:
            rstart = sr*sz*frame_number
            rend = sr*sz*(frame_number+1)
            frame = image[rstart:rend]

        return(frame)

    def get_frame_row_from_data(self, image, frame_number, row_number):
        fpf = self.frames_per_file
        sr = self.r_grid_count
        sz = self.z_grid_count

        if frame_number > fpf:
            raise ValueError('Can not process frame %i. Out of range (%i)' % (frame_number, fpf))
        elif  row_number > sr:
            raise ValueError('Can not process row $i in frame %i. Out of range (%i)' % (row_number, frame_number, sr))
        else:

            rstart = sr * sz * frame_number + row_number * sz
            rend = rstart + sz
            row = image[rstart:rend]

        return(row)

    def get_frame_point_from_data(self, image, frame_number, x, y):
        sr = self.r_grid_count
        sz = self.z_grid_count
        if  y > sz:
            raise ValueError('Can not process point x=%i y=%i in frame %i. Out of range (x=%i y=%i)' % (x, y, frame_number, sz, sr))
        row = self.get_frame_row_from_data(image, frame_number, y)

        return(row[x])


    def get_frame_col_from_data(self, image, frame_number, col_number):
        fpf = self.frames_per_file
        sr = self.r_grid_count
        sz = self.z_grid_count
        col = []

        if  col_number > sz:
            raise ValueError('Can not process column $i in frame %i. Out of range (%i)' % (col_number, frame_number, sr))
        else:
            frame = get_frame_from_data(image, frame_number)
            for i in range(0, sr):
                col.append(self.get_frame_point_from_data(frame, 0, i, col_number))

        return(col)


    def get_clim_estimation(self):
        p_factor_array = []
        factor = 1e-5 # very empirical coefitient :)
        counter = 0

        for i in self.particles:
            n_l = float(i.getElementsByTagName('left_density')[0].firstChild.data)
            n_r = float(i.getElementsByTagName('right_density')[0].firstChild.data)
            t = float(i.getElementsByTagName('temperature')[0].firstChild.data)

            n = n_l if n_l > n_r else n_r

            p_factor_array.append(factor * n * t * math.sqrt(self.step_interval))
            counter = counter+1

        self.clim_estimation = sum(p_factor_array)
