import sys
import os

from numpy import *
from pylab import *

import matplotlib.pyplot as plt
import matplotlib.animation as ani

from lib.parameters import Parameters
from lib.tinycache import TinyCache
from lib.pdp3_plot_builder import PDP3PlotBuilder

## README:
## color map reference: https://matplotlib.org/examples/color/colormaps_reference.html
## mathtext reference:  https://matplotlib.org/users/mathtext.html

class Pdp2ETViewBuilder:
    def __init__(self, cfg):
        ## define private object fields
        self._cfg = cfg
        self._plot_builder = PDP3PlotBuilder(self._cfg)

        fpf = self._cfg.frames_per_file

        ## define files data file sets range
        self.start_data_set = 0
        self.start_frame = 0
        self.end_data_set, self.end_frame = self._cfg.get_file_frame_number_by_timestamp(self._cfg.end_time)

        self.range_e_field_r = [0,0]
        self.range_e_field_z = [0,0]

        ## define public object fields
        self.data_file_e_r_pattern = 'E_r'
        self.data_file_e_z_pattern = 'E_z'

        self.x_axis_label = r'$\mathit{t (s)}$'
        self.y_r_axis_label = r'$\mathit{E_r (\frac{V}{m})}$'
        self.y_z_axis_label = r'$\mathit{E_z (\frac{V}{m})}$'

        self.use_grid = False

        self.E_r_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Radial\enspace Component}\enspace(E_r)$'
        self.E_z_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Longitudal\enspace Component}\enspace(E_z)$'

        self.__image_t_range = (self.end_data_set - self.start_data_set) * fpf

        self.radius = 0
        self.longitude = 0

    def setup_2e_view(self):
        '''
        initialize plot figure and subplots with preset object fields
        '''
        fpf = self._cfg.frames_per_file

        ## set number of ticks
        start_number_frames = self.start_data_set * fpf + self.start_frame
        end_number_frames = self.end_data_set * fpf + self.end_frame
        self.__start_time = start_number_frames * self._cfg.data_dump_interval * self._cfg.step_interval
        self.__end_time = end_number_frames * self._cfg.data_dump_interval * self._cfg.step_interval
        self.__image_t_range = end_number_frames - start_number_frames

        font_size=14

        self._plot_builder.setup_figure(width=19.2, height=10.8, font_size=font_size, dpi=100)

        ## setup E_r plot
        subplot_er = self._plot_builder.add_subplot(self.E_r_plot_name, 121)
        subplot_ez = self._plot_builder.add_subplot(self.E_z_plot_name, 122)

        subplot_er.set_title(self.E_r_plot_name)
        subplot_er.set_aspect('auto')
        subplot_er.set_xlabel(self.x_axis_label)
        subplot_er.set_ylabel(self.y_r_axis_label, rotation=45)
        subplot_er.spines['top'].set_visible(False)
        subplot_er.spines['right'].set_visible(False)
        subplot_er.ticklabel_format(style='sci', scilimits=(0,0))
        subplot_er.grid(self.use_grid)
        if self.range_e_field_r[0] != 0 or self.range_e_field_r[1] != 0: subplot_er.set_ylim(self.range_e_field_r)

        subplot_ez.set_title(self.E_z_plot_name)
        subplot_ez.set_aspect('auto')
        subplot_ez.set_xlabel(self.x_axis_label)
        subplot_ez.set_ylabel(self.y_z_axis_label, rotation=45)
        subplot_ez.spines['top'].set_visible(False)
        subplot_ez.spines['right'].set_visible(False)
        subplot_ez.ticklabel_format(style='sci', scilimits=(0,0))
        subplot_ez.grid(self.use_grid)
        if self.range_e_field_z[0] != 0 or self.range_e_field_z[1] != 0: subplot_ez.set_ylim(self.range_e_field_z)


    def create_view_with_2_plots(self):
        '''
        create movie with preset subplots and data from data files
        '''
        fpf = self._cfg.frames_per_file
        sr = self._cfg.r_grid_count
        sz = self._cfg.z_grid_count
        radius_row = self._cfg.get_row_by_radius(self.radius)
        longitude_col = self._cfg.get_col_by_longitude(self.longitude)
        cache_r_name = format('plot2d_r_%e-%e-%e#%d_%d_%d' % (
            self.__start_time,
            self.__end_time,
            self._cfg.step_interval,
            self._cfg.data_dump_interval,
            radius_row, longitude_col
        ))
        cache_z_name = format('plot2d_z_%e-%e-%e#%d_%d_%d' % (
            self.__start_time,
            self.__end_time,
            self._cfg.step_interval,
            self._cfg.data_dump_interval,
            radius_row, longitude_col
        ))

        tiny_cache = TinyCache(os.path.join(self._cfg.data_path, '.cache'))

        plot2d_timeline = linspace(self.__start_time, self.__end_time, self.__image_t_range)
        plot2d_r = tiny_cache.get_cache(cache_r_name)
        plot2d_z = tiny_cache.get_cache(cache_z_name)

        if (len(plot2d_r) == 0 or len(plot2d_z) == 0 or len(plot2d_timeline) == 0):
            data_file_e_r = os.path.join(self._cfg.data_path, self.data_file_e_r_pattern)
            data_file_e_z = os.path.join(self._cfg.data_path, self.data_file_e_z_pattern)

            # with writer.saving(self._plot_builder.figure, movie_file, self.video_dpi):
            for k in range(self.start_data_set, self.end_data_set+1):
                tstart = (k*fpf+self.start_frame) if k == self.start_data_set else k*fpf
                tend = (k*fpf+self.end_frame) if k == self.end_data_set else ((k+1)*fpf)
                i = 1;

                if not os.path.isfile(data_file_e_r + str(k)) \
                   or not os.path.isfile(data_file_e_z + str(k)):
                    print('No more data files exists. Exiting')
                    return

                print("Loading files set %d" % (k))
                ## Open data files
                fidh_e_r = open(data_file_e_r + str(k), 'r')
                fidh_e_z = open(data_file_e_z + str(k), 'r')

                h_field_e_r = fromfile(fidh_e_r, dtype=float, count=sr*sz*fpf, sep=' ')
                h_field_e_z = fromfile(fidh_e_z, dtype=float, count=sr*sz*fpf, sep=' ')

                ## Close data files
                fidh_e_r.close()
                fidh_e_z.close()

                for t in range(tstart, tend):
                    local_step = t % fpf

                    print("Processing frame %d" % (local_step))

                    plot2d_r = append(plot2d_r, self._cfg.get_frame_point_from_data(h_field_e_r, local_step, longitude_col, radius_row))
                    plot2d_z = append(plot2d_z, self._cfg.get_frame_point_from_data(h_field_e_z, local_step, longitude_col, radius_row))

        tiny_cache.update_cache(cache_r_name, plot2d_r)
        tiny_cache.update_cache(cache_z_name, plot2d_z)

        plot_er = self._plot_builder.get_subplot(self.E_r_plot_name)
        plot_ez = self._plot_builder.get_subplot(self.E_z_plot_name)

        plot_er.plot(plot2d_timeline, plot2d_r)
        plot_ez.plot(plot2d_timeline, plot2d_z)
