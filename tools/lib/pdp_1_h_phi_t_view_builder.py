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

class Pdp1HTViewBuilder:
    def __init__(self, cfg):
        ## define private object fields
        self._cfg = cfg
        self._plot_builder = PDP3PlotBuilder(self._cfg)

        fpf = self._cfg.frames_per_file

        ## define files data file sets range
        self.start_data_set = 0
        self.start_frame = 0
        self.end_data_set, self.end_frame = self._cfg.get_file_frame_number_by_timestamp(self._cfg.end_time)

        self.range_h_field_phi = [0,0]

        ## define public object fields
        self.data_file_h_phi_pattern = 'H_phi'

        self.x_axis_label = r'$\mathit{t (s)}$'
        self.y_axis_label = r'$\mathit{H_{\phi} (\frac{V}{m})}$'

        self.use_grid = False

        self.H_phi_plot_name = r'$\mathbf{Magnetic\enspace Field\enspace Rotation\enspace Component}\enspace(H_{\phi})$'

        self.__image_t_range = (self.end_data_set - self.start_data_set) * fpf

        self.radius = 0
        self.longitude = 0

    def setup_1h_view(self):
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

        ## setup h_phi plot
        subplot = self._plot_builder.add_subplot(self.H_phi_plot_name, 111)

        subplot.set_title(self.H_phi_plot_name)
        subplot.set_aspect('auto')
        subplot.set_xlabel(self.x_axis_label)
        subplot.set_ylabel(self.y_axis_label, rotation=45)
        subplot.spines['top'].set_visible(False)
        subplot.spines['right'].set_visible(False)
        subplot.ticklabel_format(style='sci', scilimits=(0,0))
        subplot.grid(self.use_grid)
        if self.range_h_field_phi[0] != 0 or self.range_h_field_phi[1] != 0: subplot.set_ylim(self.range_h_field_phi)


    def create_view_with_1_plot(self):
        '''
        create movie with preset subplots and data from data files
        '''
        fpf = self._cfg.frames_per_file
        sr = self._cfg.r_grid_count
        sz = self._cfg.z_grid_count
        radius_row = self._cfg.get_row_by_radius(self.radius)
        longitude_col = self._cfg.get_col_by_longitude(self.longitude)
        cache_name = format('plot2d_h_phi_%e-%e-%e#%d_%d_%d' % (
            self.__start_time,
            self.__end_time,
            self._cfg.step_interval,
            self._cfg.data_dump_interval,
            radius_row, longitude_col
        ))

        tiny_cache = TinyCache(os.path.join(self._cfg.data_path, '.cache'))

        plot2d_timeline = linspace(self.__start_time, self.__end_time, self.__image_t_range)
        plot2d = tiny_cache.get_cache(cache_name)

        if len(plot2d) == 0:
            data_file = os.path.join(self._cfg.data_path, self.data_file_h_phi_pattern)

            # with writer.saving(self._plot_builder.figure, movie_file, self.video_dpi):
            for k in range(self.start_data_set, self.end_data_set+1):
                tstart = (k*fpf+self.start_frame) if k == self.start_data_set else k*fpf
                tend = (k*fpf+self.end_frame) if k == self.end_data_set else ((k+1)*fpf)
                i = 1;

                if not os.path.isfile(data_file + str(k)):
                    print('No more data files exists. Exiting')
                    return

                print("Loading files set %d" % (k))
                ## Open data files
                fidh = open(data_file + str(k), 'r')

                h_field = fromfile(fidh, dtype=float, count=sr*sz*fpf, sep=' ')

                ## Close data files
                fidh.close()

                for t in range(tstart, tend):
                    local_step = t % fpf

                    print("Processing frame %d" % (local_step))

                    plot2d = append(plot2d, self._cfg.get_frame_point_from_data(h_field, local_step, longitude_col, radius_row))

        tiny_cache.update_cache(cache_name, plot2d)

        plot_h_phi = self._plot_builder.get_subplot(self.H_phi_plot_name)

        plot_h_phi.plot(plot2d_timeline, plot2d)
