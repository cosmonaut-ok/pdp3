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

class Pdp2EZTViewBuilder:
    def __init__(self, cfg):
        ## define private object fields
        self._cfg = cfg
        self._plot_builder = PDP3PlotBuilder(self._cfg)

        fpf = self._cfg.frames_per_file

        ## define files data file sets range
        self.start_data_set = 0
        self.start_frame = 0
        self.end_data_set, self.end_frame = self._cfg.get_file_frame_number_by_timestamp(self._cfg.end_time)

        self.clim_e_field_beam_scale_factor = 1
        self.clim_e_field_r = [0, 1]
        self.clim_e_field_z = [0, 1]

        ## define public object fields
        self.data_file_e_r_pattern = 'E_r'
        self.data_file_e_z_pattern = 'E_z'
        self.data_file_e_bunch_density_pattern = 'rho_beam'

        self.x_axis_label = r'$\mathit{Z (m)}$'
        self.y_axis_label = r'$\mathit{t (ns)}$'
        self.cbar_axis_label = r'$\frac{V}{m}$'
        self.cbar_bunch_density_axis_label = r'$m^{-3}$'

        self.use_grid = False

        ## define public object fields
        self.cmap = 'gray'

        self.E_r_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Radial\enspace Component}\enspace(E_r)$'
        self.E_z_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Longitudal\enspace Component}\enspace(E_z)$'

        self.image_t_range = (self.end_data_set - self.start_data_set) * fpf
        self.image_z_range = cfg.z_grid_count

        self.start_time = 0
        self.end_time = 0

        self.radius = 0

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

        self._plot_builder.x_tick_count = 4
        self._plot_builder.y_tick_count = 10
        self._plot_builder.y_tick_start = self.__start_time * 1e9 # nano- is for 1e-9
        self._plot_builder.y_tick_end = self.__end_time * 1e9 # nano- is for 1e-9

        font_size=14

        ## setup plot dimensions
        self.image_t_range = end_number_frames - start_number_frames

        self._plot_builder.x_plot_size = self.image_z_range
        self._plot_builder.y_plot_size = self.image_t_range

        self._plot_builder.setup_figure(width=40, height=30, font_size=font_size, dpi=100)

        ## setup E_r plot
        self._plot_builder.add_subplot_with_image(
            self.E_r_plot_name, 121, cmap=self.cmap, clim=self.clim_e_field_r
        )
        self._plot_builder.setup_subplot(
            self.E_r_plot_name, x_axe_label=self.x_axis_label,
            y_axe_label=self.y_axis_label, grid=self.use_grid
        )
        self._plot_builder.add_colorbar(
            self.E_r_plot_name, ticks=self.clim_e_field_r, title=self.cbar_axis_label, font_size=font_size-4, size="10%"
        )

        ## setup E_z plot
        self._plot_builder.add_subplot_with_image(
            self.E_z_plot_name, 122, cmap=self.cmap, clim=self.clim_e_field_z
        )
        self._plot_builder.setup_subplot(
            self.E_z_plot_name, x_axe_label=self.x_axis_label,
            y_axe_label=self.y_axis_label, grid=self.use_grid
        )
        self._plot_builder.add_colorbar(
            self.E_z_plot_name, ticks=self.clim_e_field_z, title=self.cbar_axis_label, font_size=font_size-4, size="10%"
        )

    def create_view_with_2_plots(self):
        '''
        create movie with preset subplots and data from data files
        '''
        fpf = self._cfg.frames_per_file
        sr = self._cfg.r_grid_count
        sz = self._cfg.z_grid_count
        radius_row = self._cfg.get_row_by_radius(self.radius)

        tiny_cache = TinyCache(os.path.join(self._cfg.data_path, '.cache'))
        cache_r_name = format('image_erz_r_%e-%e-%e#%d_%d' % (
            self.__start_time,
            self.__end_time,
            self._cfg.step_interval,
            self._cfg.data_dump_interval,
            radius_row
        ))
        cache_z_name = format('image_erz_z_%e-%e-%e#%d_%d' % (
            self.__start_time,
            self.__end_time,
            self._cfg.step_interval,
            self._cfg.data_dump_interval,
            radius_row
        ))

        image_r = tiny_cache.get_cache(cache_r_name)
        image_z = tiny_cache.get_cache(cache_z_name)

        if (len(image_r) == 0 or len(image_z) == 0):
            data_file_e_r = os.path.join(self._cfg.data_path, self.data_file_e_r_pattern)
            data_file_e_z = os.path.join(self._cfg.data_path, self.data_file_e_z_pattern)
            data_file_bunch_density = os.path.join(self._cfg.data_path, self.data_file_e_bunch_density_pattern)

            # with writer.saving(self._plot_builder.figure, movie_file, self.video_dpi):
            for k in range(self.start_data_set, self.end_data_set+1):
                tstart = (k*fpf+self.start_frame) if k == self.start_data_set else k*fpf
                tend = (k*fpf+self.end_frame) if k == self.end_data_set else ((k+1)*fpf)
                i = 1;

                if not os.path.isfile(data_file_e_r + str(k)) \
                   or not os.path.isfile(data_file_e_z + str(k)) \
                   or not os.path.isfile(data_file_bunch_density + str(k)):
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

                    image_r.extend(self._cfg.get_frame_row_from_data(h_field_e_r, local_step, radius_row))
                    image_z.extend(self._cfg.get_frame_row_from_data(h_field_e_z, local_step, radius_row))
                tiny_cache.update_cache(cache_r_name, image_r)
                tiny_cache.update_cache(cache_z_name, image_z)
        self._plot_builder.fill_image_with_data(
            self.E_z_plot_name, image_r)

        self._plot_builder.fill_image_with_data(
            self.E_r_plot_name, image_z)

        # self._plot_builder.redraw()
