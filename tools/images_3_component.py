#!/usr/bin/env python3

import sys
import os

import argparse

from numpy import *
from pylab import *

import matplotlib.pyplot as plt
import matplotlib.animation as ani

from parameters import Parameters
from pdp3_plot_builder import PDP3PlotBuilder

## README:
## color map reference: https://matplotlib.org/examples/color/colormaps_reference.html
## mathtext reference:  https://matplotlib.org/users/mathtext.html

class Pdp3Images:
    def __init__(self, cfg):
        ## define private object fields
        self.__cfg = cfg
        self.__plot_builder = PDP3PlotBuilder(self.__cfg)

        self.timestamp = 0

        ## define public object fields
        self.data_file_e_r_pattern = 'e_r'
        self.data_file_e_z_pattern = 'e_z'
        self.data_file_e_bunch_density_pattern = 'bunch_density'

        self.x_axis_label = r'$Z (m)$'
        self.y_axis_label = r'$R (m)$'
        self.cbar_axis_label = r'$V/m$'
        self.cbar_bunch_density_axis_label = r'm^{-3}$'

        self.position_e_r = [0.1, 0.70, 0.8, 0.3]
        self.position_e_z = [0.1, 0.35, 0.8, 0.3]
        self.position_bunch_density = [0.1, 0.01, 0.8, 0.3]

        ## define public object fields
        self.cmap = 'gray'

        self.E_r_plot_name = r'$E_r$'
        self.E_z_plot_name = r'$E_z$'
        self.E_bunch_density_plot_name = r'$\rho_{beam}$'

    # def get_file_frame(self):
    #     number_frames = self.timestamp / self.__cfg.step_interval / self.__cfg.data_dump_interval
    #     fpf = self.__cfg.frames_per_file
    #     file_number = int(number_frames // fpf)
    #     frame_number = int(number_frames % fpf)
    #     return [file_number, frame_number]


    def setup_plot(self, view):
        '''
        initialize plot figure and subplots with preset object fields
        '''
        self.__plot_builder.setup_figure()
        ## setup E_r plot
        self.__plot_builder.add_subplot_with_image(
            self.E_r_plot_name, 311, cmap=self.cmap, clim=self.__cfg.clim_e_field_r
        )
        self.__plot_builder.setup_subplot(
            self.E_r_plot_name, x_axe_label=self.x_axis_label,
            y_axe_label=self.y_axis_label, position=self.position_e_r
        )
        self.__plot_builder.add_colorbar(
            self.E_r_plot_name, ticks=self.__cfg.clim_e_field_r, title=self.cbar_axis_label
        )

        ## setup E_z plot
        self.__plot_builder.add_subplot_with_image(
            self.E_z_plot_name, 312, cmap=self.cmap, clim=self.__cfg.clim_e_field_z
        )
        self.__plot_builder.setup_subplot(
            self.E_z_plot_name, x_axe_label=self.x_axis_label,
            y_axe_label=self.y_axis_label, position=self.position_e_z
        )
        self.__plot_builder.add_colorbar(
            self.E_z_plot_name, ticks=self.__cfg.clim_e_field_z, title=self.cbar_axis_label
        )

        ## setup bunch_density plot
        self.__plot_builder.add_subplot_with_image(
            self.E_bunch_density_plot_name, 313, cmap=self.cmap, clim=self.__cfg.clim_e_field_beam
        )
        self.__plot_builder.setup_subplot(
            self.E_bunch_density_plot_name, x_axe_label=self.x_axis_label,
            y_axe_label=self.y_axis_label, position=self.position_bunch_density
        )
        self.__plot_builder.add_colorbar(
            self.E_bunch_density_plot_name, ticks=self.__cfg.clim_e_field_beam,
            ticklabels=[self.__cfg.bunch_density, 0], title=self.cbar_bunch_density_axis_label
        )

        if view:
            self.__plot_builder.figure.show()

    def create_images_with_3_plots(self, view=False, write=True):
        '''
        create images with preset subplots and data from data files
        '''
        # start_frame = self.start_frame
        fpf = self.__cfg.frames_per_file
        # end_frame = fpf if self.end_frame == -1 else self.end_frame
        sr = self.__cfg.r_grid_count
        sz = self.__cfg.z_grid_count

        data_file_e_r = os.path.join(self.__cfg.data_path, self.data_file_e_r_pattern)
        data_file_e_z = os.path.join(self.__cfg.data_path, self.data_file_e_z_pattern)
        data_file_bunch_density = os.path.join(self.__cfg.data_path, self.data_file_e_bunch_density_pattern)

        data_file, data_frame = self.__cfg.get_file_frame_by_timestamp(self.timestamp)

        k = data_file

        if not os.path.isfile(data_file_e_r + str(k)) \
           or not os.path.isfile(data_file_e_z + str(k)) \
           or not os.path.isfile(data_file_bunch_density + str(k)):
            print('No more data files exists. Exiting')
            return

        print("Loading files set %d" % (k))
        ## Open data files
        fidh_e_r = open(data_file_e_r + str(k), 'r')
        fidh_e_z = open(data_file_e_z + str(k), 'r')
        fidh_bunch_density = open(data_file_bunch_density + str(k), 'r')

        h_field_e_r = fromfile(fidh_e_r, dtype=float, count=sr*sz*fpf, sep=' ')
        h_field_e_z = fromfile(fidh_e_z, dtype=float, count=sr*sz*fpf, sep=' ')
        h_field_bunch_density = fromfile(fidh_bunch_density, dtype=float, count=sr*sz*fpf, sep=' ')

        ## Close data files
        fidh_e_r.close()
        fidh_e_z.close()
        fidh_bunch_density.close()

        # for t in range(tstart, tend):
        local_step = data_frame # t % fpf

        print("Processing frame %d" % (local_step))
        image_file_name = os.path.join(self.__cfg.config_path, 'image_' + str(k) + '_' + str(local_step) + '.png')

        rstart = sr*sz*local_step
        rend = sr*sz*(local_step+1)
        try:
            self.__plot_builder.fill_image_with_data(
                self.E_z_plot_name,
                h_field_e_r[rstart:rend])

            self.__plot_builder.fill_image_with_data(
                self.E_r_plot_name,
                h_field_e_z[rstart:rend])

            self.__plot_builder.fill_image_with_data(
                self.E_bunch_density_plot_name,
                h_field_bunch_density[rstart:rend])

        except ValueError: ## skip frame, when data is inconsistent
            return

        if write:
            self.__plot_builder.figure.savefig(image_file_name) # save the figure to file
        if view:
            self.__plot_builder.redraw()

def main():
    ####
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    parser.add_argument('--timestamp', type=float, help='Timestamp to generate image at', required=True)

    parser.add_argument('--dry-run', action='store_true', help='Do not write anything. For debug')

    parser.add_argument('--view', action='store_true', default=False,
                        help='View image as well as write')

    default_clim = [-1e5, 1e5]

    parser.add_argument('--clim-e-r', type=str,
                        help='Color limit range for Electrical field radial component. Default %s'
                        % ':'.join(map(str, default_clim)))
    parser.add_argument('--clim-e-z', type=str,
                        help='Color limit range for Electrical field longitudal component. Default %s'
                        % ':'.join(map(str, default_clim)))

    args = parser.parse_args()

    view=False
    write=True

    if args.view:
        view = True
    if args.dry_run:
        write = False

    clim_e_r = list(map(float, args.clim_e_r.split(':'))) if args.clim_e_r else default_clim
    clim_e_z = list(map(float, args.clim_e_z.split(':'))) if args.clim_e_z else default_clim

    # check if config file exists
    if os.path.isfile(args.properties_path):
        ## initialize config
        config = Parameters(args.properties_path, None, clim_e_r, clim_e_z)

        images = Pdp3Images(config)

        ################################################################################################
        #################### configure plot and images parameters #######################################
        ################################################################################################
        # images.start_frame = frame_range[0]
        # images.end_frame = frame_range[1]

        # images.start_data_set = data_set_range[0]
        # images.end_data_set = data_set_range[1]

        images.timestamp = args.timestamp

        images.data_file_e_r_pattern = 'E_r'
        images.data_file_e_z_pattern = 'E_z'
        images.data_file_e_bunch_density_pattern = 'rho_beam'

        images.x_axis_label = r'$\mathit{Z (m)}$'
        images.y_axis_label = r'$\mathit{R (m)}$'
        images.cbar_axis_label = r'$\frac{V}{m}$'
        images.cbar_bunch_density_axis_label = r'$m^{-3}$'

        images.position_e_r = [0.1, 0.70, 0.8, 0.3]
        images.position_e_z = [0.1, 0.35, 0.8, 0.3]
        images.position_bunch_density = [0.1, 0.01, 0.8, 0.3]

        images.cmap = 'gray'

        images.E_r_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Radial\enspace Component}\enspace(E_r)$'
        images.E_z_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Longitudal\enspace Component}\enspace(E_z)$'
        images.E_bunch_density_plot_name = r'$\mathbf{Electron\enspace Beam\enspace Density}\enspace (\rho_{beam})$'
        ################################################################################################
        ################################################################################################
        ################################################################################################

        images.setup_plot(view)
        images.create_images_with_3_plots(view, write)
        if view:
            input("Press 'Return' to exit ")
    else:
        print("Configuration file `%s' does not exists. Exiting" % args.properties_path)
        exit(1)

## call main function
if __name__ == "__main__":
    main()
