#!/usr/bin/env python3

import sys
import os

import argparse

from numpy import *
from pylab import *

import matplotlib.pyplot as plt
import matplotlib.animation as ani

from lib.parameters import Parameters
from lib.pdp3_plot_builder import PDP3PlotBuilder

from lib.pdp_3e_view_builder import Pdp3EViewBuilder

## README:
## color map reference: https://matplotlib.org/examples/color/colormaps_reference.html
## mathtext reference:  https://matplotlib.org/users/mathtext.html

class Pdp3Images(Pdp3EViewBuilder):
    def __init__(self, cfg):
        self.images_path = '.'
        super(Pdp3Images, self).__init__(cfg)

    def setup_3e_view(self, view):
        '''
        initialize plot figure and subplots with preset object fields
        '''
        super(Pdp3Images, self).setup_3e_view()

        if view:
            self._plot_builder.figure.show()

    def create_view_with_3_plots(self, view=False, write=True):
        '''
        '''


        fpf = self._cfg.frames_per_file
        sr = self._cfg.r_grid_count
        sz = self._cfg.z_grid_count

        data_file_e_r = os.path.join(self._cfg.data_path, self.data_file_e_r_pattern)
        data_file_e_z = os.path.join(self._cfg.data_path, self.data_file_e_z_pattern)
        data_file_bunch_density = os.path.join(self._cfg.data_path, self.data_file_e_bunch_density_pattern)

        for k in range(self.start_data_set, self.end_data_set+1):
            tstart = self.start_frame if k == self.start_data_set else k*fpf
            tend = self.end_frame if k == self.end_data_set else ((k+1)*fpf)
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
            fidh_bunch_density = open(data_file_bunch_density + str(k), 'r')

            h_field_e_r = fromfile(fidh_e_r, dtype=float, count=sr*sz*fpf, sep=' ')
            h_field_e_z = fromfile(fidh_e_z, dtype=float, count=sr*sz*fpf, sep=' ')
            h_field_bunch_density = fromfile(fidh_bunch_density, dtype=float, count=sr*sz*fpf, sep=' ')

            ## Close data files
            fidh_e_r.close()
            fidh_e_z.close()
            fidh_bunch_density.close()

            for t in range(tstart, tend):
                local_step = t % fpf

                print("Processing frame %d" % (local_step))

                image_file_name = os.path.join(self.images_path, 'image_' + str(k) + '_' + str(local_step) + '.png')

                rstart = sr*sz*local_step
                rend = sr*sz*(local_step+1)

                try:
                    self._plot_builder.fill_image_with_data(
                        self.E_z_plot_name,
                        h_field_e_r[rstart:rend])

                    self._plot_builder.fill_image_with_data(
                        self.E_r_plot_name,
                        h_field_e_z[rstart:rend])

                    self._plot_builder.fill_image_with_data(
                        self.E_bunch_density_plot_name,
                        h_field_bunch_density[rstart:rend])

                except ValueError: ## skip frame, when data is inconsistent
                    break

                if write:
                    self._plot_builder.figure.savefig(image_file_name) # save the figure to file
                if view:
                    self._plot_builder.redraw()


def main():
    ## configure RC properties
    plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

    ####
    parser = argparse.ArgumentParser(description='Tool for making single images and its series from PDP3 modelled data.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    default_data_set_range = [0, 10000]

    parser.add_argument('--images-path', type=str,
                        help='Full path to images. Default <path/to/parameters.xml>')

    parser.add_argument('--timestamp', type=float, help='Timestamp to generate image at')

    parser.add_argument('--time-range', type=str, help='Time range. Can be overriden by --timestamp')

    parser.add_argument('--data-set-range', type=str,
                        help='''Range of data files set (e.g. 2:10 is E_r2 to Er_10, E_z2 to E_z10 and so on).
                        Can be overriden by --time-range and --timestamp. Default %s'''
                        % ':'.join(map(str, default_data_set_range)))

    parser.add_argument('--cmap', type=str,
                        help='''Use custom colormap. Default %s.
                        Reference: https://matplotlib.org/examples/color/colormaps_reference.html''' % 'gray',
                        default='gray')

    parser.add_argument('--beam-scale-factor', type=float,
                        help='''Beam density setting automatically, but you can set scale factor to sets,
                        where initial bunch density should be placed in color range''',
                        default=2)

    parser.add_argument('--clim-e-r', type=str,
                        help='Color limit range for Electrical field longitual component')
    parser.add_argument('--clim-e-z', type=str,
                        help='Color limit range for Electrical field radial component')

    parser.add_argument('--dry-run', action='store_true', help='Do not write anything. Just for debug')

    parser.add_argument('--view', action='store_true', default=False,
                        help='View animation as well as write it to file')

    args = parser.parse_args()

    view=False
    write=True

    if args.view:
        view = True
    if args.dry_run:
        write = False

    # check if config file exists
    if os.path.isfile(args.properties_path):
        ## initialize config
        config = Parameters(args.properties_path)
        images = Pdp3Images(config)

        ################################################################################################
        #################### configure plot and view parameters #######################################
        ################################################################################################
        images.clim_e_field_r = list(map(float, args.clim_e_r.split(':'))) if args.clim_e_r else [-config.clim_estimation, config.clim_estimation]
        images.clim_e_field_z = list(map(float, args.clim_e_z.split(':'))) if args.clim_e_z else [-config.clim_estimation, config.clim_estimation]
        images.cmap = args.cmap
        images.clim_e_field_beam_scale_factor = args.beam_scale_factor


        if args.timestamp:
            images.start_data_set, images.start_frame = config.get_file_frame_number_by_timestamp(args.timestamp)
            images.end_data_set, view.end_frame = view.start_data_set, view.start_frame + 1
        elif args.time_range:
            time_range = list(map(float, args.time_range.split(':')))
            images.start_data_set, images.start_frame = config.get_file_frame_number_by_timestamp(time_range[0])
            images.end_data_set, images.end_frame = config.get_file_frame_number_by_timestamp(time_range[1])
        elif args.data_set_range:
            data_set_range = list(map(int, args.data_set_range.split(':')))
            images.start_data_set = data_set_range[0]
            images.end_data_set = data_set_range[1]
        else:
            images.start_data_set = default_data_set_range[0]
            images.end_data_set = default_data_set_range[1]


        images.images_path = config.config_path if not args.images_path else args.images_path
        ################################################################################################
        ################################################################################################
        ################################################################################################

        images.setup_3e_view(view)
        images.create_view_with_3_plots(view, write)
        if view:
            input("Press 'Return' to exit ")
    else:
        print("Configuration file `%s' does not exists. Exiting" % args.properties_path)
        exit(1)

## call main function
if __name__ == "__main__":
    main()
