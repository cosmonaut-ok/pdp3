#!/usr/bin/env python3

import sys
import os

import argparse

from numpy import *
from pylab import *

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as ani

from lib.parameters import Parameters
from lib.pdp3_plot_builder import PDP3PlotBuilder

from lib.pdp_3e_view_builder import Pdp3EViewBuilder

## README:
## color map reference: https://matplotlib.org/examples/color/colormaps_reference.html
## mathtext reference:  https://matplotlib.org/users/mathtext.html

class Pdp3View(Pdp3EViewBuilder):
    def __init__(self, cfg):
        super(Pdp3View, self).__init__(cfg)

    def setup_3e_view(self):
        '''
        initialize plot figure and subplots with preset object fields
        '''
        super(Pdp3View, self).setup_3e_view()

        mgr = self._plot_builder.figure.canvas.manager
        if matplotlib.get_backend() == 'Qt5Agg':  # 'Qt4' backend
            mgr.window.showMaximized()
        elif matplotlib.get_backend() == 'Qt4Agg':  # 'Qt4' backend
            mgr.window.showMaximized()
        elif matplotlib.get_backend() == 'WxAgg':  # 'WxAgg' backend
            mgr.frame.Maximize(True)
        elif matplotlib.get_backend() == 'TKAgg':  # 'TKAgg' backend
            matplotlib.frame.Maximize(True)

        self._plot_builder.figure.show()

def main():
    ## configure RC properties
    plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

    ####
    parser = argparse.ArgumentParser(description='Tool for visualization of data, generated by PDP3.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    default_data_set_range = [0, 10000]
    default_clim = [-1e5, 1e5]

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
                        help='Color limit range for Electrical field longitual component. Default %s'
                        % ':'.join(map(str, default_clim)))
    parser.add_argument('--clim-e-z', type=str,
                        help='Color limit range for Electrical field radial component. Default %s'
                        % ':'.join(map(str, default_clim)))

    args = parser.parse_args()

    # check if config file exists
    if os.path.isfile(args.properties_path):
        ## initialize config
        config = Parameters(args.properties_path)
        view = Pdp3View(config)

        ################################################################################################
        #################### configure plot and view parameters #######################################
        ################################################################################################
        view.clim_e_field_r = list(map(float, args.clim_e_r.split(':'))) if args.clim_e_r else default_clim
        view.clim_e_field_z = list(map(float, args.clim_e_z.split(':'))) if args.clim_e_z else default_clim
        view.cmap = args.cmap
        view.clim_e_field_beam_scale_factor = args.beam_scale_factor

        if args.timestamp:
            view.start_data_set, view.start_frame = config.get_file_frame_by_timestamp(args.timestamp)
            view.end_data_set, view.end_frame = config.get_file_frame_by_timestamp(args.timestamp)
            view.end_frame = view.end_frame + 1
        elif args.time_range:
            time_range = list(map(float, args.time_range.split(':')))
            view.start_data_set, view.start_frame = config.get_file_frame_by_timestamp(time_range[0])
            view.end_data_set, view.end_frame = config.get_file_frame_by_timestamp(time_range[1])
        elif args.data_set_range:
            data_set_range = list(map(int, args.data_set_range.split(':')))
            view.start_data_set = data_set_range[0]
            view.end_data_set = data_set_range[1]
        else:
            view.start_data_set = default_data_set_range[0]
            view.end_data_set = default_data_set_range[1]
        ################################################################################################
        ################################################################################################
        ################################################################################################

        view.setup_3e_view()
        view.create_view_with_3_plots()
        input("Press 'Return' to exit ")
    else:
        print("Configuration file `%s' does not exists. Exiting" % args.properties_path)
        exit(1)

## call main function
if __name__ == "__main__":
    main()
