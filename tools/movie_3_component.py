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

class Pdp3Movie(Pdp3EViewBuilder):
    def __init__(self, cfg):

        self.video_file = 'field_movie.avi'

        super(Pdp3Movie, self).__init__(cfg)

    def setup_3e_view(self, view):
        '''
        initialize plot figure and subplots with preset object fields
        '''
        super(Pdp3Movie, self).setup_3e_view()

        if view:
            self._plot_builder.figure.show()

    def create_view_with_3_plots(self, view=False, write=True):
        '''
        '''
        FFMpegWriter = ani.writers['ffmpeg']
        metadata = dict(title='Movie Test', artist='Matplotlib',
                        comment='Movie support!')
        writer = FFMpegWriter(fps=self._cfg.video_fps,
                              metadata=metadata,
                              codec=self._cfg.video_codec,
                              bitrate=self._cfg.video_bitrate)

        ## dirty hack to aviod real file writing
        if write:
            movie_file = self.video_file
        else:
            movie_file = '/dev/null'

        fpf = self._cfg.frames_per_file
        sr = self._cfg.r_grid_count
        sz = self._cfg.z_grid_count

        data_file_e_r = os.path.join(self._cfg.data_path, self.data_file_e_r_pattern)
        data_file_e_z = os.path.join(self._cfg.data_path, self.data_file_e_z_pattern)
        data_file_bunch_density = os.path.join(self._cfg.data_path, self.data_file_e_bunch_density_pattern)

        with writer.saving(self._plot_builder.figure, movie_file, self._cfg.video_dpi):
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
                    writer.grab_frame()
                if view:
                    self._plot_builder.redraw()

def main():
    ## configure RC properties
    plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

    ####
    parser = argparse.ArgumentParser(description='Tool for creating movie file from PDP3 modelled data.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    default_data_set_range = [0, 10000]
    default_clim = [-1e5, 1e5]

    parser.add_argument('--video-file', type=str,
                        help='Full path to output video file. Default <path/to/parameters.xml>/field_movie.avi')

    parser.add_argument('--time-range', type=str, help='Time range')

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

    ## required parameter
    # if not (args.data_set_range or args.time_range or args.timestamp):
    #     parser.error('one of --data-set-range or --time-range or --timestamp required')

    # check if config file exists
    if os.path.isfile(args.properties_path):
        ## initialize config
        config = Parameters(args.properties_path)
        movie = Pdp3Movie(config)

        ################################################################################################
        #################### configure plot and view parameters #######################################
        ################################################################################################
        movie.video_file = os.path.join(config.config_path, 'field_movie.avi') if not args.video_file else args.video_file
        movie.clim_e_field_r = list(map(float, args.clim_e_r.split(':'))) if args.clim_e_r else default_clim
        movie.clim_e_field_z = list(map(float, args.clim_e_z.split(':'))) if args.clim_e_z else default_clim
        movie.cmap = args.cmap
        movie.clim_e_field_beam_scale_factor = args.beam_scale_factor

        if args.time_range:
            time_range = list(map(float, args.time_range.split(':')))
            movie.start_data_set, movie.start_frame = config.get_file_frame_by_timestamp(time_range[0])
            movie.end_data_set, movie.end_frame = config.get_file_frame_by_timestamp(time_range[1])
        elif args.data_set_range:
            data_set_range = list(map(int, args.data_set_range.split(':')))
            movie.start_data_set = data_set_range[0]
            movie.end_data_set = data_set_range[1]
        else:
            movie.start_data_set = default_data_set_range[0]
            movie.end_data_set = default_data_set_range[1]
        ################################################################################################
        ################################################################################################
        ################################################################################################

        movie.setup_3e_view(view)
        movie.create_view_with_3_plots(view, write)
        input("Press 'Return' to exit ")
    else:
        print("Configuration file `%s' does not exists. Exiting" % args.properties_path)
        exit(1)

## call main function
if __name__ == "__main__":
    main()
