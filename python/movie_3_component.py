#!/usr/bin/env python3

import sys
import os

current_dir = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_dir)

from parameters import Parameters

from pdp3_plot_builder import PDP3PlotBuilder

import argparse

from numpy import *
from pylab import *

# from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.animation as ani

class Pdp3Movie:
    def __init__(self, cfg): # xml_config_file, video_file=None, clim_er=[0,1], clim_ez=[0,1]):
        self.__cfg = cfg # Parameters(xml_config_file, clim_er, clim_ez, video_file)

        ## define data files
        self.__data_file_e1 = os.path.join(self.__cfg.data_path, 'e1')
        self.__data_file_e3 = os.path.join(self.__cfg.data_path, 'e3')
        self.__data_file_rho_beam = os.path.join(self.__cfg.data_path, 'rho_beam')

        # %% set number of marks with names in X and Y axes
        # x_ticks_number = 10;
        # y_ticks_number = 4;
        # %% Titles and Ticks
        # self.__x_axe_title = 'Z(m)';
        # self.__y_axe_title = 'R(m)';
        # self.__x_tick_range = linspace(0, self.__cfg.z_size, x_ticks_number) # we need 10 (or x_ticks_number) ticks
        # self.__x_tick_gird_size = linspace(0, self.__cfg.z_grid_count, x_ticks_number) # from 0 to x_tick_max. it's required
        # self.__y_tick_range = linspace(0, self.__cfg.r_size, y_ticks_number) # to convert gird to real size (meters)
        # self.__y_tick_gird_size = linspace(0, self.__cfg.r_grid_count, y_ticks_number) # Same for X and Y axes

        self.cmap = 'gray'
        self.video_codec = 'mjpeg'
        self.video_fps = 30
        self.video_dpi = 100
        self.video_bitrate=32000

        self.__plot_builder = PDP3PlotBuilder(self.__cfg)
        self.__plot_builder.setup_figure()

        self.E_r_plot_name = 'E_r'
        self.E_z_plot_name = 'E_z'
        self.E_rho_beam_plot_name = 'RHO_beam'

        self.__plot_builder.add_subplot_with_image(self.E_r_plot_name, 311,
                                                   cmap=self.cmap, clim=self.__cfg.clim_e_field_r)
        self.__plot_builder.add_subplot_with_image(self.E_z_plot_name, 312,
                                                   cmap=self.cmap, clim=self.__cfg.clim_e_field_z)
        self.__plot_builder.add_subplot_with_image(self.E_rho_beam_plot_name, 313,
                                                   cmap=self.cmap, clim=self.__cfg.clim_e_field_bunch)

    def create_movie_with_3_plots(self):
        # im1, im2, im3 = self.__setup_figure()
        self.__plot_builder.figure.show()
        ###########
        fpf = self.__cfg.frames_per_file
        sr = self.__cfg.r_grid_count
        sz = self.__cfg.z_grid_count

        FFMpegWriter = ani.writers['ffmpeg']
        metadata = dict(title='Movie Test', artist='Matplotlib',
                        comment='Movie support!')
        writer = FFMpegWriter(fps=self.video_fps,
                              metadata=metadata,
                              codec=self.video_codec,
                              bitrate=self.video_bitrate)

        with writer.saving(self.__plot_builder.figure, self.__cfg.movie_file, self.video_dpi):
            for k in range(0, fpf):
                tstart = k*fpf
                tend = ((k+1)*fpf-1)
                i = 1;

                if not os.path.isfile(self.__data_file_e1 + str(k)) \
                   or not os.path.isfile(self.__data_file_e3 + str(k)) \
                   or not os.path.isfile(self.__data_file_rho_beam + str(k)):
                    print('No more data files exists. Exiting')
                    return

                print("Loading files set", k)
                ## Open data files
                fidh_e1 = open(self.__data_file_e1 + str(k), 'r')
                fidh_e3 = open(self.__data_file_e3 + str(k), 'r')
                fidh_rho_beam = open(self.__data_file_rho_beam + str(k), 'r')

                print(sr*sz*fpf)
                h_field_e1 = fromfile(fidh_e1, dtype=float, count=sr*sz*fpf, sep=' ')
                h_field_e3 = fromfile(fidh_e3, dtype=float, count=sr*sz*fpf, sep=' ')
                h_field_rho_beam = fromfile(fidh_rho_beam, dtype=float, count=sr*sz*fpf, sep=' ')

                ## Close data files
                fidh_e1.close()
                fidh_e3.close()
                fidh_rho_beam.close()

                for t in range(tstart, tend):
                    local_step = t % fpf

                    print("Processing frame %d" % (local_step))

                    rstart = sr*sz*local_step+1
                    rend = sr*sz*(local_step+1)+1
                    try:
                        self.__plot_builder.fill_image_with_data(
                            self.E_z_plot_name,
                            h_field_e1[rstart:rend])

                        self.__plot_builder.fill_image_with_data(
                            self.E_r_plot_name,
                            h_field_e3[rstart:rend])

                        self.__plot_builder.fill_image_with_data(
                            self.E_rho_beam_plot_name,
                            h_field_rho_beam[rstart:rend])

                    except ValueError: ## skip frame, when data is inconsistent
                        break

                    writer.grab_frame()

                    # self.__plot_builder.redraw()

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    parser.add_argument('--video_file', type=str,
                        help='Full path to output video file')

    # parser.add_argument('clim_e1', metavar='properties_path', type=str,
    #                     help='Full path to properties.xml')

    ## video_file=None, file_delta=100, clim_e1=[0,1], clim_e3=[0,1], clim_rho_beam=[-1e-7, 0]):


    parser.add_argument('--clim_e1', type=str,
                        help='Color limit range for e1. Default: 0:1. Not implemented')
    parser.add_argument('--clim_e3', type=str,
                        help='Color limit range for e3. Default: 0:1. Not implemented')
    parser.add_argument('--clim_rho_beam', type=str,
                        help='Color limit range for rho_beam. Default: -1e-7:0. Not implemented')

    args = parser.parse_args()
    # print(args)

    clim_e1 = list(map(int, args.clim_e1.split(':'))) if args.clim_e1 else [-1e5, 1e5]
    clim_e3 = list(map(int, args.clim_e3.split(':'))) if args.clim_e3 else [-1e5, 1e5]
    # clim_rho_beam = list(map(int, args.clim_rho_beam.split(':'))) if args.clim_rho_beam else [-1e-7, 0]
    file_delta = 100

    ## initialize config
    config = Parameters(args.properties_path, args.video_file, clim_e1, clim_e3)

    movie = Pdp3Movie(config)
    movie.create_movie_with_3_plots()



if __name__ == "__main__":
    main()
