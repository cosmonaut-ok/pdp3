#!/usr/bin/env python3

import sys
import os

current_dir = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_dir)

from parameters import Parameters

import argparse

from numpy import *
from pylab import *

# from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.animation as ani

class Pdp3Movie:
    def __init__(self, xml_config_file, video_file=None, clim_e1=[0,1], clim_e3=[0,1]):
        self.__cfg = Parameters(xml_config_file, clim_e1, clim_e3, video_file)

        ## define data files
        self.__data_file_e1 = os.path.join(self.__cfg.data_path, 'e1')
        self.__data_file_e3 = os.path.join(self.__cfg.data_path, 'e3')
        self.__data_file_rho_beam = os.path.join(self.__cfg.data_path, 'rho_beam')

        # %% set number of marks with names in X and Y axes
        x_ticks_number = 10;
        y_ticks_number = 4;
        # %% Titles and Ticks
        self.__x_axe_title = 'Z(m)';
        self.__y_axe_title = 'R(m)';
        self.__x_tick_range = linspace(0, self.__cfg.z_size, x_ticks_number) # we need 10 (or x_ticks_number) ticks
        self.__x_tick_gird_size = linspace(0, self.__cfg.z_grid_count, x_ticks_number) # from 0 to x_tick_max. it's required
        self.__y_tick_range = linspace(0, self.__cfg.r_size, y_ticks_number) # to convert gird to real size (meters)
        self.__y_tick_gird_size = linspace(0, self.__cfg.r_grid_count, y_ticks_number) # Same for X and Y axes

        self.cmap = 'gray'
        self.video_codec = 'mjpeg'
        self.video_fps = 30
        self.video_dpi = 100

        self.__figure = plt.figure(figsize=(10.50, 7.00), dpi=self.video_dpi)


    def __setup_subplot(self, subplot_number, clim, interpolation_type):
        a1 = self.__figure.add_subplot(subplot_number)
        a1.set_aspect('equal')
        ax1 = a1.get_xaxis()
        ay1 = a1.get_yaxis()
        a1.invert_yaxis()
        im1 = a1.imshow(rand(self.__cfg.r_grid_count,self.__cfg.z_grid_count),cmap=self.cmap,interpolation=interpolation_type)
        im1.set_clim(clim)
        return a1

    def __setup_figure(self):

        interpolation_type = 'nearest'

        self.__figure.set_size_inches([10,5])

        # tight_layout()

        a1 = self.__setup_subplot(311, self.__cfg.clim_e_field_r, interpolation_type)
        a2 = self.__setup_subplot(312, self.__cfg.clim_e_field_z, interpolation_type)
        a3 = self.__setup_subplot(313, self.__cfg.clim_e_field_bunch, interpolation_type)
        return a1, a2, a3

    def create_movie_with_3_plots(self):
        a1, a2, a3 = self.__setup_figure()
        self.__figure.show()
        ###########
        fpf = self.__cfg.frames_per_file
        sr = self.__cfg.r_grid_count
        sz = self.__cfg.z_grid_count


        FFMpegWriter = ani.writers['ffmpeg']
        metadata = dict(title='Movie Test', artist='Matplotlib',
                        comment='Movie support!')
        writer = FFMpegWriter(fps=self.video_fps, metadata=metadata, codec=self.video_codec)

        with writer.saving(self.__figure, self.__cfg.movie_file, self.video_dpi):
            for k in range(0, fpf):
                tstart = k*fpf # k=2 -> 200
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

                    try:
                        h_field_shot1 = h_field_e1[sr*sz*local_step+1:sr*sz*(local_step+1)+1]
                        h_field_matrix1 = flipud(reshape(h_field_shot1, (sr, sz)))

                        h_field_shot2 = h_field_e3[sr*sz*local_step+1:sr*sz*(local_step+1)+1]
                        h_field_matrix2 = flipud(reshape(h_field_shot2, (sr, sz)))

                        h_field_shot3 = h_field_rho_beam[sr*sz*local_step+1:sr*sz*(local_step+1)+1]
                        h_field_matrix3 = flipud(reshape(h_field_shot3, (sr, sz)))
                    except:
                        break

                    interpolation_type = 'nearest'
                    a1.imshow(h_field_matrix1, cmap=self.cmap, interpolation=interpolation_type)
                    a2.imshow(h_field_matrix2, cmap=self.cmap, interpolation=interpolation_type)
                    a3.imshow(h_field_matrix3, cmap=self.cmap, interpolation=interpolation_type)
                    writer.grab_frame()
                self.__figure.canvas.draw()
                self.__figure.canvas.flush_events()


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
    print(args)

    clim_e1 = list(map(int, args.clim_e1.split(':'))) if args.clim_e1 else [0,1]
    clim_e3 = list(map(int, args.clim_e3.split(':'))) if args.clim_e3 else [0,1]
    clim_rho_beam = list(map(int, args.clim_rho_beam.split(':'))) if args.clim_rho_beam else [-1e-7, 0]
    file_delta = 100

    movie = Pdp3Movie(args.properties_path, video_file=args.video_file, clim_e1=clim_e1, clim_e3=clim_e3)
    movie.create_movie_with_3_plots()



if __name__ == "__main__":
    main()
