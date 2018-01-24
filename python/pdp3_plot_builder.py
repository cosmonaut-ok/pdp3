#!/usr/bin/env python3

import sys
import os

current_dir = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_dir)

from parameters import Parameters

from numpy import *
from pylab import *

import matplotlib.pyplot as plt

class PDP3PlotBuilder:
    def __init__(self, cfg):
        self.__cfg = cfg
        self.__subplots = {}
        self.__images = {}
        self.figure = plt.figure()

        # ## define data files
        # self.__data_file_e1 = os.path.join(self.__cfg.data_path, 'e1')
        # self.__data_file_e3 = os.path.join(self.__cfg.data_path, 'e3')
        # self.__data_file_rho_beam = os.path.join(self.__cfg.data_path, 'rho_beam')

        # # %% set number of marks with names in X and Y axes
        # x_ticks_number = 10;
        # y_ticks_number = 4;
        # # %% Titles and Ticks
        # self.__x_axe_title = 'Z(m)';
        # self.__y_axe_title = 'R(m)';
        # self.__x_tick_range = linspace(0, self.__cfg.z_size, x_ticks_number) # we need 10 (or x_ticks_number) ticks
        # self.__x_tick_gird_size = linspace(0, self.__cfg.z_grid_count, x_ticks_number) # from 0 to x_tick_max. it's required
        # self.__y_tick_range = linspace(0, self.__cfg.r_size, y_ticks_number) # to convert gird to real size (meters)
        # self.__y_tick_gird_size = linspace(0, self.__cfg.r_grid_count, y_ticks_number) # Same for X and Y axes


        # self.video_codec = 'mjpeg'
        # self.video_fps = 30
        # self.video_dpi = 100

        # self.__figure = plt.figure(figsize=(10.50, 7.00), dpi=self.video_dpi)

    # def __setup_subplot(self, subplot_number, clim, interpolation_type):
    #     a1 = self.__figure.add_subplot(subplot_number)
    #     a1.set_aspect('equal')
    #     ax1 = a1.get_xaxis()
    #     ay1 = a1.get_yaxis()
    #     a1.invert_yaxis()
    #     im1 = a1.imshow(rand(self.__cfg.r_grid_count,self.__cfg.z_grid_count),cmap=self.cmap,interpolation=interpolation_type)
    #     im1.set_clim(clim)
    #     return im1

    def setup_figure(self, width=10.5, height=7, dpi=100):
        self.figure.set_size_inches([width, height])
        self.figure.dpi = dpi

    def add_subplot_with_image(self, name, subplot_number, cmap='gray', clim=[-1, 1]):
        interpolation_type = 'nearest'
        subplot = self.figure.add_subplot(subplot_number)
        self.__subplots[name] = subplot
        #
        subplot.set_aspect('equal')
        #     ax1 = a1.get_xaxis()
        #     ay1 = a1.get_yaxis()
        subplot.invert_yaxis()
        ## initialize image with random data, normalized to clim range
        # initial_data = rand(self.__cfg.r_grid_count, self.__cfg.z_grid_count)*(clim[1]-clim[0])-((clim[0]+clim[1])/2)
        initial_data = zeros([self.__cfg.r_grid_count, self.__cfg.z_grid_count])
        image = subplot.imshow(initial_data,
                               cmap=cmap,
                               interpolation=interpolation_type)
        image.set_clim(clim)

        cbar = self.figure.colorbar(image, ticks=[-1, 0, 1])
        cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar

        self.__images[name] = image

        return subplot, image

    def get_subplot(self, name):
        self.__subplots[name]

    def get_image(self, name):
        self.__images[name]

    def fill_image_with_data(self, name, data):
        sr = self.__cfg.r_grid_count
        sz = self.__cfg.z_grid_count
        data_len = len(data)
        if data_len != sr*sz:
            raise ValueError('data array length is not equal to grid dimensions multiplication: %i X %i' % (sr,sz))

        im = self.__images[name]

        data_matrix = flipud(reshape(data, (sr, sz)))
        im.set_data(data_matrix)

    def redraw(self):
        self.figure.canvas.draw_idle()
        # self.__plot_builder.canvas.start_event_loop(1)
        self.figure.canvas.flush_events()
