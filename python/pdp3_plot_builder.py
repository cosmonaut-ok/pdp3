#!/usr/bin/env python3

import sys
import os

current_dir = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_dir)

from parameters import Parameters

from numpy import *
from pylab import *

import matplotlib.pyplot as plt
# from matplotlib import rc

class PDP3PlotBuilder:
    def __init__(self, cfg):
        self.__cfg = cfg
        self.__subplots = {}
        self.__images = {}
        self.figure = plt.figure()

        self.x_tick_count = 10;
        self.y_tick_count = 4;

    def get_subplot(self, name):
        return self.__subplots[name]

    def get_image(self, name):
        return self.__images[name]

    def setup_figure(self, width=10.5, height=7, dpi=100):
        self.figure.set_size_inches([width, height])
        self.figure.dpi = dpi
        # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
        # rc('text', usetex=True)

    def add_subplot_with_image(self, name, subplot_number, x_axe_label='X', y_axe_label='Y', cmap='gray', clim=[-1, 1]):
        # x_axe_title = 'Z(m)';
        # y_axe_title = 'R(m)';
        # colorbar_title = 'V/m';

        interpolation_type = 'nearest'
        subplot = self.figure.add_subplot(subplot_number)
        #
        subplot.set_aspect('equal')
        subplot.invert_yaxis()
        ## initialize image with random data, normalized to clim range
        # initial_data = rand(self.__cfg.r_grid_count, self.__cfg.z_grid_count)*(clim[1]-clim[0])-((clim[0]+clim[1])/2)
        initial_data = zeros([self.__cfg.r_grid_count, self.__cfg.z_grid_count])
        image = subplot.imshow(initial_data,
                               cmap=cmap,
                               interpolation=interpolation_type)
        image.set_clim(clim)

        self.__subplots[name] = subplot
        self.__images[name] = image

        return subplot, image

    def setup_subplot(self, name, title=None, x_axe_label='X', y_axe_label='Y'):
        axes = self.__subplots[name]

        _title = title or name

        x_tick_range = linspace(0, self.__cfg.z_size, self.x_tick_count+1) # we need 10 (or x_ticks_number) ticks
        x_tick_grid_size = linspace(0, self.__cfg.z_grid_count, self.x_tick_count+1) # from 0 to x_tick_max. it's required
        y_tick_range = linspace(0, self.__cfg.r_size, self.y_tick_count+1) # to convert gird to real size (meters)
        y_tick_grid_size = linspace(0, self.__cfg.r_grid_count, self.y_tick_count+1) # Same for X and Y axes

        axes.set_title(_title)
        axes.set_xlabel(x_axe_label)
        axes.set_ylabel(y_axe_label)

        axes.set_xticks(x_tick_grid_size)
        axes.set_yticks(y_tick_grid_size)

        axes.set_xticklabels(x_tick_range)
        axes.set_yticklabels(y_tick_range)

        axes.spines['top'].set_visible(False)
        axes.spines['right'].set_visible(False)

    def add_colorbar(self, name, title=None, ticks=[-1, 1], ticklabels=[-1, 1]):
        image = self.get_image(name)
        cbar = self.figure.colorbar(image)

        _title = title or name

        cbar.set_ticks(ticks)
        cbar.set_ticklabels(ticklabels)
        cbar.set_label(_title)
        # cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar

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
