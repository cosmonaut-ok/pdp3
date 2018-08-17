#!/usr/bin/env python3

import sys
import os

current_dir = os.path.dirname(os.path.abspath('__file__'))
sys.path.append(current_dir)

from lib.parameters import Parameters

from numpy import *
from pylab import *

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import rc

from mpl_toolkits.axes_grid1 import make_axes_locatable

class PDP3PlotBuilder:
    def __init__(self, cfg):
        self.__cfg = cfg
        self.__subplots = {}
        self.__images = {}
        self.figure = plt.figure()

        # self.figure.subplots_adjust(hspace=1.8)

        self.x_tick_count = 20;
        self.y_tick_count = 8;

        self.x_plot_size = self.__cfg.z_grid_count
        self.y_plot_size = self.__cfg.r_grid_count

        self.x_tick_start = 0
        self.x_tick_end = self.__cfg.z_size
        self.y_tick_start = 0
        self.y_tick_end = self.__cfg.r_size

        self.aspect = 'equal'

    def get_subplot(self, name):
        ''' Object collects axes (subplots) and images for quick access
        this function allows to quick access to subplot by name
        '''
        return self.__subplots[name]

    def get_image(self, name):
        ''' Object collects axes (subplots) and images for quick access
        this function allows to quick access to image by name
        (same, as subplot, to what it bount)
        '''
        return self.__images[name]

    def setup_figure(self, width=10.5, height=7, font_size=10, dpi=100):
        ''' set figure parameters '''
        self.figure.set_size_inches([width, height])
        self.figure.dpi = dpi
        # rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':font_size})
        rc('font',**{'size':font_size})
        # rc('text', usetex=True)

    def add_subplot(self, name, subplot_number):
        ''' add subplot '''
        subplot = self.figure.add_subplot(subplot_number)
        #
        subplot.set_aspect(self.aspect)
        # subplot.invert_yaxis()

        self.__subplots[name] = subplot

        return subplot


    def add_subplot_with_image(self, name, subplot_number, x_axe_label='X', y_axe_label='Y', cmap='gray', clim=[-1, 1]):
        ''' add subplot and place there image with same name and initial data (zeros)
        subplot and image parameters configures, using 'parameters' module
        also, this function adds subplot and image objects to plot builder's dictionary
        for quick access

        FYI: cmap reference: https://matplotlib.org/examples/color/colormaps_reference.html
        '''

        interpolation_type = 'nearest'
        subplot = self.add_subplot(name, subplot_number)

        ## initialize image with random data, normalized to clim range
        # initial_data = rand(self.y_plot_size, self.x_plot_size)*(clim[1]-clim[0])-((clim[0]+clim[1])/2)
        initial_data = zeros([self.y_plot_size, self.x_plot_size])
        image = subplot.imshow(initial_data,
                               cmap=cmap,
                               interpolation=interpolation_type)
        image.set_clim(clim)

        self.__subplots[name] = subplot
        self.__images[name] = image

        return subplot, image

    def setup_subplot(self, name, title=None, x_axe_label='X', y_axe_label='Y', tickbox=False, grid=False, position=None):
        ''' wrapper for configuration subplot properties '''
        axes = self.get_subplot(name)

        _title = title or name

        # tick labels, that shows __real__ model space dimensions
        # translates from grid_size
        x_tick_range = around(linspace(self.x_tick_start, self.x_tick_end, self.x_tick_count+1), 2)
        y_tick_range = around(linspace(self.y_tick_start, self.y_tick_end, self.y_tick_count+1), 2)

        # ticks, that sets grid dimensions, required for data placement
        x_tick_grid_size = linspace(0, self.x_plot_size, self.x_tick_count+1)
        y_tick_grid_size = linspace(self.y_plot_size, 0, self.y_tick_count+1)

        # set axis properties
        axes.set_title(_title)

        axes.set_xlabel(x_axe_label)
        axes.set_ylabel(y_axe_label, rotation=45)

        axes.set_xticks(x_tick_grid_size)
        axes.set_yticks(y_tick_grid_size)

        axes.set_xticklabels(x_tick_range)
        axes.set_yticklabels(y_tick_range)

        # axes.xticks(rotation=90)
        axes.spines['top'].set_visible(tickbox)
        axes.spines['right'].set_visible(tickbox)

        # set label on every 4th grid
        for label in [x for i,x in enumerate(axes.xaxis.get_ticklabels()) if i%2 != 0]:
            label.set_visible(False)
        for label in [x for i,x in enumerate(axes.yaxis.get_ticklabels()) if i%2 != 0]:
            label.set_visible(False)

        axes.grid(grid)

        if position:
            axes.set_position(position)


    def add_colorbar(self, name, title=None, ticks=[-1, 1], ticklabels=None, font_size=6, size="2%", position="right"):
        ''' add colorbar to selected image and set it's properties '''

        image = self.get_image(name)
        ax = self.get_subplot(name)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes(position, size=size, pad=0.05)

        cbar = self.figure.colorbar(image, cax=cax) # , ax = axes)

        label_count = 3

        __title = title or name
        __ticks = linspace(ticks[0], ticks[1], label_count)
        _ticklabels = linspace(ticklabels[0], ticklabels[1], label_count) if ticklabels else __ticks

        def format_s(x):
            return('%.0e' % x)

        __ticklabels = list(map(format_s, _ticklabels))

        cbar.set_label(__title, rotation=45)
        cbar.set_ticks(__ticks)
        cbar.set_ticklabels(__ticklabels)
        cbar.ax.tick_params(labelsize=font_size)



        ## set sci limit to 100. More, than 1000 should be printed in format XeY
        # cbar.formatter.set_powerlimits([-0, 0])

        # fmt = matplotlib.ticker.ScalarFormatter(useMathText=True)
        # fmt.set_scientific(True)
        # fmt.set_powerlimits([0,0])
        # cbar.ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%e'))
        # cbar.ax.yaxis.set_major_formatter(fmt)

        # cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar

    def fill_image_with_data(self, name, data):
        ''' fill image with 1-dimenstional data array, (usually got from model files)
        sized by grid dimensions (n_grid_r, n_grid_z in parameters.xml) multiplication
        used to fill model image with data '''
        sr = self.y_plot_size
        sz = self.x_plot_size
        data_len = len(data)
        if data_len != sr*sz:
            raise ValueError('data array length is not equal to grid dimensions multiplication: %i X %i' % (sr,sz))

        im = self.get_image(name)

        data_matrix = flipud(reshape(data, (sr, sz)))
        im.set_data(data_matrix)

    def redraw(self):
        ''' Redraw figure (can be used for animation and video writing) '''
        self.figure.canvas.draw_idle()
        # self.__plot_builder.canvas.start_event_loop(1)
        self.figure.canvas.flush_events()
