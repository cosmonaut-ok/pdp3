import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

from matplotlib import cm
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data


class PlotBuiler:
    def __init__(self, parameters, fig_color=None, fig_width=10.5, fig_height=7, fig_dpi=100,
                 font_family='sans-serif', font_name='DejaVu Sans', font_size=10):
        self.__parameters__ = parameters
        self.__subplots__ = []
        self.__images__ = []
        self.__figure__ = plt.figure(figsize=[fig_width, fig_height], dpi=fig_dpi, facecolor=fig_color)
        #
        self.__font_size__ = font_size
        self.__font_family__ = font_family
        self.__font_name__ = font_name
        rc('font',**{'family':font_family, font_family:[font_name], 'size':font_size})
        rc('text', usetex=True)

        self.aspect = 'equal'

        self.number_x_ticks = 20;
        self.number_y_ticks = 8;
        self.number_y_ticks = 20;
        self.number_cbar_ticks = 3;

        self.x_plot_size = self.__parameters__.number_z_grid
        self.y_plot_size = self.__parameters__.number_r_grid
        self.z_plot_size = self.y_plot_size

        self.x_tick_start, self.y_tick_start, self.z_tick_start = 0
        self.x_tick_end = self.__parameters__.z_size
        self.y_tick_end = self.__parameters__.r_size
        self.z_tick_end = self.y_tick_end

        self.default_image_interpolation = 'nearest'

    def get_figure(self):
        return(self.__figure__)

    def get_subplot(self, name):
        ''' Object collects axes (subplots) and images for quick access
        this function allows to quick access to subplot by name
        '''
        return(self.__subplots__[name])

    def get_image(self, name):
        ''' Object collects axes (subplots) and images for quick access
        this function allows to quick access to image by name
        (same, as subplot, to what it bount)
        '''
        return(self.__images[name])

##################################################################################################

    def __add_subplot_common__(self, name, number, title, x_axe_label, y_axe_label, z_axe_label,
                               tickbox, grid, position, is_invert_y_axe, projection=None):
        ''' common method to add any'''
        subplot = self.__figure__.add_subplot(subplot_number, projection=projection)
        self.__subplots__[name] = subplot

        subplot.set_aspect(self.aspect)
        if is_invert_y_axe: subplot.invert_yaxis()

        self.__subplots[name] = subplot

        # axes = self.get_subplot(name)

        __title = title or name

        # tick labels, that shows __real__ model space dimensions
        # translates from grid_size
        x_tick_range = around(linspace(self.x_tick_start, self.x_tick_end, self.number_x_ticks + 1), 2)
        y_tick_range = around(linspace(self.y_tick_start, self.y_tick_end, self.number_y_ticks + 1), 2)
        z_tick_range = around(linspace(self.z_tick_start, self.z_tick_end, self.number_z_ticks + 1), 2)

        # ticks, that sets grid dimensions, required for data placement
        x_tick_grid_size = linspace(0, self.x_plot_size, self.number_x_ticks + 1)
        y_tick_grid_size = linspace(0, self.y_plot_size, self.number_y_ticks + 1)
        z_tick_grid_size = linspace(0, self.z_plot_size, self.number_z_ticks + 1)

        # set axis properties
        axes.set_title(__title)

        subplot.set_xlabel(x_axe_label)
        subplot.set_ylabel(y_axe_label, rotation=45)
        if projection == '3d': subplot.set_zlabel(z_axe_label, rotation=45)

        subplot.set_xticks(x_tick_grid_size)
        subplot.set_yticks(y_tick_grid_size)
        if projection == '3d': subplot.set_zticks(z_tick_grid_size)

        subplot.set_xticklabels(x_tick_range)
        subplot.set_yticklabels(y_tick_range)
        if projection == '3d': axes.set_zticklabels(z_tick_range)

        # subplot.xticks(rotation=90)
        subplot.spines['top'].set_visible(tickbox)
        subplot.spines['right'].set_visible(tickbox)

        # set label on every 4th grid
        for label in [x for i,x in enumerate(subplot.xaxis.get_ticklabels()) if i%2 != 0]:
            label.set_visible(False)
        for label in [x for i,x in enumerate(subplot.yaxis.get_ticklabels()) if i%2 != 0]:
            label.set_visible(False)
        if projection == '3d':
            for label in [x for i,x in enumerate(subplot.zaxis.get_ticklabels()) if i%2 != 0]:
                label.set_visible(False)

        subplot.grid(grid)

        if position:
            subplot.set_position(position)

        return subplot


    def add_subplot_cartesian_2d(self, name, number, title=None, x_axe_label='X', y_axe_label='Y',
                      tickbox=False, grid=False, position=None):
        ''' add 2D subplot, cartesian projection '''

        subplot = self.__add_subplot_common__(name, number, title or name,
                                              x_axe_label, y_axe_label, None,
                                              tickbox, grid, position, True,
                                              projection=None)

        return(subplot)


    def add_subplot_cartesian_3d(self, name, number, title=None, x_axe_label='X', y_axe_label='Y', z_axe_label='Z',
                      tickbox=False, grid=False, position=None):
        ''' https://matplotlib.org/gallery/mplot3d/subplot3d.html '''
        subplot = self.__add_subplot_common__(name, number, title or name,
                                              x_axe_label, y_axe_label, z_axe_label,
                                              tickbox, grid, position, True,
                                              projection='3d')

        return(subplot)

    def add_image(self, subplot_name, data, cmap='gray', clim=[-1, 1], interpolation=None):
        '''
        cmap reference: https://matplotlib.org/examples/color/colormaps_reference.html
        '''
        if not interpolation: interpolation = self.default_image_interpolation
        subplot = self.get_subplot(subplot_name)

        if subplot:
            ## initialize image with random data, normalized to clim range
            # initial_data = rand(self.y_plot_size, self.x_plot_size)*(clim[1]-clim[0])-((clim[0]+clim[1])/2)
            initial_data = zeros([self.y_plot_size, self.x_plot_size])
            image = subplot.imshow(initial_data,
                                   cmap=cmap,
                                   aspect='auto',
                                   origin='lower',
                                   interpolation=interpolation)
            image.set_clim(clim)

            ####
            sr = self.y_plot_size
            sz = self.x_plot_size
            data_len = len(data)
            data_high = len(data[0])
            if data_len != self.x_plot_size:
                raise ValueError('data array length is not equal to grid X-dimension {}. The value was {}.'.format(self.x_plot_size, data_len))
            elif data_high != self.y_plot_size:
                raise ValueError('data array height is not equal to grid Y-dimension {}. The value was {}.'.format(self.y_plot_size, data_high))
            else:
                image.set_data(np.flipud(data))

            self.__images[name] = image

            return(image)
        else:
            raise Exception('There is no subplot, named {}'.format(subplot_name))


    def add_colorbar(self, subplot_name, image_name=None, title=None, ticks=[-1, 1], ticklabels=None,
                     font_size=None, size="2%", position="right"):
        image = self.get_image(image_name)
        if not font_size: font_size = self.font_size
        if not image_name: image_name = subplot_name
        if not title: title = subplot_name
        if image:
            ax = self.get_subplot(subplot_name)
            divider = make_axes_locatable(ax)
            cax = divider.append_axes(position, size=size, pad=0.05)

            cbar = self.figure.colorbar(image, cax=cax) # , ax = axes)

            __ticks = linspace(ticks[0], ticks[1], self.number_cbar_ticks)
            __ticklabels = linspace(ticklabels[0], ticklabels[1], self.number_cbar_ticks) if ticklabels else __ticks

            def format_s(x):
                return('%.0e' % x)

            __ticklabels = list(map(format_s, _ticklabels))

            cbar.set_label(title, rotation=45)
            cbar.set_ticks(__ticks)
            cbar.set_ticklabels(__ticklabels)
            cbar.ax.tick_params(labelsize=font_size)


    def redraw(self):
        ''' Redraw figure (can be used for animation and video writing) '''
        self.figure.canvas.draw_idle()
        # self.__plot_builder.canvas.start_event_loop(1)
        self.__figure__.canvas.flush_events()
