#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np

import matplotlib.animation as ani

from lib.parameters import Parameters
from lib.plot_builder import PlotBuilder
from lib.h5_reader import H5Reader
from lib.plain_reader import PlainReader

def run(config_file, clim_e_r, clim_e_z, rho_beam_scale, image_file=None,
        timestamp=0, cmap=None, dry_run=False, view=False, use_grid=False):

    ##  configuration options
    x_axis_label = r'$\mathit{Z (m)}$'
    y_axis_label = r'$\mathit{R (m)}$'
    cbar_axis_label = r'$\frac{V}{m}$'
    cbar_bunch_density_axis_label = r'$m^{-3}$'
    
    e_r_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Radial\enspace Component}\enspace(E_r)$'
    e_z_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Longitudal\enspace Component}\enspace(E_z)$'
    rho_beam_plot_name = r'$\mathbf{Electron\enspace Bunch\enspace Density}\enspace (\rho_{bunch})$'

    ## read configfile
    cfg = Parameters(config_file)

    # check if timestamp is correct
    if timestamp > cfg.end_time: raise IndexError("Timestamp is out of simulation range {}. The value was {}".format(cfg.end_time, timestamp))
    
    # calculate/update color limit values
    el_charge = 1.6e-19
    clim_rho_beam = [-(cfg.bunch_density * el_charge * rho_beam_scale), 0]
    clim_estimation = cfg.get_clim_estimation()

    if not clim_e_r: clim_e_r = [-clim_estimation, clim_estimation]
    if not clim_e_z: clim_e_z = [-clim_estimation, clim_estimation]
    if not rho_beam_scale: rho_beam_scale = 1

    # calculate/update video file path
    if not image_file: image_file = os.path.join(os.path.dirname(config_file), "image_{}.png".format(timestamp))

    # define reader (plain reader used)
    if not cfg.use_hdf5:
        reader = PlainReader(cfg.data_path, cfg.system_state_path,
                             [cfg.number_r_grid , cfg.number_z_grid],
                             cfg.frames_per_file, False)
    else:
        raise NotImplementedError('HDF5 support still not implemented')

    # define plot builder
    plot = PlotBuilder(cfg.number_z_grid, cfg.number_r_grid,
                       fig_color=cfg.figure_color, fig_width=cfg.figure_width,
                       fig_height=cfg.figure_height, fig_dpi=cfg.figure_dpi,
                       font_family=cfg.figure_font_family,
                       font_name=cfg.figure_font_name,
                       font_size=cfg.figure_font_size,
                       
                       x_ticklabel_end=cfg.z_size, y_ticklabel_end=cfg.r_size,
                       tickbox=True, grid=use_grid, is_invert_y_axe=False,
                       aspect='equal', image_interpolation='nearest')
    
    # add subplots
    plot.add_subplot_cartesian_2d(e_r_plot_name, 311, x_axe_label=x_axis_label, y_axe_label=y_axis_label)
    plot.add_subplot_cartesian_2d(e_z_plot_name, 312, x_axe_label=x_axis_label, y_axe_label=y_axis_label)
    plot.add_subplot_cartesian_2d(rho_beam_plot_name, 313, x_axe_label=x_axis_label, y_axe_label=y_axis_label)

    frame = cfg.get_frame_number_by_timestamp(timestamp)

    data_r = reader.get_frame('E_r', frame)
    data_z = reader.get_frame('E_z', frame)
    data_beam = reader.get_frame('rho_beam', frame)
    
    # add timestamp to each frame
    plot.get_figure().suptitle("Time: {:.2e} s".format(timestamp), x=.85, y=.95)
    
    plot.add_image(e_r_plot_name, data_r, cmap=cmap, clim=clim_e_r)
    plot.add_image(e_z_plot_name, data_z, cmap=cmap, clim=clim_e_z)
    plot.add_image(rho_beam_plot_name, data_beam, cmap=cmap, clim=clim_rho_beam)

    if not dry_run: plot.get_figure().savefig(image_file)
    if view: plot.show()
    # if not dry_run: writer.grab_frame()


def main():
    ## configure RC properties
    # plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

    ####
    parser = argparse.ArgumentParser(description='Tool for creating movie file from PDP3 modelled data.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    parser.add_argument('--image-file', type=str,
                        help='Full path to output image file. Default <path/to/parameters.xml>/image_<timestamp>.png',
                        default=None)

    parser.add_argument('--timestamp', type=float, help='Timestamp for image', default=0)

    parser.add_argument('--cmap', type=str,
                        help='''Use custom colormap. Default %s.
                        Reference: https://matplotlib.org/examples/color/colormaps_reference.html''' % 'gray',
                        default=None)

    parser.add_argument('--beam-scale-factor', type=float,
                        help='''Beam density setting automatically, but you can set scale factor to sets,
                        where initial bunch density should be placed in color range''',
                        default=1)

    parser.add_argument('--clim-e-r', type=str,
                        help='Color limit range for Electrical field longitual component', default=None)

    parser.add_argument('--clim-e-z', type=str,
                        help='Color limit range for Electrical field radial component', default=None)

    parser.add_argument('--dry-run', action='store_true', help='Do not write anything. Just for debug', default=False)

    parser.add_argument('--view', action='store_true', default=False,
                        help='View animation as well as write it to file')
    
    parser.add_argument('--with-grid', action='store_true', help='Use tick grid for plots', default=False)

    args = parser.parse_args()

    # check if config file exists
    if os.path.isfile(args.properties_path):
        clim_e_r = list(map(float, args.clim_e_r.split(':'))) if args.clim_e_r else None
        clim_e_z = list(map(float, args.clim_e_z.split(':'))) if args.clim_e_r else None
        run(args.properties_path,
            clim_e_r=clim_e_r,
            clim_e_z=clim_e_z,
            rho_beam_scale=args.beam_scale_factor,
            image_file=args.image_file,
            timestamp=args.timestamp,
            cmap=args.cmap,
            dry_run=args.dry_run,
            view=args.view,
            use_grid=args.with_grid)
        if args.view:
            input("Press 'Return' to exit ")
    else:
        print("Configuration file `%s' does not exists. Exiting" % args.properties_path)
        exit(1)

## call main function
if __name__ == "__main__":
    main()
