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

def run(config_file, clim_t, video_file=None,
        time_range=None, cmap=None, frame_step=1, specie='electrons', dry_run=False, view=False, use_grid=False):

    ##  configuration options
    x_axis_label = r'$\mathit{Z (m)}$'
    y_axis_label = r'$\mathit{R (m)}$'
    cbar_axis_label = r'$\mathit{T (eV)}$'

    t_plot_name = r'$\mathbf{Temperature\enspace}\enspace(T)$'

    ## read configfile
    cfg = Parameters(config_file)

    # calculate/update color limit values
    if not clim_t: clim_t = [0, 4]

    # calculate/update video file path
    video_file = os.path.join(os.path.dirname(config_file), 'temperature_movie.avi') if not video_file else video_file

    # define reader (plain reader used)
    autoselect = True
    use_cache = False
    if not cfg.use_hdf5:
        reader = PlainReader(path = cfg.data_path,
                             data_root=cfg.data_root,
                             fullframe_size=[cfg.number_r_grid, cfg.number_z_grid],
                             fpds=cfg.frames_per_file,
                             use_cache=use_cache,
                             verbose=False)
    else:
        reader = H5Reader(str(os.path.join(cfg.data_path, 'data.h5')), use_cache=False)
        reader.verbose = True

    # define plot builder
    plot = PlotBuilder(cfg.number_z_grid, cfg.number_r_grid,
                       fig_color=cfg.figure_color, fig_width=cfg.figure_width,
                       fig_height=cfg.figure_height, fig_dpi=cfg.figure_dpi,
                       font_family=cfg.figure_font_family,
                       font_name=cfg.figure_font_name,
                       font_size=cfg.figure_font_size,

                       x_ticklabel_end=cfg.z_size, y_ticklabel_end=cfg.r_size,
                       tickbox=True, grid=use_grid, is_invert_y_axe=False,
                       aspect='equal', image_interpolation='bicubic')

    # add subplots
    plot.add_subplot_cartesian_2d(t_plot_name, 111, x_axe_label=x_axis_label, y_axe_label=y_axis_label)

    # add initial image with zeros and colorbar
    initial_image = np.zeros([cfg.number_r_grid, cfg.number_z_grid])

    plot.add_image(t_plot_name, initial_image, cmap=cmap, clim=clim_t)
    plot.add_colorbar(t_plot_name, ticks=clim_t, title=cbar_axis_label)

    if view: plot.show()

    # dirty hack
    for p in cfg.probes:
        if p.component == 'T':
            dump_interval = p.schedule
            break

    if not time_range:
        start_frame = cfg.get_frame_number_by_timestamp(cfg.start_time, dump_interval)
        end_frame = cfg.get_frame_number_by_timestamp(cfg.end_time, dump_interval)
    else:
        if time_range[0] > time_range[1]: raise ValueError("End time should not be less, than start time. The values were: {}, {}".format(time_range[0], time_range[1]))
        if time_range[1] > cfg.end_time: raise IndexError("End time is out of simulation range {}. The value was {}".format(cfg.end_time, time_range[1]))

        start_frame = cfg.get_frame_number_by_timestamp(time_range[0], dump_interval)
        end_frame = cfg.get_frame_number_by_timestamp(time_range[1], dump_interval)
    if cfg.use_hdf5:
        start_data_set = start_frame
        end_data_set = end_frame
    else:
        start_data_set, _ = reader.get_ds_frame_by_frame(start_frame)
        end_data_set, _ = reader.get_ds_frame_by_frame(end_frame)

    FFMpegWriter = ani.writers['ffmpeg']
    metadata = dict(title='Movie Test', artist='Matplotlib',
                    comment='Movie support!')
    writer = FFMpegWriter(fps=cfg.video_fps,
                          metadata=metadata,
                          codec=cfg.video_codec,
                          bitrate=cfg.video_bitrate)

    if dry_run: video_file = '/dev/null'
    fig = plot.get_figure()

    with writer.saving(fig, video_file, cfg.figure_dpi):
        for i in range(start_data_set, end_data_set):
            sys.stdout.write('Loading dataset ' + str(i) + ' ')
            sys.stdout.flush()
            if cfg.use_hdf5:
                data_t = reader.get_frame('T', i)

                # add timestamp to each frame
                timestamp = cfg.get_timestamp_by_frame_number(i, dump_interval)
                fig.suptitle("Time: {:.2e} s".format(timestamp), x=.85, y=.95)
                if i % frame_step == 0:
                    plot.add_image(t_plot_name, data_t, cmap=cmap, clim=clim_t)

                    if view: plot.redraw()
                    if not dry_run: writer.grab_frame()
                    print('done')
                else:
                    print('skip')
            else:
                shape = [0, 0, cfg.number_r_grid, cfg.number_z_grid]
                data_t = reader.get_all_frames_in_ds('T/{}'.format(specie), shape, i)

                frame = 0
                for t in data_t:
                    if (frame + i * cfg.frames_per_file) % frame_step == 0:
                        # print without newline
                        sys.stdout.write('.')
                        sys.stdout.flush()

                        # add timestamp to each frame
                        timestamp = cfg.get_timestamp_by_frame_number(frame + i * cfg.frames_per_file, dump_interval)
                        fig.suptitle("Time: {:.2e} s".format(timestamp), x=.85, y=.95)

                        if i % frame_step == 0:
                            plot.add_image(t_plot_name, t, cmap=cmap, clim=clim_t)

                            if view: plot.redraw()
                            if not dry_run: writer.grab_frame()
                    frame = frame + 1
                print('done')


def main():
    ## configure RC properties
    # plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

    ####
    parser = argparse.ArgumentParser(description='Tool for creating movie file from PDP3 modelled data.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    parser.add_argument('--video-file', type=str,
                        help='Full path to output video file. Default <path/to/parameters.xml>/field_movie.avi',
                        default=None)

    parser.add_argument('--time-range', type=str, help='Time range', default=None)

    parser.add_argument('--cmap', type=str,
                        help='''Use custom colormap. Default %s.
                        Reference: https://matplotlib.org/examples/color/colormaps_reference.html''' % 'gray',
                        default=None)

    parser.add_argument('--beam-scale-factor', type=float,
                        help='''Beam density setting automatically, but you can set scale factor to sets,
                        where initial bunch density should be placed in color range''',
                        default=1)

    parser.add_argument('--frame-step', type=int,
                        help='Use only every Nth frame for movie creation. Default: 1',
                        default=1)

    parser.add_argument('--clim-t', type=str,
                        help='Color limit range for Temperature', default=None)

    parser.add_argument('--dry-run', action='store_true', help='Do not write anything. Just for debug', default=False)

    parser.add_argument('--view', action='store_true', default=False,
                        help='View animation as well as write it to file')

    parser.add_argument('--with-grid', action='store_true', help='Use tick grid for plots', default=False)

    args = parser.parse_args()

    time_range=None

    # check if config file exists
    if os.path.isfile(args.properties_path):
        if args.time_range:
            time_range = list(map(float, args.time_range.split(':')))

        clim_t = list(map(float, args.clim_t.split(':'))) if args.clim_t else None
        run(args.properties_path,
            clim_t=clim_t,
            video_file=args.video_file,
            time_range=time_range,
            frame_step=args.frame_step,
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

# run(config_file, clim_t, clim_e_z, clim_rho_beam, video_file='field_movie.avi',
#     time_range=None, cmap=None, dry-run=False, view=False, use_grid=False):
# input("Press 'Return' to exit ")
