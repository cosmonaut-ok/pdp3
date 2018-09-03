import os
from builder.PDP_3_E_RHObeam_r_z import PDP3ERHOBeam
import argparse


def main():

    ####
    parser = argparse.ArgumentParser(description='Tool for creating movie file from PDP3 modelled data.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    default_data_set_range = [0, 10000]

    parser.add_argument('--with-grid', action='store_true', help='Use tick grid for plots')

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
                        help='Color limit range for Electrical field longitual component')
    parser.add_argument('--clim-e-z', type=str,
                        help='Color limit range for Electrical field radial component')

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

    # check if config file exists
    if os.path.isfile(args.properties_path):
        ################################################################################################
        #################### configure plot and view parameters #######################################
        ################################################################################################
#        movie.clim_e_field_r = list(map(float, args.clim_e_r.split(':'))) if args.clim_e_r else [-config.clim_estimation, config.clim_estimation]
#        movie.clim_e_field_z = list(map(float, args.clim_e_z.split(':'))) if args.clim_e_z else [-config.clim_estimation, config.clim_estimation]
#        movie.cmap = args.cmap
#        movie.clim_e_field_beam_scale_factor = args.beam_scale_factor
#        movie.use_grid = args.with_grid
#
#        if args.time_range:
#            time_range = list(map(float, args.time_range.split(':')))
#            movie.start_data_set, movie.start_frame = config.get_file_frame_number_by_timestamp(time_range[0])
#            movie.end_data_set, movie.end_frame = config.get_file_frame_number_by_timestamp(time_range[1])
#        elif args.data_set_range:
#            data_set_range = list(map(int, args.data_set_range.split(':')))
#            movie.start_data_set = data_set_range[0]
#            movie.end_data_set = data_set_range[1]
#        else:
#            movie.start_data_set = default_data_set_range[0]
#            movie.end_data_set = default_data_set_range[1]
        ################################################################################################
        ################################################################################################
        ################################################################################################

        ## initialize config
        movie = PDP3ERHOBeam(args.properties_path, args.video_file)

#        movie.setup_3e_view(view)
#        movie.create_view_with_3_plots(view, write)
        if view:
            input("Press 'Return' to exit ")
    else:
        print("Configuration file `%s' does not exists. Exiting" % args.properties_path)
        exit(1)



# call main function
if __name__ == "__main__":
    main()
