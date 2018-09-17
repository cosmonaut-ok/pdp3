import os
from builder.PDP_3_E_RHObeam_r_z import PDP3ERHOBeam
import argparse


def main():

    ####
    parser = argparse.ArgumentParser(description='Tool for creating movie file from PDP3 modelled data.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    parser.add_argument('--with-grid', action='store_true', help='Use tick grid for plots')

    parser.add_argument('--video-file', type=str,
                        help='Full path to output video file. Default <path/to/parameters.xml>/field_movie.avi')

    parser.add_argument('--time-range', type=str, help='Time range')

    parser.add_argument('--cmap', type=str,
                        help='''Use custom colormap. Default %s.
                        Reference: https://matplotlib.org/examples/color/colormaps_reference.html''' % 'terrain',
                        default='terrain')

    parser.add_argument('--beam-scale-factor', type=float,
                        help='''Beam density setting automatically, but you can set scale factor to sets,
                        where initial bunch density should be placed in color range''',
                        default=0.1)

    parser.add_argument('--clim-e-r', type=str,
                        help='Color limit range for Electrical field longitual component')
    parser.add_argument('--clim-e-z', type=str,
                        help='Color limit range for Electrical field radial component')

    parser.add_argument('--dry-run', action='store_true', help='Do not write anything. Just for debug')

    parser.add_argument('--view', action='store_true', default=False,
                        help='View animation as well as write it to file')

    args = parser.parse_args()

    # check if config file exists
    if os.path.isfile(args.properties_path):
        ############################################################################################
        #################### configure plot and view parameters ####################################
        ############################################################################################
        clim_e_r = list(map(float, args.clim_e_r.split(':'))) if args.clim_e_r else None
        clim_e_z = list(map(float, args.clim_e_z.split(':'))) if args.clim_e_z else None
        time_range = list(map(float, args.time_range.split(':'))) if args.time_range else None

        movie = PDP3ERHOBeam(args.properties_path,
                             video_file=args.video_file,
                             clim_e_r=clim_e_r, clim_e_z=clim_e_z,
                             cmap=args.cmap,
                             beam_scale_factor=args.beam_scale_factor,
                             use_grid=args.with_grid,
                             time_range=time_range,
                             view=args.view, dry_run=args.dry_run)

        movie.setup_figure()
        movie.run()

    else:
        print("Configuration file `%s' does not exists. Exiting" % args.properties_path)
        exit(1)


# call main function
if __name__ == "__main__":
    main()
