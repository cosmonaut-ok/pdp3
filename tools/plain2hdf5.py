#!/usr/bin/env python3

import os
import argparse
import numpy as np
import h5py

from lib.parameters import Parameters
from lib.h5_reader import H5Reader
from lib.plain_reader import PlainReader

def run(config_file, data_file=None, components=['E_r', 'E_z', 'H_phi', 'rho_beam'],
        time_range=None,
        main_group_name='/pdp3/result',
        enable_compression=False, compression_level=6):

    cfg = Parameters(config_file)

    # define reader (plain reader used)
    reader = PlainReader(cfg.data_path, cfg.system_state_path,
                         [cfg.number_r_grid , cfg.number_z_grid],
                         cfg.frames_per_file, False)
    reader.verbose = False

    file_path = data_file or os.path.join(cfg.data_path, "data.h5")

    f = h5py.File(file_path, 'w')

    f.create_group(main_group_name)

    for compon in components:
        f.create_group(os.path.join(main_group_name, compon))

    if not time_range:
        start_frame = cfg.get_frame_number_by_timestamp(cfg.start_time)
        end_frame = cfg.get_frame_number_by_timestamp(cfg.end_time)
    else:
        if time_range[0] > time_range[1]: raise ValueError("End time should not be less, than start time. The values were: {}, {}".format(time_range[0], time_range[1]))
        if time_range[1] > cfg.end_time: raise IndexError("End time is out of simulation range {}. The value was {}".format(cfg.end_time, time_range[1]))
        start_frame = cfg.get_frame_number_by_timestamp(time_range[0])
        end_frame = cfg.get_frame_number_by_timestamp(time_range[1])

    start_data_set, _ = reader.get_ds_frame_by_frame(start_frame)
    end_data_set, _ = reader.get_ds_frame_by_frame(end_frame)

    for i in range(start_data_set, end_data_set):
        for compon in components:
            data_compon = reader.get_all_frames_in_ds(compon, i)
            for frame in range(0, cfg.frames_per_file):
                abs_frame_number = i * cfg.frames_per_file + frame
                name = os.path.join(main_group_name, compon, str(abs_frame_number))
                ds = data_compon[frame]
                if enable_compression:
                    f.create_dataset(name, data=ds, compression='gzip', compression_opts=compression_level)
                else:
                    f.create_dataset(name, data=ds)
                print("Imported {}/{}".format(compon,  abs_frame_number))
    f.close()
    print('done')


def main():
    parser = argparse.ArgumentParser(description='Tool for creating movie file from PDP3 modelled data.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    parser.add_argument('--data-file', type=str,
                        help='Full path to output video file. Default <path/to/parameters.xml>/data.h5',
                        default=None)

    parser.add_argument('--components', type=str,
                        help='components separated by ",". E.g. --components=E_r,E_z,J_phi',
                        default='E_r,E_z,H_phi,rho_beam')

    parser.add_argument('--main-group', type=str,
                        help='Name of main HDF5 group. Default /pdp3/result',
                        default='/pdp3/result')

    parser.add_argument('--time-range', type=str, help='Time range', default=None)

    parser.add_argument('--compress', action='store_true', default=False,
                        help='Should HDF5 compression be enabled. Default false')

    parser.add_argument('--compression-level', type=float,
                        help='Set compression level (0 to 9). Default 6',
                        default=6)

    # parser.add_argument('--dry-run', action='store_true', help='Do not write anything. Just for debug', default=False)

    args = parser.parse_args()

    time_range=None

    # check if config file exists
    if os.path.isfile(args.properties_path):
        if args.time_range:
            time_range = list(map(float, args.time_range.split(':')))
        components = list(args.components.split(',')) if args.components else None

        run(args.properties_path,
            data_file=args.data_file,
            components=components,
            time_range=time_range,
            main_group_name=args.main_group,
            enable_compression=args.compress,
            compression_level=args.compression_level)

        # run(args.properties_path,
        #     clim_e_r=clim_e_r,
        #     clim_e_z=clim_e_z,
        #     rho_beam_scale=args.beam_scale_factor,
        #     video_file=args.video_file,
        #     time_range=time_range,
        #     cmap=args.cmap,
        #     dry_run=args.dry_run,
        #     view=args.view,
        #     use_grid=args.with_grid)

    else:
        print("Configuration file `%s' does not exists. Exiting" % args.properties_path)
        exit(1)

## call main function
if __name__ == "__main__":
    main()
