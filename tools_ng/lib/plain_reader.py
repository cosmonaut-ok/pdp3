import os
import sys
import re
import fnmatch

import numpy as np
from os.path import join # to use "join" for namespaces

class PlainReader:
    '''
    ds: dataset
    space: direct dataset space (like 'E_z')
    frame: frame in single dataset (determined by fpds - frame-per-dataset)
    row: row in frame (which consists of rows and columns)
    col: column in frame
    '''
    def __init__(self, data_path, dump_path=None, shape=[0, 0], fpds=1):
        self.__data_path__ = data_path
        self.__dump_path__ = dump_path or data_path
        self.__shape__ = [shape[0], shape[1]]
        self.__fpds__ = fpds
        self.__data_set_ranges__ = {}
        self.__frame_set_ranges__ = {}

    #### special/service functions

    def __del__(self):
        return self


    def __enter__(self):
        '''
        to use as
        ```
        with H5Reader('path/to/db.h5') as h5:
            ...
        ```
        '''
        return self


    def __exit__(self, exception_type, exception_value, traceback):
        self.__del__()


    def get_ds_frame_by_frame(self, frame):
        ''' get dataset and frame in this file by absolute frame number '''
        frame_file = frame // self.__fpds__
        frame_in_last_file = frame % self.__fpds__ # - 1

        return frame_file, frame_in_last_file


    def __get_path__(self, space, ds=''):
        path = join(self.__data_path__, "{}{}".format(space, str(ds)))

        return path


    def __get_file_range__(self, space):
        files = fnmatch.filter(os.listdir(self.__data_path__), '{}*'.format(space))

        return max(map(lambda x: int(re.sub(space, '', x)), files))

    def __get_ds_range__(self, space):
        try:
            frange = self.__data_set_ranges__[space]
        except KeyError:
            frange = self.__get_file_range__(space)
            self.__data_set_ranges__[space] = frange

        return frange


    #### validators

    def __validate_ds__(self, space, number):
        frange = self.__get_ds_range__(space)
        if frange < number:
            raise IndexError('Dataset should be less, than {}. The value was {}.'
                             .format(frange, number))
        else:

            return True

    def __validate_frame__(self, space, frame):
        if frame < 0:
            raise IndexError('frame should not be less, than 0. The value was: {}'.format(frame))
        else:
            frame_length = self.__get_ds_range__(space) * self.__fpds__
            if frame_length <= frame:
                raise IndexError('frame should be less, than {}. The value was: {}'
                                .format(frame_length - 1, frame))
            else:

                return True


    def __validate_row__(self, _space, row):
        if row < 0:
            raise IndexError('row should not be less than 0. The value was {}'
                            .format(row))
        else:
            if self.__shape__[0] <= row:
                raise IndexError('Out of range: row should be less, than {}. The value was {}.'
                                .format(row, self.__shape__[0]))
            else:

                return True


    def __validate_col__(self, _space, col):
        if col < 0:
            raise IndexError('column should not be less than 0. The value was {}'
                            .format(col))
        else:
            if self.__shape__[1] <= col:
                raise IndexError('Out of range: column should be less, than {}. The value was {}.'
                                .format(col, self.__shape__[1]))
            else:

                return True


    def __validate_frame_range__(self, space, from_frame, to_frame):
        if to_frame < from_frame:
            raise IndexError('from_frame should be less or equal, than to_frame. The values was: {} and {}'
                            .format(from_frame, to_frame))
        else:

            return self.__validate_frame__(space, to_frame) and self.__validate_frame__(space, to_frame)


    def __validate_row_range__(self, space, from_row, to_row):
        if to_row < from_row:
            raise IndexError('from_row should be less than to_row. The values were {} and {}'
                            .format(from_row, to_row))
        else:

            return self.__validate_row__(space, from_row) and self.__validate_row__(space, to_row)


    def __validate_col_range__(self, space, from_col, to_col):
        if to_col < from_col:
            raise IndexError('from_col should be less than to_col. The values were {} and {}'
                            .format(from_col, to_col))
        else:

            return self.__validate_col__(space, from_col) and self.__validate_col__(space, to_col)


    #### external interface functions

    def get_all_frames_in_ds(self, space, ds_number):
        ''' get all frames in dataset.
        helper function only for plain_reader
        to optimize reading from files set for
        external user '''
        self.__validate_ds__(space, ds_number)

        path = self.__get_path__(space, ds_number)
        with open(path, 'r', encoding='utf-8') as datafile:
            frames = np.fromfile(datafile, dtype=float,
                                 count=self.__shape__[0] * self.__shape__[1] * self.__fpds__,
                                 sep=' ')

        real_shape = [self.__fpds__, self.__shape__[0], self.__shape__[1]]

        return np.reshape(frames, real_shape)


    def get_frame(self, space, number):
        ''' get frame by number. Find required dataset automatically '''
        self.__validate_frame__(space, number)

        frameds, frame_in_ds = self.get_ds_frame_by_frame(number)
        self.__validate_ds__(space, frameds)
        frame = self.get_all_frames_in_ds(space, frameds)[frame_in_ds]

        return frame


    def get_row(self, space, number, row_number):
        self.__validate_row__(space, row_number)
        row = self.get_frame(space, number)[row_number]

        return row


    def get_col(self, space, number, col_number):
        self.__validate_col__(space, col_number)
        col = self.get_frame(space, number)[:, col_number]

        return col


    def get_point(self, space, number, row_number, col_number):
        self.__validate_row__(space, row_number)
        self.__validate_col__(space, col_number)
        point = self.get_frame(space, number)[row_number, col_number]

        return point


    ###################################################################

    def get_frame_range(self, space, from_frame=0, to_frame=None):
        ''' get frame range to 3D array '''
        self.__validate_frame_range__(space, from_frame, to_frame)

        if not to_frame:
            to_frame = self.__get_ds_range__(space)

        frame_range = to_frame - from_frame
        from_ds, from_frame_in_ds = self.get_ds_frame_by_frame(from_frame)
        to_ds, to_frame_in_ds = self.get_ds_frame_by_frame(to_frame)
        frames = np.empty([frame_range, self.__shape__[0], self.__shape__[1]])

        # first ds
        frames[0:self.__fpds__ - from_frame_in_ds - 1] = self.get_all_frames_in_ds(space, from_ds)[from_frame_in_ds:self.__fpds__ - 1]
        # last ds
        frames[frame_range - to_frame_in_ds - 1:frame_range - 1] = self.get_all_frames_in_ds(space, to_ds)[:to_frame_in_ds]

        shifted_frame = self.__fpds__ - from_frame_in_ds + 1
        for i in range(from_ds + 1, to_ds):
            i_shifted = i - from_ds - 1
            k = i_shifted * self.__fpds__
            k_1 = (i_shifted + 1) * self.__fpds__
            frames[shifted_frame + k:shifted_frame + k_1 - 1] = self.get_all_frames_in_ds(space, i)[0:self.__fpds__ - 1]

        return frames


    def get_frame_range_col(self, space, number, from_frame=0, to_frame=None):
        self.__validate_frame_range__(space, from_frame, to_frame)

        if not to_frame:
            to_frame = self.__get_ds_range__(space)

        frame_range = to_frame - from_frame
        from_ds, from_frame_in_ds = self.get_ds_frame_by_frame(from_frame)
        to_ds, to_frame_in_ds = self.get_ds_frame_by_frame(to_frame)

        frames = np.empty([frame_range, self.__shape__[0]])

        # first ds
        frames[0:self.__fpds__ - from_frame_in_ds - 1] = self.get_all_frames_in_ds(space, from_ds)[from_frame_in_ds:self.__fpds__ - 1, :, number]
        # last ds
        frames[frame_range - to_frame_in_ds - 1:frame_range - 1] = self.get_all_frames_in_ds(space, to_ds)[:to_frame_in_ds, :, number]

        shifted_frame = self.__fpds__ - from_frame_in_ds + 1
        for i in range(from_ds + 1, to_ds):
            i_shifted = i - from_ds - 1
            k = i_shifted * self.__fpds__
            k_1 = (i_shifted + 1) * self.__fpds__
            frames[shifted_frame + k:shifted_frame + k_1 - 1] = self.get_all_frames_in_ds(space, i)[0:self.__fpds__ - 1, :, number]

        return frames


    def get_frame_range_row(self, space, number, from_frame=0, to_frame=None):
        self.__validate_frame_range__(space, from_frame, to_frame)

        if not to_frame:
            to_frame = self.__get_ds_range__(space)

        frame_range = to_frame - from_frame
        from_ds, from_frame_in_ds = self.get_ds_frame_by_frame(from_frame)
        to_ds, to_frame_in_ds = self.get_ds_frame_by_frame(to_frame)

        frames = np.empty([frame_range, self.__shape__[1]])

        # first ds
        frames[0:self.__fpds__ - from_frame_in_ds - 1] = self.get_all_frames_in_ds(space, from_ds)[from_frame_in_ds:self.__fpds__ - 1, number]
        # last ds
        frames[frame_range - to_frame_in_ds - 1:frame_range - 1] = self.get_all_frames_in_ds(space, to_ds)[:to_frame_in_ds, number]

        shifted_frame = self.__fpds__ - from_frame_in_ds + 1
        for i in range(from_ds + 1, to_ds):
            i_shifted = i - from_ds - 1
            k = i_shifted * self.__fpds__
            k_1 = (i_shifted + 1) * self.__fpds__
            frames[shifted_frame + k:shifted_frame + k_1 - 1] = self.get_all_frames_in_ds(space, i)[0:self.__fpds__ - 1, number]

        return frames


    def get_frame_range_dot(self, space, row_number, col_number, from_frame=0, to_frame=None):
        self.__validate_frame_range__(space, from_frame, to_frame)

        if not to_frame:
            to_frame = self.__get_ds_range__(space)

        frame_range = to_frame - from_frame
        from_ds, from_frame_in_ds = self.get_ds_frame_by_frame(from_frame)
        to_ds, to_frame_in_ds = self.get_ds_frame_by_frame(to_frame)

        frames = np.empty(frame_range)

        # first ds
        frames[0:self.__fpds__ - from_frame_in_ds - 1] = self.get_all_frames_in_ds(space, from_ds)[from_frame_in_ds:self.__fpds__ - 1, row_number, col_number]
        # last ds
        frames[frame_range - to_frame_in_ds - 1:frame_range - 1] = self.get_all_frames_in_ds(space, to_ds)[:to_frame_in_ds, row_number, col_number]

        shifted_frame = self.__fpds__ - from_frame_in_ds + 1
        for i in range(from_ds + 1, to_ds):
            i_shifted = i - from_ds - 1
            k = i_shifted * self.__fpds__
            k_1 = (i_shifted + 1) * self.__fpds__
            frames[shifted_frame + k:shifted_frame + k_1 - 1] = self.get_all_frames_in_ds(space, i)[0:self.__fpds__ - 1, row_number, col_number]

        return frames
