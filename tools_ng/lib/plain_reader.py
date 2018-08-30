import os
import sys
import re
import fnmatch

import numpy as np
from os.path import join # to use "join" for namespaces

class PlainReader:
    def __init__(self, data_path, dump_path=None, shape=[0, 0], fpf=1):
        self.__data_path__ = data_path
        self.__dump_path__ = dump_path or data_path
        self.__shape__ = shape
        self.__fpf__ = fpf


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


    def __exit__(self):
        self.__del__()

    def __get_file_frame_by_frame__(self, frame):
        ''' get file and frame in this file by absolute frame number '''
        frame_file = frame // fpf
        frame_in_last_file = frame % fpf
        return(frame_file, frame_in_last_file)

    def __get_path__(self, space, ds=''):
        path = join(self.__data_path__, "{}_{}".format(path, str(ds)))
        return(path)


    def __get_file_range__(self, space):
        files = fnmatch.filter(os.listdir(self.__data_path__), '{}*'.format(space))
        return(max(map(lambda x: int(re.sub(space, '', x)), files)))


    def __check_frame__(self, space, frame):
        if frame < 0:
            raise Exception('frame should not be less, than 0. The value was: {}'.format(frame))
        else:
            path = self.__get_path__(space)
            space_length = self.__get_file_range__(space) * self.__fpf__
            if space_length < frame:
                raise Exception('frame should be less, than {}. The value was: {}'.format(space_length, frame))
            else:
                return(True)


    def __check_row__(self, space, row):
        if row < 0:
            raise Exception('row should not be less than 0. The value was {}'.format(row))
        else:
            if self.__shape__[0] < row:
                raise Exception('Out of range: row should be less, than {}. The value was {}.'.format(row, self.__shape__[0]))
            else:
                return(True)


    def __check_col__(self, space, col):
        if col < 0:
            raise Exception('column should not be less than 0. The value was {}'.format(col))
        else:
            if self.__shape__[1] < col:
                raise Exception('Out of range: column should be less, than {}. The value was {}.'.format(col, self.__shape__[1]))
            else:
                return(True)


    def __check_frame_range__(self, space, from_frame, to_frame):
        if to_frame < from_frame:
            raise Exception('from_frame should be less or equal, than to_frame. The values was: {} and {}'.format(from_frame, to_frame))
        elif from_frame < 0:
            raise Exception('from_frame should not be less, than 0. The value was: {}'.format(from_frame))
        else:
            return(self.__check_frame__(space, to_frame))


    def __check_row_range__(self, space, from_row, to_row):
        if from_row < 0:
            raise Exception('from_row should not be less than 0. The value was {}'.format(from_row))
        elif to_row < from_row:
            raise Exception('from_row should be less than to_row. The values were {} and {}'.format(from_row, to_row))
        else:
            return(self.__check_row__(space, to_row))


    def __check_col_range__(self, space, from_col, to_col):
        if from_col < 0:
            raise Exception('from_col should not be less than 0. The value was {}'.format(from_col))
        elif to_col < from_col:
            raise Exception('from_col should be less than to_col. The values were {} and {}'.format(from_col, to_col))
        else:
            return(self.__check_col__(space, to_col))

    #################################################################################################
    #################################################################################################
    #################################################################################################


    def get_all_frames_in_file(self, space, filenumber):
        if self.__check_frame__(space, number):

            path = self.__get_path__(space, framefile)
            with open(path, 'r', encoding='utf-8') as datafile:
                frames = np.fromfile(datafile, dtype=float, count=self.__shape__[0] * self.__shape__[1] * self.__fpf__, sep=' ')

            real_shape = []
            real_shape.append(self.__fpf__)
            real_shape.extend(self.__shape__)

            return(np.reshape(frames, real_shape))


    def get_frame(self, space, number):
        if self.__check_frame__(space, number):
            framefile, frame_in_file = self.__get_file_frame_by_frame__(number)
            frame = self.get_all_frames_in_file(space, framefile)[frame_in_file]

            return frame


    def get_row(self, space, number, row_number):
        if self.__check_frame__(space, number) and self.__check_row__(space, row_number):
            row = self.get_frame(space, number)[row_number]
            return(row)


    def get_col(self, space, number, col_number):
        if self.__check_frame__(space, number) and self.__check_col__(space, col_number):
            col = self.get_frame(space, number)[:,col_number]
            return(col)


    def get_point(self, space, number, row_number, col_number):
        if self.__check_frame__(space, number) and self.__check_row__(space, row_number) and self.__check_col__(space, col_number):
            point = self.get_frame(space, number)[row_number,col_number]
            return(point)

        ## ^^ meaning DONE ^^

###################################################################

    def get_frame_range(self, space, from_frame=0, to_frame=None):

        # path = self.__get_path__(space, 0)
        # frame_length = len(self.file[path][:])
        # frame_height = len(self.file[path][0])

        # if not to_frame:
        #     path = self.__get_path__(space)
        #     space_length = len(self.file[path])
        #     to_frame = space_length-1

        # if self.__check_frame_range__(space, from_frame, to_frame):
        #     frames = np.empty([to_frame - from_frame + 1, frame_length, frame_height])

        #     for i in range(from_frame, to_frame + 1):
        #         path = self.__get_path__(space, i)
        #         frames[i-from_frame] = self.file[path][:]
        #     return(frames)


    def get_frame_range_row(self, space, from_frame=0, to_frame=None, row_number=0):
        path = self.__get_path__(space, 0)
        frame_length = len(self.file[path][:])
        frame_height = len(self.file[path][0])

        if not to_frame:
            space_length = len(self.file[path])
            to_frame = space_length-1

        if self.__check_frame_range__(space, from_frame, to_frame) and self.__check_row__(space, row_number):
            rows = []
            for i in range(from_frame, to_frame+1):
                path = self.__get_path__(space, i)
                rows.append(self.file[path][row_number])
            return(rows)


    def get_frame_range_col(self, space, from_frame=0, to_frame=None, col_number=0):
        if not to_frame:
            path = self.__get_path__(space, 0)
            space_length = len(self.file[path])
            to_frame = space_length-1

        if self.__check_frame_range__(space, from_frame, to_frame) and self.__check_col__(space, col_number):
            cols = []
            for i in range(from_frame, to_frame+1):
                path = self.__get_path__(space, i)
                cols.append(self.file[path][:,col_number])
            return(cols)
