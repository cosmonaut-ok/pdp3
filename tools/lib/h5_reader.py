import numpy as np
import h5py
from os.path import join # to use "join" for namespaces

class H5Reader:
    def __init__(self, h5_path, data_keyspace='/pdp3/result', dump_keyspace='/pdp3/dump', use_cache=False):

        self.file = h5py.File(h5_path)
        self.__data_keyspace__ = data_keyspace
        self.__dump_keyspace__ = dump_keyspace
        self.verbose = False


    def __del__(self):
        self.file.close()


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


    def __get_path__(self, space, ds=''):
        path = join(self.__data_keyspace__, space, str(ds))
        return path


    def __check_frame__(self, space, frame):
        if frame < 0:
            raise Exception('frame should not be less, than 0. The value was: {}'.format(frame))
        else:
            path = self.__get_path__(space)
            space_length = len(self.file[path])
            if space_length < frame:
                raise Exception('frame should be less, than {}. The value was: {}'.format(space_length, frame))
            else:
                return True


    def __check_row__(self, space, row):
        if row < 0:
            raise Exception('row should not be less than 0. The value was {}'.format(row))
        else:
            path = self.__get_path__(space, 0)
            frame_length = len(self.file[path])
            if frame_length < row:
                raise Exception('Out of range: row should be less, than {}. The value was {}.'.format(row, frame_length))
            else:
                return True


    def __check_col__(self, space, col):
        if col < 0:
            raise Exception('column should not be less than 0. The value was {}'.format(col))
        else:
            path = self.__get_path__(space, 0)
            frame_height = len(self.file[path][0])
            if frame_height < col:
                raise Exception('Out of range: column should be less, than {}. The value was {}.'.format(col, frame_height))
            else:
                return True


    def __check_frame_range__(self, space, from_frame, to_frame):
        if to_frame < from_frame:
            raise Exception('from_frame should be less or equal, than to_frame. The values was: {} and {}'.format(from_frame, to_frame))
        elif from_frame < 0:
            raise Exception('from_frame should not be less, than 0. The value was: {}'.format(from_frame))
        else:
            return self.__check_frame__(space, to_frame)


    def __check_row_range__(self, space, from_row, to_row):
        if from_row < 0:
            raise Exception('from_row should not be less than 0. The value was {}'.format(from_row))
        elif to_row < from_row:
            raise Exception('from_row should be less than to_row. The values were {} and {}'.format(from_row, to_row))
        else:
            return self.__check_row__(space, to_row)


    def __check_col_range__(self, space, from_col, to_col):
        if from_col < 0:
            raise Exception('from_col should not be less than 0. The value was {}'.format(from_col))
        elif to_col < from_col:
            raise Exception('from_col should be less than to_col. The values were {} and {}'.format(from_col, to_col))
        else:
            return self.__check_col__(space, to_col)

#################################################################################################
#################################################################################################
#################################################################################################

    def get_frame(self, space, number):
        if self.__check_frame__(space, number):
            path = self.__get_path__(space, number)
            frame = self.file[path][:]
            return frame


    def get_row(self, space, number, row_number):
        if self.__check_frame__(space, number) and self.__check_row__(space, row_number):
            path = self.__get_path__(space, number)
            row = self.file[path][row_number]
            return row


    def get_col(self, space, number, col_number):
        if self.__check_frame__(space, number) and self.__check_col__(space, col_number):
            path = self.__get_path__(space, number)
            col = self.file[path][:,col_number]
            return col


    def get_point(self, space, number, row_number, col_number):
        if self.__check_frame__(space, number) and self.__check_row__(space, row_number) and self.__check_col__(space, col_number):
            path = self.__get_path__(space, number)
            point = self.file[path][row_number][col_number]
            return point


###################################################################

    def get_frame_range(self, space, from_frame=0, to_frame=None):
        path = self.__get_path__(space, 0)
        frame_length = len(self.file[path])
        frame_height = len(self.file[path][0])

        if not to_frame:
            path = self.__get_path__(space)
            space_length = len(self.file[path])
            to_frame = space_length-1

        if self.__check_frame_range__(space, from_frame, to_frame):
            frames = np.empty([to_frame - from_frame + 1, frame_length, frame_height])

            for i in range(from_frame, to_frame + 1):
                path = self.__get_path__(space, i)
                frames[i-from_frame] = self.file[path]
            return frames


    def get_frame_range_row(self, space, row_number=0, from_frame=0, to_frame=None):
        path = self.__get_path__(space, 0)
        frame_length = len(self.file[path])
        frame_height = len(self.file[path][0])

        if not to_frame:
            space_length = len(self.file[path])
            to_frame = space_length-1

        if self.__check_frame_range__(space, from_frame, to_frame) and self.__check_row__(space, row_number):
            rows = []
            for i in range(from_frame, to_frame+1):
                path = self.__get_path__(space, i)
                rows.append(self.file[path][row_number])
            return rows


    def get_frame_range_col(self, space, col_number=0, from_frame=0, to_frame=None):
        if not to_frame:
            path = self.__get_path__(space, 0)
            space_length = len(self.file[path])
            to_frame = space_length-1

        if self.__check_frame_range__(space, from_frame, to_frame) and self.__check_col__(space, col_number):
            cols = []
            for i in range(from_frame, to_frame+1):
                path = self.__get_path__(space, i)
                cols.append(self.file[path][:,col_number])
            return cols


    def get_frame_range_dot(self, space, row_number=0, col_number=0, from_frame=0, to_frame=None):
        if not to_frame:
            path = self.__get_path__(space, 0)
            space_length = len(self.file[path])
            to_frame = space_length-1

        if self.__check_frame_range__(space, from_frame, to_frame) and self.__check_row__(space, row_number) and self.__check_col__(space, col_number):
            dots = []
            for i in range(from_frame, to_frame+1):
                path = self.__get_path__(space, i)
                dots.append(self.file[path][row_number,col_number])
            return dots
