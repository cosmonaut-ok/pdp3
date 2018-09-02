from lib.parameters import Parameters
from lib.h5_reader import H5Reader
from lib.plain_reader import PlainReader
from lib.plot_builder import PlotBuilder

class Builder:
    def __init__(self,
                 parameters_file,
                 fig_color=None,
                 fig_width=10.5,
                 fig_height=7,
                 fig_dpi=100,
                 font_family='sans-serif',
                 font_name='DejaVu Sans',
                 font_size=10):

        self.__cfg__ = Parameters(parameters_file)

        if self.use_hdf5:
            self.__reader__ = H5Reader(self.__cfg__.data_path,
                                       base_keyspace='/pdp3/result',
                                       dump_keyspace='/pdp3/dump')
        else:
            self.__reader__ = PlainReader(self.__cfg__.data_path,
                                          self.__cfg__.dump_path,
                                          shape=[self.__cfg__.number_r_grid,
                                                 self.__cfg__.number_z_grid],
                                          fpds=self.__cfg__.frames_per_file,
                                          use_cache=False):

        self.__plotter__ = PlotBuilder(self.__cfg__,
                                       fig_color=fig_color,
                                       fig_width=fig_width,
                                       fig_heigt=fig_height,
                                       fig_dpi=fig_dpi,
                                       font_family=font_family,
                                       font_name=font_name,
                                       font_sizse=font_size)
