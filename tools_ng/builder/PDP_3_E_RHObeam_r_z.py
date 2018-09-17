from os import path

# from lib.plot_builder import PlotBuilder
import matplotlib.animation as ani
import matplotlib.pyplot as plt

from builder.builder import Builder

# README:
# color map reference: https://matplotlib.org/examples/color/colormaps_reference.html
# mathtext reference:  https://matplotlib.org/users/mathtext.html

class PDP3ERHOBeam(Builder):
    def __init__(self,
                 parameters_file, video_file,
                 cmap=None, clim_e_r=None, clim_e_z=None,
                 beam_scale_factor=0.1,
                 use_grid=False, time_range=None,
                 view=False, dry_run=False):

        self.__fig_color__=None
        self.__fig_width__=10.5
        self.__fig_height__=7
        self.__fig_dpi__=100
        self.__font_family__='sans-serif'
        self.__font_name__='DejaVu Sans'
        self.__font_size__=10

        super(PDP3ERHOBeam, self).__init__(parameters_file,
                                           self.__fig_color__,
                                           self.__fig_width__,
                                           self.__fig_height__,
                                           self.__fig_dpi__,
                                           self.__font_family__,
                                           self.__font_name__,
                                           self.__font_size__)

        # configure animation parameters
        plt.rcParams['animation.ffmpeg_path'] = '/usr/bin/ffmpeg'

        # get video file name
        config_file_name = 'field_movie.avi'
        self.__video_file__ = path.join(self.__cfg__.config_path, config_file_name) if not video_file else video_file

        # get clims
        self.__clim_e_r__ = clim_e_r or [-self.__cfg__.clim_estimation, self.__cfg__.clim_estimation]
        self.__clim_e_z__ = clim_e_z or [-self.__cfg__.clim_estimation, self.__cfg__.clim_estimation]

        self.__cmap__ = cmap or 'terrain'
        self.__beam_scale_factor__ = beam_scale_factor or 0.1
        self.__view__ = view
        self.__time_range__ = time_range or [self.__cfg__.start_time, self.__cfg__.end_time]

    def run(self):
        # some code to run
        if self.__view__:
            input("Press 'Return' to exit ")
        return True
