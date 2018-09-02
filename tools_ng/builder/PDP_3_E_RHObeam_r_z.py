# from lib.plot_builder import PlotBuilder
from builder.builder import Builder

class PDP3ERHOBeam(Builder):
    def __init__(self,
                 parameters_file, video_file,
                 cmap='terrain', clim_e_r=[0, 1], clim_e_z=[0, 1],
                 beam_scale_factor=0.1,
                 view=True, dry_run=False):

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


# (self, parameters_xml, video_file,
#                  time_range=None,
#                  cmap='terrain', clim_e_r=, clim_e_z, beam_scale_factor, view, dry_run):
