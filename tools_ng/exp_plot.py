from lib.parameters import Parameters
from lib.plot_builder import PlotBuilder
from lib.h5_reader import H5Reader
from lib.plain_reader import PlainReader

import sys

##  configuration options
#  set color limit for E_r, E_z and rho_beam
clim_e_r = clim_e_z = [-5e5, 5e5]

#  set timestamp
timestamp = 1e-9

# set color map
cmap = 'terrain'

# define configuration file
# config_file = '/home/cosmonaut/pdp3_modeling/model32_single/parameters.xml'
config_file = '../testdir/parameters.xml'

## read configfile
cfg = Parameters(config_file)

clim_rho_beam = [-(cfg.bunch_density * 1.6e-19), 0]

# define reader (plain reader used)
reader = PlainReader(cfg.data_path, cfg.system_state_path,
                     [cfg.number_r_grid , cfg.number_z_grid],
                     cfg.frames_per_file, False)

# define plot builder
plot = PlotBuilder(cfg.number_z_grid, cfg.number_r_grid,
                   fig_color=cfg.figure_color, fig_width=cfg.figure_width,
                   fig_height=cfg.figure_height, fig_dpi=cfg.figure_dpi,
                   font_family=cfg.figure_font_family,
                   font_name=cfg.figure_font_name,
                   font_size=cfg.figure_font_size,

                   x_ticklabel_end=cfg.z_size, y_ticklabel_end=cfg.r_size,
                   tickbox=True, grid=True, is_invert_y_axe=False,
                   aspect='equal', image_interpolation='nearest')

# add subplots
plot.add_subplot_cartesian_2d(r'$E_r$', 311)

plot.add_subplot_cartesian_2d(r'$E_z$', 312)
plot.add_subplot_cartesian_2d(r'$\rho_{beam}$', 313)
plot.show()

start_frame = cfg.get_frame_number_by_timestamp(cfg.start_time)
end_frame = cfg.get_frame_number_by_timestamp(cfg.end_time)

start_data_set, _ = reader.get_ds_frame_by_frame(start_frame)
end_data_set, _ = reader.get_ds_frame_by_frame(end_frame)

initial_image = zeros([cfg.number_r_grid, cfg.number_r_grid])

plot.add_image(r'$E_r$', initial_image, cmap=cmap, clim=clim_e_r)
plot.add_colorbar(r'$E_r$', ticks=clim_e_r)

plot.add_image(r'$E_z$', initial_image, cmap=cmap, clim=clim_e_z)
plot.add_colorbar(r'$E_z$', ticks=clim_e_z)

plot.add_image(r'$\rho_{beam}$', initial_image, cmap=cmap, clim=clim_rho_beam)
plot.add_colorbar(r'$\rho_{beam}$', ticks=clim_rho_beam)

for i in range(start_data_set, end_data_set):
    sys.stdout.write('Loading dataset ' + str(i) + ' ')
    sys.stdout.flush()

    data_r = reader.get_all_frames_in_ds('E_r', i)
    data_z = reader.get_all_frames_in_ds('E_z', i)
    data_beam = reader.get_all_frames_in_ds('rho_beam', i)

    for r, z, beam in zip(data_r, data_z, data_beam):
        # print without newline
        sys.stdout.write('.')
        sys.stdout.flush()

        plot.add_image(r'$E_r$', r, cmap=cmap, clim=clim_e_r)
        plot.add_image(r'$E_z$', z, cmap=cmap, clim=clim_e_z)
        plot.add_image(r'$\rho_{beam}$', beam, cmap=cmap, clim=clim_rho_beam)
        plot.redraw()
    print()

input("Press 'Return' to exit ")
