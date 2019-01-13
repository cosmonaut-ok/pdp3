#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# define matplotlibplotting backend
# %matplotlib -l shows all available backends


# In[ ]:


import os
import numpy as np

from lib.parameters import Parameters
from lib.plot_builder import PlotBuilder
from lib.h5_reader import H5Reader
from lib.plain_reader import PlainReader


# In[ ]:


##  configuration options
el_charge = 1.6e-19
rho_beam_scale = 1
config_file = '../parameters.xml'

cfg = Parameters(config_file)
clim_estimation = cfg.get_clim_estimation()

# update shape to get custom shaped images
# useful to get just part of frame
# or get frames, which has smaller shape than full frame
# shape=[0, 0, cfg.number_r_grid, cfg.number_z_grid]
shape=[0,0,180,1800]

timestamp=1e-9
show_grid=False
use_cache=True
cmap='terrain'

# color limits (WARNING: clim_estimation may works incorrectly)
clim_e_r = [-clim_estimation, clim_estimation]
clim_e_z = [-clim_estimation, clim_estimation]
clim_rho_beam = [-(cfg.bunch_density * el_charge * rho_beam_scale), 0]

autoselect = True

x_axis_label = r'$\mathit{Z (m)}$'
y_axis_label = r'$\mathit{R (m)}$'
cbar_axis_label = r'$\frac{V}{m}$'
cbar_bunch_density_axis_label = r'$m^{-3}$'

e_r_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Radial\enspace Component}\enspace(E_r)$'
e_z_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Longitudal\enspace Component}\enspace(E_z)$'
rho_beam_plot_name = r'$\mathbf{Electron\enspace Bunch\enspace Density}\enspace (\rho_{bunch})$'


# In[ ]:


# define reader (plain reader used)
if not cfg.use_hdf5:
    reader = PlainReader(path = cfg.data_path,
                         data_root=cfg.data_root,
                         fullframe_size=[cfg.number_r_grid , cfg.number_z_grid],
                         fpds=cfg.frames_per_file, 
                         use_cache=use_cache,
                         verbose=True)
else:
    reader = H5Reader(str(os.path.join(cfg.data_path, 'data.h5')), use_cache=use_cache)
    reader.verbose = True

r_scale = (shape[2] - shape[0]) / cfg.number_r_grid
z_scale = (shape[3] - shape[1]) / cfg.number_z_grid

# define plot builder
plot = PlotBuilder(shape[3] - shape[1], shape[2] - shape[0],
                   fig_color=cfg.figure_color, 
                   fig_width=cfg.figure_width,
                   fig_height=cfg.figure_height, 
                   fig_dpi=cfg.figure_dpi,
                   font_family=cfg.figure_font_family,
                   font_name=cfg.figure_font_name,
                   font_size=cfg.figure_font_size,
                   
                   x_ticklabel_end=cfg.z_size * z_scale, y_ticklabel_end=cfg.r_size * r_scale,
                   tickbox=True, grid=show_grid, is_invert_y_axe=False,
                   aspect='equal', image_interpolation='nearest', guess_number_ticks=20)


# In[ ]:


# add subplots
plot.add_subplot_cartesian_2d(e_r_plot_name, 311, x_axe_label=x_axis_label, y_axe_label=y_axis_label)
plot.add_subplot_cartesian_2d(e_z_plot_name, 312, x_axe_label=x_axis_label, y_axe_label=y_axis_label)
plot.add_subplot_cartesian_2d(rho_beam_plot_name, 313, x_axe_label=x_axis_label, y_axe_label=y_axis_label)

# add initial image with zeros and colorbar
initial_image = np.zeros([shape[2] - shape[0], shape[3] - shape[1]])

# add dummy images
plot.add_image(e_r_plot_name, initial_image, cmap=cmap, clim=clim_e_r)
plot.add_image(e_z_plot_name, initial_image, cmap=cmap, clim=clim_e_z)
plot.add_image(rho_beam_plot_name, initial_image, cmap=cmap, clim=clim_rho_beam)

# add colorbars
plot.add_colorbar(e_r_plot_name, ticks=clim_e_r, title=cbar_axis_label)
plot.add_colorbar(e_z_plot_name, ticks=clim_e_z, title=cbar_axis_label)
plot.add_colorbar(rho_beam_plot_name, ticks=clim_rho_beam, title=cbar_bunch_density_axis_label)


plot.show()


# In[ ]:


# get data
data_r = data_z = data_beam = []

for probe in cfg.probes:
    frame = cfg.get_frame_number_by_timestamp(timestamp, probe.schedule)
    if (probe.type == 'frame') and (probe.r_start == shape[0]) and (probe.z_start == shape[1]) and(probe.r_end == shape[2]) and(probe.z_end == shape[3]):
        if probe.component == 'E_r': data_r = reader.get_frame('E_r', shape, frame)
        if probe.component == 'E_z': data_z = reader.get_frame('E_z', shape, frame)
        if probe.component == 'rho_beam': data_beam = reader.get_frame('rho_beam', shape, frame)

# try bigger frames, if autoselect enabled
if len(data_r) == 0 and autoselect:
    for probe in cfg.probes:
        probe_shape = [probe.r_start, probe.z_start, probe.r_end, probe.z_end]
        frame = cfg.get_frame_number_by_timestamp(timestamp, probe.schedule)
        if (probe.type == 'frame') and (probe_shape[0] <= shape[0]) and (probe_shape[1] <= shape[1]) and(probe_shape[2] >= shape[2]) and(probe_shape[3] >= shape[3]):
            if probe.component == 'E_r' and len(data_r) == 0: data_r = reader.get_frame('E_r', probe_shape, frame)[shape[0]:shape[2], shape[1]:shape[3]]
            if probe.component == 'E_z' and len(data_z) == 0: data_z = reader.get_frame('E_z', probe_shape, frame)[shape[0]:shape[2], shape[1]:shape[3]]
            if probe.component == 'rho_beam' and len(data_beam) == 0: data_beam = reader.get_frame('rho_beam', probe_shape, frame)[shape[0]:shape[2], shape[1]:shape[3]]


# In[ ]:


# add timestamp to plot
plot.get_figure().suptitle("Time: {:.2e} s".format(timestamp), x=.85, y=.95)

# add images
plot.add_image(e_r_plot_name, data_r, cmap=cmap, clim=clim_e_r)
plot.add_image(e_z_plot_name, data_z, cmap=cmap, clim=clim_e_z)
plot.add_image(rho_beam_plot_name, data_beam, cmap=cmap, clim=clim_rho_beam)

plot.redraw()


input("Please press RETURN to exit ")
