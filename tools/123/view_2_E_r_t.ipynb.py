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
config_file = '../workdir/parameters.xml'

cfg = Parameters(config_file)
clim_estimation = cfg.get_clim_estimation()

longitude=0.025
time_range=[cfg.start_time, cfg.end_time]
use_grid=True
cmap='terrain'
clim_e_r = [-clim_estimation, clim_estimation]
clim_e_z = [-clim_estimation, clim_estimation]
use_cache = False
verbose = True
autoselect = True

x_axis_label = r'$\mathit{R (m)}$'
y_axis_label = r'$\mathit{t (ns)}$'
cbar_axis_label = r'$\frac{V}{m}$'

e_r_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Radial\enspace Component}\enspace(E_r)$'
e_z_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Longitudal\enspace Component}\enspace(E_z)$'

cbar_axis_label = r'$\frac{V}{m}$'


# In[ ]:


# define reader (plain reader used)
if not cfg.use_hdf5:
    reader = PlainReader(path = cfg.data_path,
                         data_root=cfg.data_root,
                         fullframe_size=[cfg.number_r_grid, cfg.number_z_grid],
                         fpds=cfg.frames_per_file, 
                         use_cache=use_cache,
                         verbose=verbose)
else:
    reader = H5Reader(str(os.path.join(cfg.data_path, 'data.h5')),
                      shape=[cfg.number_r_grid , cfg.number_z_grid],
                      use_cache=use_cache)
    reader.verbose = True


# In[ ]:


# get data
start_frame = cfg.get_frame_number_by_timestamp(time_range[0])
end_frame = 460 # cfg.get_frame_number_by_timestamp(time_range[1])
col_number = cfg.get_col_by_longitude(longitude)

data_r = data_z = []

for probe in cfg.probes:
    if (probe.type == 'col') and (probe.z_start == col_number):
        if probe.component == 'E_r': data_r = reader.get_frame_range_col('E_r', col_number, start_frame, end_frame)
        if probe.component == 'E_z': data_z = reader.get_frame_range_col('E_z', col_number, start_frame, end_frame)
        
if (len(data_r) == 0 or len(data_z) == 0) and autoselect:
    for probe in cfg.probes:
        if (probe.type == 'frame') and (probe.z_start <= col_number) and (probe.r_start == 0) and (probe.z_end >= col_number) and (probe.r_end == cfg.number_r_grid):
            shape = [probe.r_start, probe.z_start, probe.r_end, probe.z_end]
            if probe.component == 'E_r': data_r = reader.get_frame_range('E_r', shape, start_frame, end_frame)[col_number]
            if probe.component == 'E_z': data_z = reader.get_frame_range('E_z', shape, start_frame, end_frame)[col_number]


# In[ ]:


# define plot builder
plot = PlotBuilder(cfg.number_r_grid, len(data_r),
                   fig_color=cfg.figure_color, 
                   fig_width=cfg.figure_width,
                   fig_height=cfg.figure_height, 
                   fig_dpi=cfg.figure_dpi,
                   font_family=cfg.figure_font_family,
                   font_name=cfg.figure_font_name,
                   font_size=cfg.figure_font_size,
                   tickbox=True, grid=use_grid, is_invert_y_axe=False,
                   aspect='auto', guess_number_ticks=10,
                   # number_x_ticks=10, number_y_ticks=10
                   x_ticklabel_end=1e-9, y_ticklabel_end=1e-9
                  )

# add subplots
plot_r = plot.add_subplot_cartesian_2d(e_r_plot_name, 121, x_axe_label=x_axis_label, y_axe_label=y_axis_label)
plot_z = plot.add_subplot_cartesian_2d(e_z_plot_name, 122, x_axe_label=x_axis_label, y_axe_label=y_axis_label)


# In[ ]:


# add images
plot.add_image(e_r_plot_name, data_r, cmap=cmap, clim=clim_e_r)
plot.add_image(e_z_plot_name, data_z, cmap=cmap, clim=clim_e_z)

plot.add_colorbar(e_r_plot_name, ticks=clim_e_r, title=cbar_axis_label, size='10%')
plot.add_colorbar(e_z_plot_name, ticks=clim_e_z, title=cbar_axis_label, size='10%')

plot.redraw()


input("Please press RETURN to exit ")
