#!/usr/bin/env python
# coding: utf-8

# In[1]:


# define matplotlibplotting backend
# %matplotlib -l shows all available backends


# In[2]:


import os
import numpy as np

from lib.parameters import Parameters
from lib.plot_builder import PlotBuilder
from lib.h5_reader import H5Reader
from lib.plain_reader import PlainReader


# In[3]:


##  configuration options
config_file = '../workdir/parameters.xml'

cfg = Parameters(config_file)
clim_estimation = cfg.get_clim_estimation()

radius=0.01
time_range=[cfg.start_time, cfg.end_time]
use_grid=True
cmap='terrain'
clim = [-1e3, 1e3]
use_cache = False
verbose = True
autoselect = True

x_axis_label = r'$\mathit{Z (m)}$'
y_axis_label = r'$\mathit{t (ns)}$'
cbar_axis_label = r'$\frac{A}{m}$'
plot_name = r'$\mathbf{Magneitc\enspace Field\enspace Rotational\enspace Component}\enspace(H_{\phi})$'


# In[4]:


# define reader (plain reader used)
if not cfg.use_hdf5:
    reader = PlainReader(path = cfg.data_path,
                         data_root=cfg.data_root,
                         fullframe_size=[cfg.number_r_grid, cfg.number_z_grid],
                         fpds=cfg.frames_per_file, 
                         use_cache=use_cache,
                         verbose=verbose)
else:
    reader = H5Reader(str(os.path.join(cfg.data_path, 'data.h5')), use_cache=use_cache)
    reader.verbose = True


# In[5]:


# get data
start_frame = cfg.get_frame_number_by_timestamp(time_range[0])
end_frame = cfg.get_frame_number_by_timestamp(time_range[1])
row_number = 1 # cfg.get_row_by_radius(radius)

data = []

# data = reader.get_frame_range_row('H_phi', row_number, start_frame, end_frame - 1)

for probe in cfg.probes:
    if (probe.type == 'row') and (probe.r_start == row_number):
        if probe.component == 'H_phi': data = reader.get_frame_range_row('H_phi', row_number, start_frame, end_frame)
        
if len(data) == 0 and autoselect:
    for probe in cfg.probes:
        if (probe.type == 'frame') and (probe.r_start <= row_number) and (probe.z_start == 0) and (probe.r_end >= row_number) and (probe.z_end == cfg.number_z_grid):
            shape = [probe.r_start, probe.z_start, probe.r_end, probe.z_end]
            if probe.component == 'H_phi': data = reader.get_frame_range('H_phi', shape, start_frame, end_frame)[:, row_number]


# In[6]:


data_timeline = np.linspace(time_range[0], time_range[1], len(data))

# define plot builder
plot = PlotBuilder(0, 0, # let the system detects sizes automatically
                   fig_color=cfg.figure_color, 
                   fig_width=cfg.figure_width,
                   fig_height=cfg.figure_height, 
                   fig_dpi=cfg.figure_dpi,
                   font_family=cfg.figure_font_family,
                   font_name=cfg.figure_font_name,
                   font_size=cfg.figure_font_size,
                   tickbox=True, grid=use_grid, is_invert_y_axe=False,
                   aspect='auto', guess_number_ticks=20,
                   # number_x_ticks=10, number_y_ticks=10
                   x_ticklabel_end=1e-9, y_ticklabel_end=1e-9
                  )

# add subplots
the_plot = plot.add_subplot_cartesian_2d(plot_name, 111, x_axe_label=x_axis_label, y_axe_label=y_axis_label)


# In[7]:


# add images
plot.add_image(plot_name, data, cmap=cmap, clim=clim)
plot.add_colorbar(plot_name, ticks=clim, title=cbar_axis_label)

plot.show()


input("Please press RETURN to exit ")
