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

radius=0.01
longitude=0.01
time_range=[cfg.start_time, cfg.end_time]
use_grid=True
ylim=None
use_cache = False
verbose = True
autoselect = True

x_axis_label = r'$\mathit{t (s)}$'
y_axis_label = r'$\mathit{H_{\phi} (\frac{A}{m})}$'

plot_name = r'$\mathbf{Magnetic\enspace Field\enspace Rotational\enspace Component}\enspace(H_{\phi})$'


# In[ ]:


# define reader (plain reader used)
if not cfg.use_hdf5:
    reader = PlainReader(path = cfg.data_path,
                         data_root=cfg.data_root,
                         fullframe_size=[cfg.number_r_grid , cfg.number_z_grid],
                         fpds=cfg.frames_per_file, 
                         use_cache=use_cache,
                         verbose=verbose)
else:
    reader = H5Reader(str(os.path.join(cfg.data_path, 'data.h5')), use_cache=use_cache)
    reader.verbose = True


# In[ ]:


# get data
start_frame = cfg.get_frame_number_by_timestamp(time_range[0])
end_frame = cfg.get_frame_number_by_timestamp(time_range[1])
row_number = cfg.get_row_by_radius(radius)
col_number = cfg.get_col_by_longitude(longitude)

# data = reader.get_frame_range_dot('H_phi', row_number, col_number, start_frame, end_frame - 1)

data = []

for probe in cfg.probes:
    if (probe.type == 'dot') and (probe.r_start == row_number) and (probe.z_start == col_number):
        if probe.component == 'H_phi': data = reader.get_frame_range_dot('H_phi', row_number, col_number, start_frame, end_frame)
        
if len(data) == 0 and autoselect:
    for probe in cfg.probes:
        if (probe.type == 'col') and (probe.z_start == col_number):
            if probe.component == 'H_phi': data = reader.get_frame_range_col('H_phi', col_number, start_frame, end_frame)[:, row_number]
            break
        elif (probe.type == 'row') and (probe.r_start == row_number):
            if probe.component == 'H_phi': data = reader.get_frame_range_row('H_phi', row_number, start_frame, end_frame)[:, col_number]
            break
        elif (probe.type == 'frame') and (probe.r_start <= row_number) and (probe.z_start <= col_number) and (probe.r_end >= row_number) and (probe.z_end >= col_number):
            shape = [probe.r_start, probe.z_start, probe.r_end, probe.z_end]
            if probe.component == 'H_phi': data = reader.get_frame_range('H_phi', shape, start_frame, end_frame)[:, row_number, col_number]


# In[ ]:


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
plot_h_phi = plot.add_subplot_cartesian_2d(plot_name, 111, x_axe_label=x_axis_label, y_axe_label=y_axis_label)

# set y-limits
if ylim is not None:
    plot_r.subplot_er.set_ylim(ylim_r)


# In[ ]:


# add data
plot_h_phi.plot(data_timeline, data)

plot.show()


input("Please press RETURN to exit ")
