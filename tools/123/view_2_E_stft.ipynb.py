#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# define matplotlibplotting backend
# %matplotlib -l shows all available backends


# In[ ]:


import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.signal import stft

from mpl_toolkits.axes_grid1 import make_axes_locatable

from lib.parameters import Parameters
from lib.plot_builder import PlotBuilder
from lib.h5_reader import H5Reader
from lib.plain_reader import PlainReader


# In[ ]:


##  configuration options
config_file = '../workdir/parameters.xml'

cfg = Parameters(config_file)
clim_estimation = cfg.get_clim_estimation()

radius=0.01
longitude=0.05
time_range=[cfg.start_time, cfg.end_time]
use_grid=True
cmap='terrain'
# clim_e_r = [-clim_estimation, clim_estimation]
# clim_e_z = [-clim_estimation, clim_estimation]
clim_e_r = [0,0]
clim_e_z = [0,0]
use_cache = False
verbose = True
autoselect = True

cmap = 'terrain'

# N dots per segment
nperseg=254

x_axis_label = r'$\mathit{Time (by window)}$'
y_axis_label = r'$\mathit{F (Hz)}$'
cbar_axis_label = r'$Amplitude (A.U.)$'

e_r_plot_name = r'$\mathbf{E_r\enspace Short\enspace Time\enspace Spectra}$'
e_z_plot_name = r'$\mathbf{E_z\enspace Short\enspace Time\enspace Spectra}$'

cbar_axis_label = r'$\frac{V}{m}$'


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
    reader = H5Reader(str(os.path.join(cfg.data_path, 'data.h5')),
                      shape=[cfg.number_r_grid , cfg.number_z_grid],
                      use_cache=use_cache)
    reader.verbose = True


# In[ ]:


# get data
start_frame = cfg.get_frame_number_by_timestamp(time_range[0])
end_frame = cfg.get_frame_number_by_timestamp(time_range[1])
row_number = cfg.get_row_by_radius(radius)
col_number = cfg.get_col_by_longitude(longitude)

data_r = data_z = []

for probe in cfg.probes:
    if (probe.type == 'dot') and (probe.r_start == row_number) and (probe.z_start == col_number):
        if probe.component == 'E_r': data_r = reader.get_frame_range_dot('E_r', row_number, col_number, start_frame, end_frame)
        if probe.component == 'E_z': data_z = reader.get_frame_range_dot('E_z', row_number, col_number, start_frame, end_frame)
        
if (len(data_r) == 0 or len(data_z) == 0) and autoselect:
    for probe in cfg.probes:
        if (probe.type == 'col') and (probe.z_start == col_number):
            if probe.component == 'E_r': data_r = reader.get_frame_range_col('E_r', col_number, start_frame, end_frame)[:, row_number]
            if probe.component == 'E_z': data_z = reader.get_frame_range_col('E_z', col_number, start_frame, end_frame)[:, row_number]
            break
        elif (probe.type == 'row') and (probe.r_start == row_number):
            if probe.component == 'E_r': data_r = reader.get_frame_range_row('E_r', row_number, start_frame, end_frame)[:, col_number]
            if probe.component == 'E_z': data_z = reader.get_frame_range_row('E_z', row_number, start_frame, end_frame)[:, col_number]
            break
        elif (probe.type == 'frame') and (probe.r_start <= row_number) and (probe.z_start <= col_number) and (probe.r_end >= row_number) and (probe.z_end >= col_number):
            shape = [probe.r_start, probe.z_start, probe.r_end, probe.z_end]
            if probe.component == 'E_r': data_r = reader.get_frame_range('E_r', shape, start_frame, end_frame)[:, row_number, col_number]
            if probe.component == 'E_z': data_z = reader.get_frame_range('E_z', shape, start_frame, end_frame)[:, row_number, col_number]
                
# data_r = reader.get_frame_range_dot('E_r', row_number, col_number, start_frame, end_frame - 1)
# data_z = reader.get_frame_range_dot('E_z', row_number, col_number, start_frame, end_frame - 1)


# In[ ]:


## sampling frequency (chastota dyskretyzacii)
fs = 1 / (cfg.step_interval * cfg.data_dump_interval)

f_r, t_r, Zxx_r = stft(data_r, fs, nperseg=nperseg)
f_z, t_z, Zxx_z = stft(data_z, fs, nperseg=nperseg)
color_map = plt.get_cmap(cmap)
print(len(t_r))


# In[ ]:


# define plot builder
plot = PlotBuilder(0, 0,
                   fig_color=cfg.figure_color,
                   fig_width=cfg.figure_width,
                   fig_height=cfg.figure_height,
                   fig_dpi=cfg.figure_dpi,
                   font_family=cfg.figure_font_family,
                   font_name=cfg.figure_font_name,
                   font_size=cfg.figure_font_size,
                   tickbox=True, grid=use_grid, is_invert_y_axe=False,
                   aspect='auto', guess_number_ticks=5
                   # number_x_ticks=10, number_y_ticks=10,
                   # x_ticklabel_end=1e8, y_ticklabel_end=1e9
                  )

# add subplots
plot_r = plot.add_subplot_cartesian_2d(e_r_plot_name, 121, x_axe_label=x_axis_label, y_axe_label=y_axis_label)
plot_z = plot.add_subplot_cartesian_2d(e_z_plot_name, 122, x_axe_label=x_axis_label, y_axe_label=y_axis_label)


# In[ ]:


if clim_e_r == [0, 0]:
    im0 = plot_r.pcolormesh(t_r, f_r, np.abs(Zxx_r), cmap=color_map)
else:
    im0 = plot_r.pcolormesh(t_r, f_r, np.abs(Zxx_r),
                                cmap=color_map, vmin=clim_e_r[0],
                                vmax=clim_e_r[1])

divider = make_axes_locatable(plot_r)
cax = divider.append_axes(position='right', size='10%', pad=0.05)
plot.get_figure().colorbar(im0, cax=cax, format='%.2e') # , ax = axes)

if clim_e_z == [0, 0]:
    im1 = plot_z.pcolormesh(t_z, f_z, np.abs(Zxx_z), cmap=color_map)
else:
    im1 = plot_z.pcolormesh(t_z, f_z, np.abs(Zxx_z),
                                cmap=color_map, vmin=clim_e_z[0],
                                vmax=clim_e_z[1])

divider = make_axes_locatable(plot_z)
cax = divider.append_axes(position='right', size='10%', pad=0.05)
plot.get_figure().colorbar(im0, cax=cax, format='%.2e') # , ax = axes)


input("Please press RETURN to exit ")
