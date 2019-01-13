#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# define matplotlibplotting backend
# %matplotlib -l shows all available backends


# In[ ]:


import os
import numpy as np
from numpy import pi
from scipy.fftpack import fft

from lib.parameters import Parameters
from lib.plot_builder import PlotBuilder
from lib.h5_reader import H5Reader
from lib.plain_reader import PlainReader


# In[ ]:


##  configuration options
config_file = '../testdir/parameters.xml'

cfg = Parameters(config_file)

radius=0.01
longitude=0.01
time_range=[cfg.start_time, cfg.end_time]
use_grid=True
ylim=None
use_cache=False
verbose = True
autoselect = True

data_sets = ['E_r', 'E_z', 'H_phi']

x_axis_label = r'$\mathit{\omega}$'
y_axis_label = r'$\mathit{A.U.}$'
plot_name = r'$\mathbf{Multicomponent\enspace Spectrum}$'


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
start_frame = False
end_frame = False
row_number = cfg.get_row_by_radius(radius)
col_number = cfg.get_col_by_longitude(longitude)


# In[ ]:


# define plot builder
plot = PlotBuilder(0, 0, # let the system detects sizes automatically
                   fig_color=cfg.figure_color, 
                   fig_width=cfg.figure_width,
                   fig_height=cfg.figure_height, 
                   fig_dpi=cfg.figure_dpi,
                   font_family=cfg.figure_font_family,
                   font_name=cfg.figure_font_name,
                   font_size=12,
                   tickbox=True, grid=use_grid, is_invert_y_axe=False,
                   aspect='auto', guess_number_ticks=20,
                   # number_x_ticks=10, number_y_ticks=10
                   # x_ticklabel_end=1e-9, y_ticklabel_end=1e-9
                  )

# add subplots
the_plot = plot.add_subplot_cartesian_2d(plot_name, 111, x_axe_label=x_axis_label, y_axe_label=y_axis_label)

# set y-limits
if ylim is not None:
    the_plot.set_ylim(ylim)

plot.show()


# In[ ]:


data_dump_interval = 0

for ds in data_sets:
    # get data
    data = []
    for probe in cfg.probes:
        if (probe.type == 'dot') and (probe.r_start == row_number) and (probe.z_start == col_number):
            if probe.component == ds:
                print(start_frame, end_frame)
                start_frame = start_frame or cfg.get_frame_number_by_timestamp(time_range[0], probe.schedule)
                end_frame = end_frame or cfg.get_frame_number_by_timestamp(time_range[1], probe.schedule)
                data = reader.get_frame_range_dot(ds, row_number, col_number, start_frame, end_frame)
                data_dump_interval = probe.schedule
        

    if len(data) == 0 and autoselect:
        for probe in cfg.probes:
            if (probe.type == 'col') and (probe.z_start == col_number):
                if probe.component == ds:
                    start_frame = start_frame or cfg.get_frame_number_by_timestamp(time_range[0], probe.schedule)
                    end_frame = end_frame or cfg.get_frame_number_by_timestamp(time_range[1], probe.schedule)
                    print(start_frame, end_frame)
                    data = reader.get_frame_range_col(ds, col_number, start_frame, end_frame)[:, row_number]
                    data_dump_interval = probe.schedule
                break
            elif (probe.type == 'row') and (probe.r_start == row_number):
                if probe.component == ds:
                    start_frame = start_frame or cfg.get_frame_number_by_timestamp(time_range[0], probe.schedule)
                    end_frame = end_frame or cfg.get_frame_number_by_timestamp(time_range[1], probe.schedule)
                    data = reader.get_frame_range_row(ds, row_number, start_frame, end_frame)[:, col_number]
                    data_dump_interval = probe.schedule
                break
            elif (probe.type == 'frame') and (probe.r_start <= row_number) and (probe.z_start <= col_number) and (probe.r_end >= row_number) and (probe.z_end >= col_number):
                shape = [probe.r_start, probe.z_start, probe.r_end, probe.z_end]
                if probe.component == ds:
                    start_frame = start_frame or cfg.get_frame_number_by_timestamp(time_range[0], probe.schedule)
                    end_frame = end_frame or cfg.get_frame_number_by_timestamp(time_range[1], probe.schedule)
                    data = reader.get_frame_range(ds, shape, start_frame, end_frame)[:, row_number, col_number]
                    data_dump_interval = probe.schedule

    N = len(data)
    sampling_rate = 1 / (cfg.step_interval * data_dump_interval)

    # Nyquist Sampling Criteria
    T = 1 / sampling_rate # inverse of the sampling rate
    x_f = np.linspace(0.0, 1.0/(2.0*T), int(N/2))

    # FFT algorithm
    yr = fft(data) # "raw" FFT with both + and - frequencies
    if ds == 'H_phi':
        y_f = 2/N * np.abs(yr[0:np.int(N/2)]) * 100 # positive freqs only
    else:
        y_f = 2/N * np.abs(yr[0:np.int(N/2)]) # positive freqs only

    if ds == 'H_phi':
        _ds = 'H_\phi'
    else:
        _ds = ds
    the_plot.plot(x_f, y_f, label=r'${}$'.format(_ds))
the_plot.legend(loc='upper left')


input("Please press RETURN to exit ")
