#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# define matplotlibplotting backend
# %matplotlib -l shows all available backends


# In[ ]:


import numpy as np
import math

from lib.parameters import Parameters
from lib.plot_builder import PlotBuilder
# from lib.h5_reader import H5Reader
# from lib.plain_reader import PlainReader


# In[ ]:


config_file = '/media/nfs/PDP3/model35r1/parameters.xml'

cfg = Parameters(config_file)

shape=[130,180,160,280]

use_grid=True
ylim_r=None
ylim_z=None
use_cache=False
verbose=True
set_ylim=False
component='electrons'
# time_range = [0, cfg.end_time]
time_range = [0, 1e-8]

x_axis_label = r'$\mathit{t (s)}$'
y_axis_label = r'$\mathit{T ( K^\circ )}$'

plot_name = r'$\mathbf{Temperature-time \enspace Dependency}\enspace(T)$'


# In[ ]:


start_frame = None
end_frame = None
dump_interval = None
temperature=[]

temperature_r=[]
temperature_phi=[]
temperature_z=[]

for probe in cfg.probes:
    if (probe.type == 'mpframe') and (probe.component == component) and (probe.r_start == shape[0]) and (probe.z_start == shape[1]) and (probe.r_end == shape[2]) and (probe.z_end == shape[3]):
        dump_interval = probe.schedule
        start_frame = cfg.get_frame_number_by_timestamp(time_range[0], dump_interval)
        end_frame = cfg.get_frame_number_by_timestamp(time_range[1], dump_interval)
        for i in range(start_frame, end_frame):
            el_sum_r=0
            el_sum_phi=0
            el_sum_z=0
            # calculating element
            r = np.fromfile("{}/{}/{}/mpframe_{}:{}_{}:{}/{}_vel_r.dat".format(cfg.config_path, cfg.data_root, probe.component, shape[0], shape[1], shape[2], shape[3], i), dtype='float', sep=' ')
            phi = np.fromfile("{}/{}/{}/mpframe_{}:{}_{}:{}/{}_vel_phi.dat".format(cfg.config_path, cfg.data_root, probe.component, shape[0], shape[1], shape[2], shape[3], i), dtype='float', sep=' ')
            z = np.fromfile("{}/{}/{}/mpframe_{}:{}_{}:{}/{}_vel_z.dat".format(cfg.config_path, cfg.data_root, probe.component, shape[0], shape[1], shape[2], shape[3], i), dtype='float', sep=' ')
            print("processing frame", i)
            for i in range(0, len(r)-1):
                # el_sum_sq += r[i]*r[i]+phi[i]*phi[i]+z[i]*z[i]
                el_sum_r += abs(r[i])
                el_sum_phi += abs(phi[i])
                el_sum_z += abs(z[i])
                
            # temperature.append(math.sqrt(el_sum_sq / len(r)))
            temperature_r.append(el_sum_r / len(r))
            temperature_phi.append(el_sum_phi / len(r))
            temperature_z.append(el_sum_z / len(r))


# In[ ]:


data_timeline = np.linspace(time_range[0], time_range[1], len(temperature_r))

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
                   x_ticklabel_start=time_range[0],
                   x_ticklabel_end=time_range[1]
                  )

# add subplots
plot_t = plot.add_subplot_cartesian_2d(plot_name, 111, x_axe_label=x_axis_label, y_axe_label=y_axis_label)

# set y-limits
if ylim_r is not None:
    plot_t.subplot_er.set_ylim(ylim_r)


# In[ ]:


plot_t.plot(data_timeline, temperature_r, label="r")
plot_t.plot(data_timeline, temperature_phi, label="phi")
plot_t.plot(data_timeline, temperature_z, label="z")
plot_t.legend(loc='upper left')

plot.show()


input("Please press RETURN to exit ")
