#!/usr/bin/env python

import math

# import sys
# import os

import argparse

# from numpy import *
# from pylab import *

# import matplotlib.pyplot as plt
# import matplotlib.animation as ani

from parameters import Parameters
from pdp3_plot_builder import PDP3PlotBuilder

def langmur_freq(density):
    # pi = 3.1415
    charge_el = -1.6e-19 # coul
    mass_el = 9.1e-31 # kg
    epsilon_0 = 8.8e-12
    ##
    w_p = math.sqrt(density * math.pow(charge_el, 2) / (mass_el * epsilon_0))
    return w_p

def debye_length(density, temperature):
    charge_el = -1.6e-19 # coul
    mass_el = 9.1e-31 # kg
    epsilon_0 = 8.8e-12
    boltzman_const = 1.38e-23

    l_d = epsilon_0 * boltzman_const * temperature / (4 * math.pi * density * math.pow(charge_el, 2))
    l_d = math.sqrt(l_d)
    return l_d

def eV2kelvin(ev):
    res = ev * 1.16045221e4
    return res

def wake_length(density, beam_velocity):
    w_p = langmur_freq(density)
    lambda_ = 2* math.pi * beam_velocity / w_p
    return lambda_

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('properties_path', metavar='properties_path', type=str,
                        help='Full path to properties.xml')

    args = parser.parse_args()

    config = Parameters(args.properties_path)

    for i in config.particles:
        p_name = i.getAttribute('name').lower()
        if p_name == 'electrons':
            left_el_density = float(i.getElementsByTagName('left_density')[0].firstChild.data)
            right_el_density = float(i.getElementsByTagName('right_density')[0].firstChild.data)
            el_density = (left_el_density + right_el_density) / 2

            el_temperature = float(i.getElementsByTagName('temperature')[0].firstChild.data)

    if not el_density:
        print("Incorrect parameters.xml file (no particles->particle_kind->electrons section)")
    else:
        bunch_duration = config.bunch_length / config.beam_initial_velocity
        beam_duration = bunch_duration * config.number_bunches + config.bunches_distance * (config.number_bunches -1)
        w_p = langmur_freq(el_density)
        wake_len = wake_length(el_density, config.beam_initial_velocity)
        el_temperature = eV2kelvin(el_temperature)
        debye = debye_length(el_density, el_temperature)
        beam_part_number = math.pi * math.pow(config.bunch_radius, 2) * bunch_duration * config.beam_initial_velocity * config.bunch_density

        print("Expected plasma frequency is %.4g Hz"%(w_p/(2 * math.pi)))
        print("Expected wake wavelength is %.2g m"%(wake_len))
        print("Expected Debye length is %.4g m"%(debye))
        print("Initial beam velocity is %.4g m/s"%(config.beam_initial_velocity))
        print("Single particles bunch duration is %.4g m"%bunch_duration)
        print("Whole particles beam duration is %.4g m"%beam_duration)
        print("Number of particles in beam is %.2g"%(beam_part_number))
        print()
        print("Estimated beam passage time is %.4g s"%(config.z_size / config.beam_initial_velocity + beam_duration / config.beam_initial_velocity))
        print("Number of calculation steps is %.0g"%((config.end_time - config.start_time) / config.step_interval))

## call main function
if __name__ == "__main__":
    main()
