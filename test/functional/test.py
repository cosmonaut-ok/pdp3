#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import shutil
import subprocess
import gzip
import argparse
from scipy.fftpack import fft

import jinja2 as j2

import colorama
from colorama import Fore, Back, Style
LOGLEVEL=0

me = os.path.realpath(os.path.dirname(__file__))

regression_dir = os.path.join(os.sep, 'home', 'cosmonaut', 'regression')

sys.path.append(me + '/../../tools/')

from lib.parameters import Parameters
from lib.plot_builder import PlotBuilder
from lib.h5_reader import H5Reader
from lib.plain_reader import PlainReader

## end import path

class Util ():
    def __init__(self):
        True

    def cliexec (self, cmd, cwd=None, view=False, wait=False):
        if not cwd:
            cwd = me
        p = subprocess.Popen(cmd, cwd=cwd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if view:
            output = p.stdout.readline()
            while True:
                output = p.stdout.readline()
                if output.decode() == '' and p.poll() is not None:
                    break
                if output:
                    print(output.decode().strip())
            p.poll()

        if wait:
            p.wait()
        return p.returncode


    def mkdir (self, path, force=False):
        if not os.path.exists(path):
            os.makedirs(path)
        else:
            if force:
                try:
                    shutil.rmtree(path, ignore_errors=True)
                except:
                    os.remove(path)
                os.makedirs(path)
            else:
                raise OSError("Directory {} already exists".format(path))

    def cp (self, src, dst, recursively=True):
        if os.path.isfile(src):
            _dst = os.path.join(dst, os.path.basename(src))
            shutil.copyfile(src, _dst, follow_symlinks=True)
            shutil.copymode(src, _dst)
        else:
            if os.path.isdir(dst):
                _dst = os.path.join(dst, os.path.basename(src))
            else:
                _dst = dst
            shutil.copytree(src, _dst, symlinks=True, ignore=None,
                            ignore_dangling_symlinks=False)

class bootstrap ():
    def __init__(self, testdir='testdir', parameters_template_name='parameters.xml.tmpl',
                 result_path='.', data_root='pdp3_result', keep_working_dir=False,
                 accept_ieee=True, verbose=False, debug=False):
        self.rootdir = os.path.realpath(os.path.dirname(os.path.dirname(me)))
        self.testdir = os.path.join(self.rootdir, testdir)
        self.keep_working_dir = keep_working_dir
        self.accept_ieee = accept_ieee
        self.verbose = verbose or debug
        self.debug = debug
        tmpldir = me

        utils = Util()

        print("Preparing the code")
        utils.cliexec('make distclean', cwd=self.rootdir, view=self.verbose, wait=True)
        utils.cliexec('./autogen.sh', cwd=self.rootdir, view=self.verbose, wait=True)
        if self.accept_ieee:
            utils.cliexec('./configure --enable-ieee --enable-testmode', cwd=self.rootdir, view=self.verbose, wait=True)
        else:
            utils.cliexec('./configure --enable-testmode', cwd=self.rootdir, view=self.verbose, wait=True)
        utils.cliexec('make', cwd=self.rootdir, view=self.verbose, wait=True)

        pdp3_file = os.path.join(self.rootdir, 'pdp3')
        tools_dir = os.path.join(self.rootdir, 'tools')

        print("Copying files to testing directory")
        print(self.testdir)
        utils.mkdir(self.testdir, True)
        utils.cp(pdp3_file, self.testdir)
        utils.cp(tools_dir, self.testdir)

        tmpl_path = os.path.join(tmpldir, parameters_template_name)

        j2_env = j2.Environment(loader=j2.FileSystemLoader(tmpldir), trim_blocks=True)
        rendered_tmpl = j2_env.get_template(parameters_template_name).render(
            result_path=result_path,
            data_root=data_root
            )

        with open(os.path.join(self.testdir, 'parameters.xml'), "w") as f:
            f.write(rendered_tmpl)

        print("Launging application to prepare data for testing")

        if self.verbose: print("\nApplication Output:\n===================\n")
        utils.cliexec(os.path.join(self.testdir, 'pdp3'), cwd=self.testdir, view=self.verbose, wait=True)
        if self.verbose: print("\nEnd of application Output.\n==========================\n")


    def __del__(self):
        utils = Util()
        print("Clearing testing data")
        if not self.keep_working_dir:
            shutil.rmtree(self.testdir, ignore_errors=True)
        utils.cliexec('make clean', cwd=self.rootdir, view=False, wait=True)


class pdp3Test ():
    def __init__(self, config_file, rel_tolerance=1e-3, abs_tolerance=0):
        self.cfg = Parameters(config_file)
        self.components = {}
        self.rootdir = os.path.realpath(os.path.dirname(os.path.dirname(me)))
        self.testdir = os.path.join(self.rootdir, os.path.dirname(config_file))
        self.r_tol = rel_tolerance
        self.a_tol = abs_tolerance
        self.verbose = False

        self.true_data_dir = os.path.join(me, 'true_data')
        # self.frame_name = 'frame_0:0_254:2046'

    def compare(self, component_name, filename):
        frame_name = self.components[component_name]
        true_path = os.path.join(me, 'true_data', component_name, frame_name, filename)
        test_path = os.path.join(self.testdir, self.cfg.data_path, self.cfg.data_root,
                                 component_name, frame_name, filename)

        with gzip.open(true_path + '.gz', 'r') as myfile:
            true_data = np.fromstring(myfile.read(), dtype=float, sep=' ')

        test_data = np.fromfile(test_path, dtype=float, sep=' ')

        ## clear from NaNs. Required to temperature comparation
        index_true = np.argwhere(np.isnan(true_data))
        index_test = np.argwhere(np.isnan(test_data))
        true_data = np.delete(true_data, index_true)
        test_data = np.delete(test_data, index_test)

        isc = np.allclose(test_data, true_data, rtol=self.r_tol, atol=self.a_tol)

        if isc:
            print("Data matching for " + component_name + ":" + filename,
                  Fore.BLUE + "PASSED" + Style.RESET_ALL)
        else:
            print("Data matching for " + component_name + ":" + filename,
                  Fore.RED + "FAILED" + Style.RESET_ALL)

        if not isc and self.verbose:
            print("Difference in component {}:".format(component_name))
            for i in range(0, len(true_data)):
                if abs(test_data[i] - true_data[i]) > self.a_tol + self.r_tol * abs(true_data[i]):
                    print("number: {}, true data: {}, real data: {}".format(i, true_data[i], test_data[i]))

        return isc


def test_example(template_name, number, accept_ieee=True,
                 rel_tolerance=0.01, verbose=False, debug=False):
    status=True

    b = bootstrap(testdir='testingdir',
                  parameters_template_name=template_name,
                  keep_working_dir=debug, # if debug, keep working dir for analysis
                  accept_ieee=accept_ieee,
                  verbose=verbose,
                  debug=debug)

    t = pdp3Test(os.path.join(b.testdir, 'parameters.xml'),
                 rel_tolerance=rel_tolerance,
                 abs_tolerance=1) # use abs tolerance to avoid comparing of very small numbers
    t.verbose = verbose
    t.debug = debug

    t.components['E_r'] = 'frame_0:0_32:256'
    t.components['E_phi'] = 'frame_32:256_64:512'
    t.components['E_z'] = 'frame_64:512_96:768'

    t.components['H_r'] = 'frame_96:768_128:1024'
    t.components['H_phi'] = 'frame_128:1024_160:1280'
    t.components['H_z'] = 'frame_160:1280_192:2046'

    t.components['J_r'] = 'frame_192:1536_244:1792'
    t.components['J_phi'] = 'frame_244:1792_255:2047'
    t.components['J_z'] = 'frame_0:0_32:256'

    t.components['rho_beam'] = 'frame_0:0_254:2046'

    # compare positions and velocities
    t.components['T/Electrons'] = 'frame_100:100_150:150'
    t.components['T/Ions'] = 'frame_100:100_150:150'

    components = ['E_r', 'E_phi', 'E_z', 'H_r', 'H_phi', 'H_z', 'J_r', 'J_phi', 'J_z', 'rho_beam', 'T/Electrons', 'T/Ions']
    for i in components:
        s = t.compare(i, '{}.dat'.format(number))
        if status: status = s

    # for i in ['Electrons', 'Ions']:
    #     for j in ['vel', 'pos']:
    #         for k in ['r', 'z']:
    #             s = t.compare(i, "{}_{}_{}.dat".format(number, j, k))
    #             if status: status = s
    #     s = t.compare(i, "{}_vel_phi.dat".format(number))
    #     if status: status = s

    return status


def regression_test_example(template_name, accept_ieee=True, verbose=False, debug=False):
    status=True

    # regression_dir

    b = bootstrap(testdir='testingdir',
                  parameters_template_name=template_name,
                  keep_working_dir=verbose, # if verbose, keep working dir for analysis
                  accept_ieee=accept_ieee,
                  verbose=verbose,
                  debug=debug)


    utils = Util()

    cfg = Parameters(os.path.join(b.testdir, 'parameters.xml'))

    el_charge = 1.6e-19
    rho_beam_scale = 1
    clim_estimation = cfg.get_clim_estimation()
    shape=[0,0,100,500]
    temperature_shape=[90,1000,130,1200] # 90,130,1000,1200]

    timestamp=1e-8
    time_range=[cfg.start_time, cfg.end_time]

    use_cache=False
    use_grid=True
    cmap='terrain'
    temperature_component='electrons'

    clim_e_r = [-clim_estimation, clim_estimation]
    clim_e_z = [-clim_estimation, clim_estimation]
    clim_rho_beam = [-(cfg.bunch_density * el_charge * rho_beam_scale), 0]
    autoselect = True
    x_axis_label = r'$\mathit{Z (m)}$'
    y_axis_label = r'$\mathit{R (m)}$'
    temperature_x_axis_label = r'$\mathit{t (s)}$'
    temperature_y_axis_label = r'$\mathit{T ( K^\circ )}$'
    e_x_axis_label = r'$\mathit{t (s)}$'
    e_y_r_axis_label = r'$\mathit{E_r (\frac{V}{m})}$'
    e_y_z_axis_label = r'$\mathit{E_z (\frac{V}{m})}$'
    cbar_axis_label = r'$\frac{V}{m}$'
    cbar_bunch_density_axis_label = r'$m^{-3}$'

    e_r_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Radial\enspace Component}\enspace(E_r)$'
    e_z_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Longitudal\enspace Component}\enspace(E_z)$'
    rho_beam_plot_name = r'$\mathbf{Electron\enspace Bunch\enspace Density}\enspace (\rho_{bunch})$'
    temperature_plot_name = r'$\mathbf{Temperature-time \enspace Dependency}\enspace(T)$'
    e_e_r_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Radial\enspace Component}\enspace(E_r)$'
    e_e_z_plot_name = r'$\mathbf{Electrical\enspace Field\enspace Longitudal\enspace Component}\enspace(E_z)$'

    e_shape = [34, 341]

    r_scale = (shape[2] - shape[0]) / cfg.number_r_grid
    z_scale = (shape[3] - shape[1]) / cfg.number_z_grid

    reader = PlainReader(path = cfg.data_path,
                         data_root=cfg.data_root,
                         fullframe_size=[cfg.number_r_grid , cfg.number_z_grid],
                         fpds=cfg.frames_per_file,
                         use_cache=use_cache,
                         verbose=verbose)

    ############################################################################

    plot = PlotBuilder(shape[3] - shape[1], shape[2] - shape[0],
                       fig_color=cfg.figure_color,
                       fig_width=cfg.figure_width,
                       fig_height=cfg.figure_height,
                       fig_dpi=cfg.figure_dpi,
                       font_family=cfg.figure_font_family,
                       font_name=cfg.figure_font_name,
                       font_size=cfg.figure_font_size,
                       x_ticklabel_end=cfg.z_size * z_scale, y_ticklabel_end=cfg.r_size * r_scale,
                       tickbox=True, grid=use_grid, is_invert_y_axe=False,
                       aspect='equal', image_interpolation='nearest', guess_number_ticks=20)

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

    # plot.show()

    data_r = data_z = data_beam = []

    for probe in cfg.probes:
        frame = cfg.get_frame_number_by_timestamp(timestamp, probe.schedule)
        if (probe.type == 'frame') and (probe.r_start == shape[0]) and (probe.z_start == shape[1]) and(probe.r_end == shape[2]) and(probe.z_end == shape[3]):
            if probe.component == 'E_r': data_r = reader.get_frame('E_r', shape, frame)
            if probe.component == 'E_z': data_z = reader.get_frame('E_z', shape, frame)
            if probe.component == 'rho_beam': data_beam = reader.get_frame('rho_beam', shape, frame)

    # add timestamp to plot
    plot.get_figure().suptitle("Time: {:.2e} s".format(timestamp), x=.85, y=.95)

    # add images
    plot.add_image(e_r_plot_name, data_r, cmap=cmap, clim=clim_e_r)
    plot.add_image(e_z_plot_name, data_z, cmap=cmap, clim=clim_e_z)
    plot.add_image(rho_beam_plot_name, data_beam, cmap=cmap, clim=clim_rho_beam)

    plot.save(os.path.join(b.rootdir, 'image.png'))

    ############################################################################

    start_frame = None
    end_frame = None
    dump_interval = None
    temperature=[]

    temperature_r=[]
    temperature_phi=[]
    temperature_z=[]

    for probe in cfg.probes:
        if (probe.type == 'mpframe') and (probe.component == temperature_component) and (probe.r_start == temperature_shape[0]) and (probe.z_start == temperature_shape[1]) and (probe.r_end == temperature_shape[2]) and (probe.z_end == temperature_shape[3]):
            dump_interval = probe.schedule
            start_frame = cfg.get_frame_number_by_timestamp(time_range[0], dump_interval)
            end_frame = cfg.get_frame_number_by_timestamp(time_range[1], dump_interval) - 1
            for i in range(start_frame, end_frame):
                el_sum_r=0
                el_sum_phi=0
                el_sum_z=0
                # calculating element
                r = np.fromfile("{}/{}/{}/mpframe_{}:{}_{}:{}/{}_vel_r.dat".format(cfg.config_path, cfg.data_root, probe.component, temperature_shape[0], temperature_shape[1], temperature_shape[2], temperature_shape[3], i), dtype='float', sep=' ')
                phi = np.fromfile("{}/{}/{}/mpframe_{}:{}_{}:{}/{}_vel_phi.dat".format(cfg.config_path, cfg.data_root, probe.component, temperature_shape[0], temperature_shape[1], temperature_shape[2], temperature_shape[3], i), dtype='float', sep=' ')
                z = np.fromfile("{}/{}/{}/mpframe_{}:{}_{}:{}/{}_vel_z.dat".format(cfg.config_path, cfg.data_root, probe.component, temperature_shape[0], temperature_shape[1], temperature_shape[2], temperature_shape[3], i), dtype='float', sep=' ')
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
    plot_t = plot.add_subplot_cartesian_2d(temperature_plot_name, 111,
                                           x_axe_label=temperature_x_axis_label,
                                           y_axe_label=temperature_y_axis_label)

    plot_t.plot(data_timeline, temperature_r, label="r")
    plot_t.plot(data_timeline, temperature_phi, label="phi")
    plot_t.plot(data_timeline, temperature_z, label="z")
    plot_t.legend(loc='upper left')

    plot.save(os.path.join(b.rootdir, 'image_temperature.png'))

    ############################################################################

    # get data
    start_frame = None # cfg.get_frame_number_by_timestamp(time_range[0])
    end_frame = None # 460 # cfg.get_frame_number_by_timestamp(time_range[1])

    data_r = data_z = []

    for probe in cfg.probes:
        if (probe.type == 'dot') and (probe.r_start == e_shape[0]) and (probe.z_start == e_shape[1]):
            start_frame = cfg.get_frame_number_by_timestamp(time_range[0], probe.schedule)
            end_frame = cfg.get_frame_number_by_timestamp(time_range[1], probe.schedule) - 1
            if probe.component == 'E_r': data_r = reader.get_frame_range_dot('E_r', e_shape[0], e_shape[1], start_frame, end_frame)
            if probe.component == 'E_z': data_z = reader.get_frame_range_dot('E_z', e_shape[0], e_shape[1], start_frame, end_frame)

    data_timeline = np.linspace(time_range[0], time_range[1], len(data_r))

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
    plot_r = plot.add_subplot_cartesian_2d(e_e_r_plot_name, 121, x_axe_label=e_x_axis_label, y_axe_label=e_y_r_axis_label)
    plot_z = plot.add_subplot_cartesian_2d(e_e_z_plot_name, 122, x_axe_label=e_x_axis_label, y_axe_label=e_y_z_axis_label)

    plot_r.plot(data_timeline, data_r)
    plot_z.plot(data_timeline, data_z)

    plot.save(os.path.join(b.rootdir, 'image_e.png'))

    return status



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', type=str, help='smoke, ext', default='smoke')
    parser.add_argument('--fastmath', action='store_true',
                        help='Use fast math, which are not compatible with IEEE calculation standard',
                        default=False)
    parser.add_argument('-v', '--verbose', action='count', default=0)
    parser.add_argument('-d', '--debug', action='count', default=0)

    args = parser.parse_args()

    ieee = not args.fastmath
    status = True

    if args.type == 'smoke':
        if ieee:
            rtol = 0.001
        else:
            rtol = 0.1
        status = test_example('parameters.xml.tmpl', 4,
                              accept_ieee=ieee, rel_tolerance=rtol,
                              verbose=args.verbose,
                              debug=args.debug)
    elif args.type == 'ext':
        if ieee:
            rtol = 0.001
        else:
            rtol = 0.25
        status = test_example('parameters_ext.xml.tmpl', 11,
                              accept_ieee=ieee, rel_tolerance=rtol,
                              verbose=args.verbose,
                              debug=args.debug)
    elif args.type == 'regression':
        status = regression_test_example('parameters_regression.xml.tmpl',
                                         accept_ieee=ieee,
                                         verbose=args.verbose,
                                         debug=args.debug)
    else:
        raise Exception("there is no test type {}".format(args.type))

    if status:
        print(Fore.BLUE + "Test PASSED" + Style.RESET_ALL)
    else:
        print(Fore.RED + "Test FAILED" + Style.RESET_ALL)
