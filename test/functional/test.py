#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
import shutil
import subprocess
import gzip
import argparse

import jinja2 as j2

import colorama
from colorama import Fore, Back, Style
LOGLEVEL=0

me = os.path.realpath(os.path.dirname(__file__))

sys.path.append(me + '/../../tools/')

from lib.parameters import Parameters

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
                 accept_ieee=True, verbose=False):
        self.rootdir = os.path.realpath(os.path.dirname(os.path.dirname(me)))
        self.testdir = os.path.join(self.rootdir, testdir)
        self.keep_working_dir = keep_working_dir
        self.accept_ieee = accept_ieee
        self.verbose = verbose
        tmpldir = me

        utils = Util()

        print("Preparing the code")
        utils.cliexec('make distclean', cwd=self.rootdir, view=self.verbose, wait=True)
        utils.cliexec('./autogen.sh', cwd=self.rootdir, view=self.verbose, wait=True)
        if self.accept_ieee:
            utils.cliexec('./configure --enable-ieee', cwd=self.rootdir, view=self.verbose, wait=True)
        else:
            utils.cliexec('./configure', cwd=self.rootdir, view=self.verbose, wait=True)
        utils.cliexec('make', cwd=self.rootdir, view=self.verbose, wait=True)

        pdp3_file = os.path.join(self.rootdir, 'pdp3')
        tools_dir = os.path.join(self.rootdir, 'tools')

        print("Copying files to testing directory")
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

        isc = np.allclose(test_data, true_data, rtol=self.r_tol, atol=self.a_tol)

        if not isc and self.verbose:
            print("Difference in component {}:".format(component_name))
            for i in range(0, len(true_data)):
                if abs(test_data[i] - true_data[i]) > self.a_tol + self.r_tol * abs(true_data[i]):
                    print("number: {}, true data: {}, real data: {}".format(i, true_data[i], test_data[i]))
        if isc:
            print("Data matching for " + component_name + ":" + filename,
                  Fore.BLUE + "PASSED" + Style.RESET_ALL)
        else:
            print("Data matching for " + component_name + ":" + filename,
                  Fore.RED + "FAILED" + Style.RESET_ALL)
        return isc


def test_example(template_name, number, accept_ieee=True, rel_tolerance=0.01):
    verb=False
    status=True
    if LOGLEVEL > 0: verb=True

    b = bootstrap(testdir='testingdir',
                  parameters_template_name=template_name,
                  keep_working_dir=verb, # if verbose, keep working dir for analysis
                  accept_ieee=accept_ieee,
                  verbose=verb)

    t = pdp3Test(os.path.join(b.testdir, 'parameters.xml'),
                 rel_tolerance=rel_tolerance,
                 abs_tolerance=1) # use abs tolerance to avoid comparing of very small numbers
    t.verbose = verb

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
    t.components['Electrons'] = 'mpframe_100:100_150:150'
    t.components['Ions'] = 'mpframe_100:100_150:150'

    components = ['E_r', 'E_phi', 'E_z', 'H_r', 'H_phi', 'H_z', 'J_r', 'J_phi', 'J_z', 'rho_beam']
    for i in components:
        s = t.compare(i, '{}.dat'.format(number))
        if status: status = s

    for i in ['Electrons', 'Ions']:
        for j in ['vel', 'pos']:
            for k in ['r', 'z']:
                s = t.compare(i, "{}_{}_{}.dat".format(number, j, k))
                if status: status = s
        s = t.compare(i, "{}_vel_phi.dat".format(number))
        if status: status = s

    return status


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', type=str, help='smoke, ext', default='smoke')
    parser.add_argument('--fastmath', action='store_true',
                        help='Use fast math, which are not compatible with IEEE calculation standard',
                        default=False)
    parser.add_argument('-v', '--verbose', action='count', default=0)

    args = parser.parse_args()

    ieee = not args.fastmath
    LOGLEVEL = args.verbose
    status = True

    if args.type == 'smoke':
        if ieee:
            rtol = 0.001
        else:
            rtol = 0.1
        status = test_example('parameters.xml.tmpl', 4,
                              accept_ieee=ieee, rel_tolerance=rtol)
    elif args.type == 'ext':
        if ieee:
            rtol = 0.001
        else:
            rtol = 0.25
        status = test_example('parameters_ext.xml.tmpl', 11,
                              accept_ieee=ieee, rel_tolerance=rtol)
    else:
        raise Exception("there is no test type {}".format(args.type))

    if status:
        print(Fore.BLUE + "Test PASSED" + Style.RESET_ALL)
    else:
        print(Fore.RED + "Test FAILED" + Style.RESET_ALL)
