# -*- coding: utf-8 -*-

# Copyright 2018, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

import os
import platform
from distutils.command.build import build
from subprocess import call

from setuptools.dist import Distribution
from setuptools import setup, find_packages

VERSION_PATH = os.path.join(os.path.dirname(__file__), "VERSION.txt")
with open(VERSION_PATH, "r") as version_file:
    VERSION = version_file.read().strip()

class NWQSimSimulatorBuild(build):
    def run(self):
        super().run()
        # Store the current working directory, as invoking cmake involves
        # an out of source build and might interfere with the rest of the steps.
        current_directory = os.getcwd()
        BUILD_DIR = "build/lib/qiskit_nwqsim_provider"
        try:
            supported_platforms = ['Linux', 'Darwin', 'Windows']
            current_platform = platform.system()
            if current_platform not in supported_platforms:
                # TODO: stdout is silenced by pip if the full setup.py invocation is
                # successful, unless using '-v' - hence the warnings are not printed.
                print('WARNING: NWQSim simulator is meant to be built with these '
                      'platforms: {}.'
                      .format(supported_platforms))
                return

            cmd_cmake = ['cmake', '']
            if 'USER_LIB_PATH' in os.environ:
                cmd_cmake.append('-DUSER_LIB_PATH={}'.format(os.environ['USER_LIB_PATH']))
            if current_platform == 'Windows':
                # We only support MinGW so far
                cmd_cmake.append("-GMinGW Makefiles")
            cmd_cmake.append('.')
            cmd_cmake.append('-B{}'.format(BUILD_DIR))

            make_exe = 'make'
            if current_platform == 'Windows':
                make_exe = 'mingw32-make'

            cmd_make = [make_exe, '-C', BUILD_DIR]

            def compile_simulator():
                print(" ".join(cmd_cmake))
                call(cmd_cmake)
                call(cmd_make)
                if current_platform == 'Windows':
                    copy_dll_files() 

            def copy_dll_files():
                try:
                    from shutil import copyfile
                    build_dir = os.path.abspath(BUILD_DIR)
                    #print("-- Checking for enviornment variable MPFRDIR... {}".format("MPFRDIR" in os.environ))
                    #mpfr_dll = os.path.join(os.environ["MPFRDIR"], "mpfr.dll")
                    #mpfr_dll_dst = os.path.join(build_dir, "mpfr.dll")
                    #print("-- Checking for enviornment variable GMPDIR... {}".format("GMPDIR" in os.environ))
                    #mpir_dll = os.path.join(os.environ["GMPDIR"], "mpir.dll")
                    #mpir_dll_dst = os.path.join(build_dir, "mpir.dll")
                    #print("-- Copying {} to {}".format(mpfr_dll, mpfr_dll_dst))
                    #copyfile(mpfr_dll, mpfr_dll_dst)
                    #print("-- Copying {} to {}".format(mpir_dll, mpir_dll_dst))
                    #copyfile(mpir_dll, mpir_dll_dst)

                except Exception as copy_exception:
                    print("WARNING: DLL files copy failed: {}".format(copy_exception))


            self.execute(compile_simulator, [], 'Compiling NWQSim Simulator')
        except Exception as e:
            print(str(e))
            print("WARNING: Seems like NWQSim simulator can't be built.")

        # Restore working directory.
        os.chdir(current_directory)

# This is for creating wheel specific platforms
class BinaryDistribution(Distribution):
    def has_ext_modules(self):
        return True
    
setup(
    name="qiskit-nwqsim-provider",
    python_requires='>3.5.0',
    version=VERSION,
    author="Ang Li",
    author_email="ang.li@pnnl.gov",
    description="Qiskit simulator whose backend is based on the NWQSim simulator",
    long_description = "This module contains [Qiskit](https://www.qiskit.org/) simulator whose backend is written in NWQSim simulator. This simulator simulate a Quantum circuit on a classical computer.",
    url="https://github.com/pnnl/svsim/qiskit-nwqsim-provider",
    license="MIT",
    classifiers=[
        "Environment :: Console",
        "License :: OSI Approved :: MIT Software License",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering",
    ],
    install_requires=['qiskit-terra>=0.8'],
    keywords="qiskit quantum nwqsim dmsim svsim",
    packages=find_packages(exclude=['test*']),
    include_package_data=True,
    cmdclass={
        'build': NWQSimSimulatorBuild,
    },
    distclass=BinaryDistribution
)
