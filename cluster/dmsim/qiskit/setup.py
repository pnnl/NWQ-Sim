# -*- coding: utf-8 -*-

# Copyright 2018, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

import os
import platform
from distutils.command.build import build
from subprocess import call
from distutils.sysconfig import get_python_lib

from setuptools.dist import Distribution
from setuptools import setup, find_packages

class NWQSimSimulatorBuild(build):
    def run(self):
        super().run()
        # Store the current working directory, as invoking cmake involves
        # an out of source build and might interfere with the rest of the steps.
        current_directory = os.getcwd()
        BUILD_DIR = "build/lib/"
        try:
            cmd_cmake = ['cmake', '']
            cmd_cmake.append('.')
            cmd_cmake.append('-B{}'.format(BUILD_DIR))

            make_exe = 'make'
            cmd_make = [make_exe, '-C', BUILD_DIR]

            def compile_simulator():
                print(" ".join(cmd_cmake))
                call(cmd_cmake)
                call(cmd_make)
                copy_so_files() 

            def copy_so_files():
                try:
                    from shutil import copyfile
                    build_dir = os.path.abspath(BUILD_DIR)
                    svsim_so_src = os.path.join(build_dir, "libdmsim.so")
                    #print ("current path")
                    #print (os.getcwd())
                    #print ("src path")
                    #print (svsim_so_src)
                    #print ("dst path")
                    svsim_so_dst = os.path.join(build_dir, "./qiskit_nwqsim_provider/libdmsim.so")
                    #print (svsim_so_dst)
                    try:
                        copyfile(svsim_so_src, svsim_so_dst)
                    except:
                        print("Error copying libdmsim.so.")

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
    python_requires='>3.8.0',
    version=0.1,
    author="Ang Li",
    author_email="ang.li@pnnl.gov",
    description="Qiskit simulator whose backend is based on the NWQSim simulator",
    long_description = "This module contains [Qiskit](https://www.qiskit.org/) simulator whose backend is written in NWQSim simulator. This simulator simulate a Quantum circuit on a classical computer.",
    url="https://github.com/pnnl/dmsim/qiskit-nwqsim-provider",
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
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering",
    ],
    install_requires=['qiskit-terra>=0.8'],
    keywords="qiskit quantum nwqsim dmsim svsim",
    packages=find_packages(exclude=['test*']),
    include_package_data=True,
    options={'build':{'build_lib':'build/lib'}},

    cmdclass={
        'build': NWQSimSimulatorBuild,
    },
    distclass=BinaryDistribution
)
