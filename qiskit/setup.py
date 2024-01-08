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
        current_directory = os.getcwd()
        BUILD_DIR = "./build"
        try:
            cmd_cmake = ['cmake', '']
            cmd_cmake.append('../')
            cmd_cmake.append('-B{}'.format(BUILD_DIR))
            make_exe = 'make'
            cmd_make = [make_exe, '-C', BUILD_DIR]

            def compile_simulator():
                print(" ".join(cmd_cmake))
                call(cmd_cmake)
                call(cmd_make)
                copy_nwqsim_exe() 

            def copy_nwqsim_exe():
                try:
                    from shutil import copyfile
                    build_dir = os.path.abspath('.')
                    #print (build_dir)
                    #print (nwq_src)
                    #print (nwq_dst)
                    nwq_src = os.path.join(build_dir, "build/qasm/nwq_qasm")
                    nwq_dst = os.path.join(build_dir, "build/qiskit_nwqsim_provider/nwq_qasm")
                    try:
                        copyfile(nwq_src, nwq_dst)
                    except:
                        print("Error copying nwq_qasm")
                except Exception as copy_exception:
                    print("WARNING: EXE files copy failed: {}".format(copy_exception))
            self.execute(compile_simulator, [], 'Compiling NWQSim Simulator')
        except Exception as e:
            print(str(e))
            print("WARNING: it seems NWQSim simulator can't be built.")

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
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering",
    ],
    install_requires=['qiskit-terra>=0.8'],
    keywords="qiskit quantum nwqsim dmsim svsim",
    packages=find_packages(exclude=['test*']),
    include_package_data=True,
    options={'build':{'build_lib':'build'}},

    cmdclass={
        'build': NWQSimSimulatorBuild,
    },
    distclass=BinaryDistribution,
    zip_safe=False
)
