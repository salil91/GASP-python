#!/usr/bin/env python

import os
import glob

from setuptools import setup, find_packages

module_dir = os.path.dirname(os.path.abspath(__file__))

if __name__ == "__main__":
    setup(
        name='GASP',
        version='1.0',
        description='Genetic algorithm for structure and phase prediction',
        long_description=open(os.path.join(module_dir, 'README.rst')).read(),
        url='https://github.com/salil91/GASP-python',
        author='Benjamin Revard, Venkata Surya Chaitanya Kolluru, Salil Bavdekar',
        author_email='salil.bavdekar@ufl.edu',
        license='MIT',
        packages=find_packages(),
        package_data={},
        zip_safe=False,
        install_requires=['pymatgen<=2020.12.31'],
        classifiers=['Programming Language :: Python :: 2.7',
                     "Programming Language :: Python :: 3",
                     "Programming Language :: Python :: 3.7",
                     'Development Status :: 4 - Beta',
                     'Intended Audience :: Science/Research',
                     'Operating System :: OS Independent',
                     'Topic :: Other/Nonlisted Topic',
                     'Topic :: Scientific/Engineering'],
        test_suite='nose.collector',
        tests_require=['nose'],
        scripts=glob.glob(os.path.join(module_dir, "gasp", "scripts", "*")),
        entry_points={
            'console_scripts': ['run_gasp = gasp.scripts.run:main']
        }   
    )
