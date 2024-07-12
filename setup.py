#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

#try:
#  from setuptools import setup, find_packages
#  setup
#except ImportError:
#  from distutils.core import setup
#  setup

#from distutils.core import setup

from __future__ import unicode_literals
from setuptools import setup, find_packages
from core import VERSION


setup(
  name = "KASTREDUX",
  version = VERSION,
  packages = find_packages(),
#  packages = find_packages(exclude=['docs','tests']),

  # Project uses reStructuredText, so ensure that the docutils get
  # installed or upgraded on the target machine
  install_requires = [
    'astropy',
    'matplotlib',
    'numpy',
    'pandas',
    'scipy'
  ],

  package_dir = {'kastredux': 'kastredux'},    
  package_data = {'kastredux': ['docs/*','resources/*','training/*']},
  include_package_data=True,

  zip_safe = True,
  use_2to3 = False,
  classifiers=[
      'Development Status :: 3 - Alpha',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: MIT License',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 2',
      'Programming Language :: Python :: 2.6',
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3',
      'Programming Language :: Python :: 3.3',
      'Programming Language :: Python :: 3.4',
      'Programming Language :: Python :: 3.5',
      'Topic :: Scientific/Engineering :: Astronomy',
      'Topic :: Scientific/Engineering :: Physics'
  ],

  # metadata for upload to PyPI
  author = "Adam Burgasser",
  author_email = "aburgasser@ucsd.edu",
  description = "KAST Optical Reduction Toolkit",
#  long_description = long_description,
  license = "MIT",
#    download_url='%s/astropy-%s.tar.gz' % (DOWNLOAD_BASE_URL, VERSION),
  keywords = ['kast','spectroscopy', 'spectral reduction', 'astronomy','astrophysics'],
  #url = "http://www.browndwarfs.org/splat/",   # project home page, if any


)
