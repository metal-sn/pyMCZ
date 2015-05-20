#!/usr/bin/env python
from distutils.core import setup

setup(name='mcz',
      version='1.0',
      author='Federica Bianco',
      author_email='fb55@nyu.edu',
      description='MC simulation for getting metallicity and confidence regions from line flux values',
      url='https://github.com/nyusngroup/MC_Metalicity',
      packages=['pyMCZ'],
      scripts=['pyMCZ/mcz.py', 'pyMCZ/metallicity.py','pyMCZ/metscales.py','pyMCZ/pylabsetup.py']
     )
