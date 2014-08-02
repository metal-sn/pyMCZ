import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='MCZ',
      version='1.0',
      author='Seung Man Oh',
      author_email='smo304@nyu.edu',
      description='MC simulation for getting metallicity confidence region',
      url='https://github.com/smanoh/MC_Metalicity',
      packages=['MCZ'],
     )
