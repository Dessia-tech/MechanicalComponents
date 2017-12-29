# -*- coding: utf-8 -*-
"""
@author: Steven Masfaraud
"""

from setuptools import setup
import re

def readme():
    with open('README.rst') as f:
        return f.read()

with open('mechanical_components/__init__.py','r') as f:
    metadata = dict(re.findall("__([a-z]+)__\s*=\s*'([^']+)'", f.read()))

#print(metadata)

#import powertransmission
setup(name='mechanical_components',
      version=metadata['version'],
      description="Design of elementary components by AI",
      long_description='',
      keywords='',
      url='',
#      cmdclass['register']=None,
      author='Steven Masfaraud',
      author_email='masfaraud@dessia.tech',
      packages=['mechanical_components'],
      install_requires=['numpy','scipy','volmdlr'])

