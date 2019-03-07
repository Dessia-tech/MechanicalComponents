# -*- coding: utf-8 -*-
"""
@author: Steven Masfaraud
"""

from setuptools import setup

import re

def readme():
    with open('README.rst') as f:
        return f.read()

from os.path import dirname, isdir, join
import re
from subprocess import CalledProcessError, check_output


tag_re = re.compile(r'\btag: %s([0-9][^,]*)\b')
version_re = re.compile('^Version: (.+)$', re.M)


def get_version():
    # Return the version if it has been injected into the file by git-archive
    version = tag_re.search('$Format:%D$')
    if version:
        return version.group(1)

    d = dirname(__file__)
    
    if isdir(join(d, '.git')):
        cmd = 'git describe --tags  --dirty'
        try:
            version = check_output(cmd.split()).decode().strip()[:]
        except CalledProcessError:
            raise RuntimeError('Unable to get version number from git tags')
        if version[0]=='v':
            version = version[1:]
        # PEP 440 compatibility
        if '-' in version:
            if version.endswith('-dirty'):
                version = '.dev'.join(version.split('-')[:-1][:2])+'-dirty'
        else:
            version = '.dev'.join(version.split('-')[:2])

    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)
            
    # Writing to file
    with open('mechanical_components/version.py', 'w+') as vf:
        vf.write("# -*- coding: utf-8 -*-\nversion = '{}'".format(version))
                 
    return version

setup(name='mechanical_components',
      version=get_version(),
      description="Design of elementary components by AI",
      long_description='',
      keywords='',
      url='',
      zip_safe=False,# To ensure static files can be loaded
      author='DessIA Technologies',
      author_email='root@dessia.tech',
      packages=['mechanical_components', 'mechanical_components.catalogs',
                'mechanical_components.optimization'],
      setup_requires=['numpy'],
      install_requires=['scipy','volmdlr','numpy', 'pandas', 'dectree',
                        'networkx', 'matplotlib'],
      data_files=[('mechanical_components/catalogs',['mechanical_components/catalogs/ferroflex.csv',
                                                     'mechanical_components/catalogs/schaeffler.json'])]
      )

