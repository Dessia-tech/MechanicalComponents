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
        cmd = 'git describe --tags'
        try:
            version = check_output(cmd.split()).decode().strip()[:]
        except CalledProcessError:
            raise RuntimeError('Unable to get version number from git tags')
        if version[0]=='v':
            version = version[1:]
        # PEP 440 compatibility
        if '-' in version:
            future_version = version.split('-')[0].split('.')
            if 'post' in future_version[-1]:
                future_version = future_version[:-1]
            future_version[-1] = str(int(future_version[-1])+1)
            future_version = '.'.join(future_version)
            number_commits = version.split('-')[1]
            version = '{}.dev{}'.format(future_version, number_commits)
            return version

    else:
        # Extract the version from the PKG-INFO file.
        with open(join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)

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
      install_requires=['dessia-common', 'scipy', 'volmdlr', 'numpy', 'pandas', 'dectree>=0.0.4',
                        'networkx', 'matplotlib', 'genmechanics>=0.0.7',
                        'dessia_common>=0.0.3'],
      include_package_data=True,
      data_files=[('mechanical_components/catalogs',['mechanical_components/catalogs/ferroflex.csv',
                                                     'mechanical_components/catalogs/schaeffler.json'])]
      )

