#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""
import os
import sys
import re
from subprocess import CalledProcessError, check_output
from setuptools import setup
from Cython.Distutils import build_ext
from distutils.extension import Extension
import time
from datetime import timedelta
import calendar
import shutil

from distutils import log as logger

import hashlib
import netifaces

import wheel.bdist_wheel
from wheel.wheelfile import WheelFile


protected_files = ['mechanical_components/optimization/bearings_protected.py',
                   'mechanical_components/optimization/common.py',
                   'mechanical_components/optimization/meshes_protected.py',
                   'mechanical_components/optimization/wires_protected.py',
                   'mechanical_components/shafts.py',
                   #'mechanical_components/shafts_assembly.py'
                   ]


ext_modules = []
for file in protected_files:
    module = file.replace('/', '.')
    module = module[:-3]
    file_to_compile = file[:-3] + '_protected.pyx'
    ext_modules.append(Extension(module,  [file_to_compile]))
    
print('ext_modules', ext_modules)


class ClientWheelDist(wheel.bdist_wheel.bdist_wheel):
    description = 'Creating client distribution with compiled packages and license'
    user_options = [
        # The format is (long option, short option, description).
        ('exp-year=', None, 'Year of license expiration'),
        ('exp-month=', None, 'Month of year of license expiration'),
        ('exp-day=', None, 'Day of month of license expiration'),
        ('getnodes=', None, 'UUID given by getnode (comma-separated list)'),
        ('macs=', None, 'MACS of machine (comma-separated list)'),
        ('detect-macs', None, 'using this machine macs'),
        ('delete-pyx=', None, 'delete compilation files')
    ]

    addresses_determination_lines = ['import netifaces\n',
                                     'addrs = []\n',
                                     'for i in netifaces.interfaces():\n',
                                     "    if i != 'lo':\n",
                                     "        for k, v in netifaces.ifaddresses(i).items():\n",
                                     "            for v2 in v:\n",
                                     "                if 'addr' in v2:\n",
                                     "                    a = v2['addr']\n",
                                     "                    if len(a) == 17:\n",
                                     "                        addrs.append(a.replace(':',''))\n",
                                     "addrs = set(addrs)\n\n"
                                     ]

    def get_machine_macs(self):
        addrs = []
        for i in netifaces.interfaces():
            if i != 'lo':
                for k, v in netifaces.ifaddresses(i).items():
                    for v2 in v:
                        if 'addr' in v2:
                            a = v2['addr']
                            if len(a) == 17:
                                addrs.append(a.replace(':', ''))
        return list(set(addrs))

    def initialize_options(self):
        """Set default values for options."""

        wheel.bdist_wheel.bdist_wheel.initialize_options(self)

        # Each user option must be listed here with their default value.
        self.exp_year = 2000
        self.exp_month = 1
        self.exp_day = 1
        self.formats = 'gztar'
        self.getnodes = None
        self.macs = None
        self.detect_macs = None
        self.dist_dir = 'client_dist'
        self.delete_pyx = True

    def finalize_options(self):
        """Post-process options."""

        wheel.bdist_wheel.bdist_wheel.finalize_options(self)

        self.detect_macs = self.detect_macs is not None
        if not self.detect_macs:

            macs = []
            if self.macs is not None:
                self.macs = self.macs.replace('[', '')
                self.macs = self.macs.replace(']', '')
                self.macs = self.macs.replace(':', '')

                for mac in self.macs.split(','):
                    if len(mac) != 12:
                        raise ValueError(
                            'A mac address must be 12 digits long, got: {} instead'.format(
                                mac))
                    macs.append(mac)
                self.macs = macs
                print('\nCompiling for macs: {}'.format(self.macs))

            getnodes = []
            if self.getnodes is not None:
                for getnode in self.getnodes.split(','):
                    getnodes.append(getnode)
                self.getnodes = getnodes
                print('\nCompiling for getnodes: {}'.format(
                    [str(g for g in self.getnodes)]))

        else:
            self.macs = self.get_machine_macs()
            print('Using detected mac of this machine: {}'.format(self.macs))

        if not self.detect_macs and self.getnodes is None and self.macs is None:
            print('Warning: no macs protection')
            # raise ValueError(
            #     'Define either a mac or a getnode to protect the code or use detect-macs option')

        if self.exp_year is not None:
            self.exp_year = int(self.exp_year)
        if self.exp_month is not None:
            self.exp_month = int(self.exp_month)
        if self.exp_day is not None:
            self.exp_day = int(self.exp_day)

        self.expiration = int(calendar.timegm(time.struct_time((self.exp_year,
                                                                self.exp_month,
                                                                self.exp_day,
                                                                0, 0, 0, 0, 0,
                                                                0))))
        self.not_before = int(time.time())
        number_days = timedelta(seconds=self.expiration - self.not_before).days
        if number_days < 0:
            raise ValueError('Negative number of days')

        print('\nExpiration date: {}/{}/{}: duration of {} days'.format(
            self.exp_day,
            self.exp_month,
            self.exp_year,
            number_days))

        formats = []
        for s in self.formats.split(','):
            if not s in ['zip', 'tar', 'gztar', 'bztar']:
                raise NotImplementedError(
                    'Unsupported file type: {}'.format(s))
            formats.append(s)
        self.formats = formats

        self.delete_pyx = self.delete_pyx == 'True'

    def protection_lines(self):
        protection_lines = []
        if self.macs is not None:
            physical_token = hashlib.sha256(
                str(self.macs[0]).encode()).hexdigest()
        elif self.getnodes is not None:
            physical_token = hashlib.sha256(
                str(self.getnodes[0]).encode()).hexdigest()
        else:
            physical_token=''

        error_msg_time_before = 'Invalid license Please report this error to DessIA support with this traceback token: TB{}'.format(
            physical_token)
        error_msg_time_after = 'Invalid license Please report this error to DessIA support with this traceback token: TA{}'.format(
            physical_token)
        error_msg_mac = 'Invalid license. Please report this error to DessIA support with this traceback token: M{}'.format(
            physical_token)
        protection_lines = ['valid_license = True\n',
                            't_execution = time_package.time()\n',
                            'if t_execution > {}:\n'.format(self.expiration),
                            '    print("{}")\n'.format(error_msg_time_after),
                            '    raise RuntimeError\n\n',
                            'if t_execution < {}:\n'.format(self.not_before),
                            '    print("{}")\n'.format(error_msg_time_before),
                            '    raise RuntimeError\n\n']

        if self.macs is not None:
            protection_lines.extend(
                ['if addrs.isdisjoint(set({})):\n'.format(self.macs),
                 '    print("{}")\n'.format(error_msg_mac),
                 '    raise RuntimeError\n\n'])
        elif self.getnodes is not None:
            protection_lines.extend(
                ['if getnode() not in {}:\n'.format(self.getnodes),
                 '    print("{}")\n'.format(error_msg_mac),
                 '    raise RuntimeError\n\n'])

        return protection_lines

    def write_pyx_files(self):
        self.files_to_compile = []
        for file in protected_files:
            with open(file, 'r', encoding="utf-8") as f:
                # File Parsing
                new_file_lines = []
                lines = f.readlines()
                line_index = 0
                # Inserting imports for uuid and time
                line = lines[line_index]
                while line.startswith('#!') or line.startswith('# -*-') or (
                        'cython' in line):
                    new_file_lines.append(line)
                    line_index += 1
                    line = lines[line_index]

                new_file_lines.append('import time as time_package\n')
                if self.getnodes is not None:
                    new_file_lines.append('from uuid import getnode\n')
                if self.macs is not None:
                    new_file_lines.extend(self.addresses_determination_lines)

                while line_index < len(lines):
                    line = lines[line_index]
                    new_file_lines.append(line)
                    if line.startswith('def '):  # Function
                        # counting parenthesis
                        op = line.count('(')
                        cp = line.count(')')
                        while op != cp:
                            line_index += 1
                            line = lines[line_index]
                            new_file_lines.append(line)
                            op += line.count('(')
                            cp += line.count(')')
                        # list of args is finished.
                        # Now trying to see if lines are docstrings
                        line_index += 1
                        line = lines[line_index]

                        if line.startswith('    """'):
                            new_file_lines.append(line)
                            line_index += 1
                            line = lines[line_index]

                            while not line.startswith('    """'):
                                new_file_lines.append(line)
                                line_index += 1
                                line = lines[line_index]
                            new_file_lines.append(line)
                            line_index += 1
                            line = lines[line_index]

                        for protection_line in self.protection_lines():
                            new_file_lines.append(
                                '    {}'.format(protection_line))

                        new_file_lines.append(line)


                    elif line.startswith('    def ') and not line.startswith(
                            '    def __init__('):  # Method of class, not init
                        # counting parenthesis
                        op = line.count('(')
                        cp = line.count(')')
                        while op != cp:
                            line_index += 1
                            line = lines[line_index]
                            op += line.count('(')
                            cp += line.count(')')
                            new_file_lines.append(line)
                        # list of args is finished.
                        # Now trying to see if lines are docstrings
                        line_index += 1
                        line = lines[line_index]

                        if line.startswith('        """'):
                            new_file_lines.append(line)
                            line_index += 1
                            line = lines[line_index]
                            while not line.startswith('        """'):
                                new_file_lines.append(line)
                                line_index += 1
                                line = lines[line_index]
                            new_file_lines.append(line)
                            line_index += 1
                            line = lines[line_index]

                        for protection_line in self.protection_lines():
                            new_file_lines.append(
                                '        {}'.format(protection_line))

                        new_file_lines.append(line)

                    line_index += 1

                new_file_name = file[:-3] + '_protected.pyx'
                self.files_to_compile.append(new_file_name)
                with open(new_file_name, 'w+', encoding="utf-8") as nf:
                    nf.writelines(new_file_lines)

    def delete_compilation_files(self):
        # Remove _protected files and .c
        for file in self.files_to_compile:
            if os.path.exists(file):
                os.remove(file)

            file = file[:-3] + 'c'
            if os.path.exists(file):
                os.remove(file)

    def run(self):
        self.write_pyx_files()

        build_scripts = self.reinitialize_command('build_scripts')
        build_scripts.executable = 'python'
        build_scripts.force = True

        build_ext = self.reinitialize_command('build_ext')
        build_ext.inplace = False

        if not self.skip_build:
            self.run_command('build')

        install = self.reinitialize_command('install',
                                            reinit_subcommands=True)
        install.root = self.bdist_dir
        install.compile = False
        install.skip_build = self.skip_build
        install.warn_dir = False

        # A wheel without setuptools scripts is more cross-platform.
        # Use the (undocumented) `no_ep` option to setuptools'
        # install_scripts command to avoid creating entry point scripts.
        install_scripts = self.reinitialize_command('install_scripts')
        install_scripts.no_ep = True

        # Use a custom scheme for the archive, because we have to decide
        # at installation time which scheme to use.
        for key in ('headers', 'scripts', 'data', 'purelib', 'platlib'):
            setattr(install,
                    'install_' + key,
                    os.path.join(self.data_dir, key))

        basedir_observed = ''

        if os.name == 'nt':
            # win32 barfs if any of these are ''; could be '.'?
            # (distutils.command.install:change_roots bug)
            basedir_observed = os.path.normpath(
                os.path.join(self.data_dir, '..'))
            self.install_libbase = self.install_lib = basedir_observed

        setattr(install,
                'install_purelib' if self.root_is_pure else 'install_platlib',
                basedir_observed)

        logger.info("installing to %s", self.bdist_dir)

        self.run_command('install')

        impl_tag, abi_tag, plat_tag = self.get_tag()
        archive_basename = "{}-{}-{}-{}".format(self.wheel_dist_name, impl_tag,
                                                abi_tag, plat_tag)
        if not self.relative:
            archive_root = self.bdist_dir
        else:
            archive_root = os.path.join(
                self.bdist_dir,
                self._ensure_relative(install.install_base))

        self.set_undefined_options('install_egg_info',
                                   ('target', 'egginfo_dir'))
        distinfo_dirname = '{}-{}.dist-info'.format(
            wheel.bdist_wheel.safer_name(self.distribution.get_name()),
            wheel.bdist_wheel.safer_version(self.distribution.get_version()))
        distinfo_dir = os.path.join(self.bdist_dir, distinfo_dirname)

        package_name = self.distribution.get_name()
        distdir = os.path.join(self.bdist_dir, package_name)

        for root, dirs, files in os.walk(distdir):
            for file in files:
                file_path = os.path.join(root, file)
                rel_fp = os.path.join(package_name,
                                      os.path.relpath(file_path, distdir))
                if rel_fp in protected_files:
                    os.remove(file_path)

        self.egg2dist(self.egginfo_dir, distinfo_dir)

        self.write_wheelfile(distinfo_dir)

        # Make the archive
        if not os.path.exists(self.dist_dir):
            os.makedirs(self.dist_dir)

        wheel_path = os.path.join(self.dist_dir, archive_basename + '.whl')
        with WheelFile(wheel_path, 'w', self.compression) as wf:
            wf.write_files(archive_root)

        # Add to 'Distribution.dist_files' so that the "upload" command works
        getattr(self.distribution, 'dist_files', []).append(
            ('bdist_wheel',
             '{}.{}'.format(*sys.version_info[:2]),  # like 3.7
             wheel_path))

        if self.delete_pyx:
            self.delete_compilation_files()

        if not self.keep_temp:
            logger.info('removing %s', self.bdist_dir)
            if not self.dry_run:
                shutil.rmtree(self.bdist_dir,
                              onerror=wheel.bdist_wheel.remove_readonly)


tag_re = re.compile(r'\btag: %s([0-9][^,]*)\b')
version_re = re.compile('^Version: (.+)$', re.M)
    
def version_from_git_describe(version):
    if version[0]=='v':
            version = version[1:]

    # PEP 440 compatibility
    number_commits_ahead = 0
    if '-' in version:
        version, number_commits_ahead, commit_hash = version.split('-')
        number_commits_ahead = int(number_commits_ahead)

    print('number_commits_ahead', number_commits_ahead)

    split_versions = version.split('.')
    if 'post' in split_versions[-1]:
        suffix = split_versions[-1]
        split_versions = split_versions[:-1]
    else:
        suffix = None

    for pre_release_segment in ['a', 'b', 'rc']:
        if pre_release_segment in split_versions[-1]:
            if number_commits_ahead > 0:
                split_versions[-1] = str(split_versions[-1].split(pre_release_segment)[0])
                if len(split_versions) == 2:
                    split_versions.append('0')
                if len(split_versions) == 1:
                    split_versions.extend(['0', '0'])

                split_versions[-1] = str(int(split_versions[-1])+1)
                future_version = '.'.join(split_versions)
                return '{}.dev{}'.format(future_version, number_commits_ahead)
            else:
                return '.'.join(split_versions)

    if number_commits_ahead > 0:
        if len(split_versions) == 2:
            split_versions.append('0')
        if len(split_versions) == 1:
            split_versions.extend(['0', '0'])
        split_versions[-1] = str(int(split_versions[-1])+1)
        split_versions = '.'.join(split_versions)
        return '{}.dev{}'.format(split_versions, number_commits_ahead)
    else:
        if suffix is not None:
            split_versions.append(suffix)

        return '.'.join(split_versions)


def readme():
    with open('README.md') as f:
        return f.read()
    

def get_version():
    # Return the version if it has been injected into the file by git-archive
    version = tag_re.search('$Format:%D$')
    if version:
        return version.group(1)

    d = os.path.dirname(__file__)

    if os.path.isdir(os.path.join(d, '.git')):
        cmd = 'git describe --tags'
        try:
            version = check_output(cmd.split()).decode().strip()[:]

        except CalledProcessError:
            raise RuntimeError('Unable to get version number from git tags')

        return version_from_git_describe(version)
    else:
        # Extract the version from the PKG-INFO file.
        with open(os.path.join(d, 'PKG-INFO')) as f:
            version = version_re.search(f.read()).group(1)

    # print('version', version)
    return version



setup(name='mechanical_components',
      version=get_version(),
      description="Design of elementary components by AI",
      long_description='',
      keywords='',
      url='',
      zip_safe=False,# To ensure static files can be loaded
      author='DessiA Technologies',
      author_email='root@dessia.tech',
      packages=['mechanical_components', 'mechanical_components.catalogs',
                'mechanical_components.models',
                'mechanical_components.optimization'],
      setup_requires=['numpy'],
      install_requires=['dessia-common>=0.4.0', 'scipy', 'volmdlr>=0.2.1', 'numpy', 'pandas', 'dectree>=0.0.4',
                        'networkx', 'matplotlib', 'genmechanics>=0.1.5'
                        ],
      include_package_data=True,
      data_files=[('mechanical_components/catalogs',['mechanical_components/catalogs/ferroflex.csv',
                                                     'mechanical_components/catalogs/schaeffler.json'])],
      cmdclass = {'build_ext': build_ext,
                  'cdist_wheel': ClientWheelDist},
      ext_modules = ext_modules 
      
      )
