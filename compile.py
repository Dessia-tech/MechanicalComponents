#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 00:44:13 2018

@author: steven
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os
import time
import calendar
#from shutil import copyfile

import hashlib
import sys

args_delete = []
for arg in sys.argv:
    if arg.startswith('--mac='):
        mac = arg[6:]
        print('Compiling for MAC adress: {}'.format(mac))
        args_delete.append(arg)
    if arg.startswith('--exp_year='):
        year = int(arg[11:])
        args_delete.append(arg)
    if arg.startswith('--exp_month='):
        month = int(arg[12:])
        args_delete.append(arg)
for arg in args_delete:
    sys.argv.remove(arg)        

try:
    month
except NameError:
    print('Month of expiration undefined: pass with --exp_month= option')
    raise NameError

try:
    year
except NameError:
    print('Year of expiration undefined: pass with --exp_year= option')
    raise NameError
    
try:
    mac
except NameError:
    print('MAC address undefined: pass with --mac= option')
    raise NameError

expiration = int(calendar.timegm(time.struct_time((year, month, 1, 0, 0, 0 ,0, 0, 0))))
print('Expiration date: {}/{}: {}'.format(month, year, expiration))
not_before = int(time.time())

protected_files = ['mechanical_components/optimization/bearings.py',
                   'mechanical_components/optimization/common.py',
                   'mechanical_components/optimization/meshes.py',
                   'mechanical_components/optimization/wires.py',
                   ]

error_msg = 'Error, report this error to DessIA support with this traceback token: {}'.format(hashlib.sha256(str(mac).encode()).hexdigest())
protection_lines = ['valid_license = True\n',
                    't_execution = time_package.time()\n',
                    'if t_execution > {}:\n'.format(expiration), 
                    '    valid_license = False\n',
                    'if t_execution < {}:\n'.format(not_before),
                    '    valid_license = False\n',
                    'if getnode() != {}:\n'.format(mac),
                    '    valid_license = False\n',
                    'if not valid_license:\n',
                    '    print("{}")\n\n'.format(error_msg),
                    '    raise RuntimeError\n']


files_to_compile = []
for file in protected_files:
    with open(file, 'r') as f:
        # File Parsing
        new_file_lines = []        
        lines = f.readlines()
        line_index = 0
        # Inserting imports for uuid and time
        line = lines[line_index]
        while line.startswith('#!') or line.startswith('# -*-'):
            new_file_lines.append(line)
            line_index += 1
            line = lines[line_index]
            
        time_imported = False
        uuid_imported = False
        for line2 in lines:
            if 'import time as time_package' in line2:
                time_imported = True
            if 'import uuid' in line2:
                uuid_imported = True
        if not time_imported:
            new_file_lines.append('import time as time_package\n')
        if not uuid_imported:
            new_file_lines.append('from uuid import getnode\n')
        

        while line_index < len(lines):
            line = lines[line_index]
            new_file_lines.append(line)
            if line.startswith('def '):# Function
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
                    

                for protection_line in protection_lines:
                    new_file_lines.append('    {}'.format(protection_line))

                new_file_lines.append(line)


            elif line.startswith('    def ') and not line.startswith('    def __init__('):# Method of class, not init
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
                                

                for protection_line in protection_lines:
                    new_file_lines.append('        {}'.format(protection_line))

                new_file_lines.append(line)
               
            line_index+=1

                    
        new_file_name = file[:-3]+'_protected.pyx'
        files_to_compile.append(new_file_name)
        with open(new_file_name, 'w+') as nf:
            nf.writelines(new_file_lines)
                
        

ext_modules = []
for file, file_to_compile in zip(protected_files, files_to_compile):
    module = file.replace('/','.')
#    module. '_compiled'
    module = module[:-3]
    ext_modules.append(Extension(module,  [file_to_compile]))

print(ext_modules)

setup(
    name = 'powertransmission',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
            
# Remove _protected files and .c
for file in files_to_compile:
    os.remove(file)
    file = file[:-3]+'c'
    os.remove(file)