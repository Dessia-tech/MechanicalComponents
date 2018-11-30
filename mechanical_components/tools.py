#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

def StringifyDictKeys(d):
    if type(d) == list or type(d) == tuple:
        new_d = []
        for di in d:
            new_d.append(StringifyDictKeys(di))
        
    elif type(d) ==dict:
        new_d = {}
        for k,v in d.items():
            new_d[str(k)] = StringifyDictKeys(v)
    else:
        return d
    return new_d