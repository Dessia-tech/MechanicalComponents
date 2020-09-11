#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 11:00:11 2019

@author: Pierrem
"""

import sys as sys
import numpy as npy
import volmdlr as vm
import copy

class StiffnessMatrix:
    def __init__(self, connection, tab_pos, fix, load):
        self.connection_all = connection
        self.tab_pos = tab_pos
        self.matrix, self.sd_member = self.Matrix(connection, fix, load)
        
    def Matrix(self, connection, fix, load):
        vect = list(connection.keys())
        matrix = npy.zeros((len(vect), len(vect)))
        sd_member = npy.zeros(len(vect))
        for node1, li_data in connection.items():
            pos1 = vect.index(node1)
            for node2, (k, j, f0) in li_data.items():
                pos2 = vect.index(node2)
                matrix[pos1, pos1] = matrix[pos1, pos1] + k
                matrix[pos1, pos2] = matrix[pos1, pos2] - k
                if self.tab_pos[node1] < self.tab_pos[node2]:
                    sd_member[pos1] = sd_member[pos1] + k*j
                else:
                    sd_member[pos1] = sd_member[pos1] - k*j
                
        for i in fix:
            if i in vect:
                pos_i = vect.index(i)
                matrix[pos_i, :] = [0]*len(vect)
                matrix[pos_i, pos_i] = 1
                sd_member[pos_i] = 0
            
        for node, l in load.items():
            pos = vect.index(node)
            sd_member[pos] = sd_member[pos] + l
            
        return matrix, sd_member
    
    def Update(self, connection, fix, load):
        self.matrix, self.sd_member = self.Matrix(connection, fix, load)
    
    def Solve(self):
        x = npy.linalg.solve(self.matrix, self.sd_member)
        return x
    
    def AnalyzeLoad(self):
        vect = list(self.connection_init.keys())
        for node1, li_data in self.connection_init.items():
            pos1 = vect.index(node1)
            for node2, (k, j, f0) in li_data.items():
                pos2 = vect.index(node2)
                if tab_pos[node1] < tab_pos[node2]:
                    load = k*(x[pos1] - x[pos2]) - k*j
                    dist = x[pos1] - x[pos2]
                else:
                    load = -k*(x[pos1] - x[pos2]) - k*j
                    dist = -(x[pos1] - x[pos2])
    
    