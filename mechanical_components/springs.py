#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 17:02:33 2018

@author: jezequel
"""
import math
import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

from scipy.optimize import minimize

class Spring():
    def __init__(self, displacement = 0.005, D = 0.0050, d = 0.0010, n = 5,
                 F0 = 10, F1 = 100,
                 young_modulus = 210*10**9, poisson_ratio = 0.33):
        self.displacement = displacement
        self.D = D
        self.d = d
        self.F0 = F0
        self.F1 = F1
        self.young_modulus = young_modulus
        self.poisson_ratio = poisson_ratio
        self.n = n
        
    def Update(self, values):
        for key,value in values.items():
            str_split = key.split('.')
            if len(str_split) == 2:
                setattr(getattr(self, str_split[0]), str_split[1], value)
            elif len(str_split) == 1:
                setattr(self,key,value)
                
    def Lengths(self):
        l0 = self.displacement*self.F0/(self.F1-self.F0)
        l1 = l0 + self.displacement
        
        return l0, l1
                
    def TargetStiffness(self):
        l0, l1 = self.Lengths()
        
        k = self.F1/l1
        
        return k
    
    def Stiffness(self):
        E = self.young_modulus
        nu = self.poisson_ratio
        
        G = E/(2*(1+nu))
        
        k = G*self.d**4/(8*self.n*self.D**3)
        return k
    
    def EngagedResultingForce(self):
        l0, l1 = self.Lengths()
        F = -self.Stiffness()*l1
        
        return F
    
    def LoseResultingForce(self):
        l0, l1 = self.Lengths()
        F = -self.Stiffness()*l0
        
        return F
        
class SpringOptimizer():
    def __init__(self, spring, specs):
        self.spring = spring
        self.specs = specs
        self.bounds = []
        self.attributes = []
        self.fixed_values = {}
        
        for k,v in self.specs.items():
            tv = type(v)
            if tv == tuple:
                self.attributes.append(k)
                self.bounds.append(v)
            else:
                self.fixed_values[k] = v

        self.n = len(self.attributes)
        self.spring.Update(self.fixed_values)
        
    def Optimize(self):
        def Objective(xa):
            print('xa',xa)
            values = {}
            for xai, attribute, bounds in zip(xa, self.attributes, self.bounds):                    
                values[attribute] = bounds[0] + (bounds[1] - bounds[0])*xai
            self.spring.Update(values)
            print('fun',(self.spring.Stiffness() - self.spring.TargetStiffness())**2)
            return (1e-3*abs(self.spring.Stiffness() - self.spring.TargetStiffness()))
        
        def DiameterRatioConstraintMin(xa):
            values = {}
            for xai, attribute, bounds in zip(xa, self.attributes, self.bounds):                    
                values[attribute] = bounds[0] + (bounds[1] - bounds[0])*xai
            self.spring.Update(values)
            print('min',self.spring.D/self.spring.d-5)
            return self.spring.D/self.spring.d - 5
        
        def DiameterRatioConstraintMax(xa):
            values = {}
            for xai, attribute, bounds in zip(xa, self.attributes, self.bounds):                    
                values[attribute] = bounds[0] + (bounds[1] - bounds[0])*xai
            self.spring.Update(values)
            print('max',13 - self.spring.D/self.spring.d)
            return 13 - self.spring.D/self.spring.d
        
        def ParametersBoundsMin(xa):
            return min(xa)
        
        def ParametersBoundsMax(xa):
            return 1-max(xa)
        
#        def DiameterConstructionBoundsMin(xa):
#            return (self.spring.D+self.spring.d/2)/(self.spring.D-self.spring.d/2) - 1.4
#        
#        def DiameterConstructionBoundsMax(xa):
#            return 2 - (self.spring.D+self.spring.d/2)/(self.spring.D-self.spring.d/2)
        
        fun_constraints = [{'type' : 'ineq', 'fun' : DiameterRatioConstraintMin},
                           {'type' : 'ineq', 'fun' : DiameterRatioConstraintMax}]
#                           {'type' : 'ineq', 'fun' : ParametersBoundsMin},
#                           {'type' : 'ineq', 'fun' : ParametersBoundsMax}]#,
#                           {'type' : 'ineq', 'fun' : DiameterConstructionBoundsMin},
#                           {'type' : 'ineq', 'fun' : DiameterConstructionBoundsMax}]
        
        for i in range(1000):
            xra0 = npy.random.random(self.n)
            if DiameterRatioConstraintMin(xra0)>0:
                if DiameterRatioConstraintMax(xra0)>0:
                    print('valid')
                    break
                
        res = minimize(Objective, xra0, constraints = fun_constraints, bounds = [(0., 1.)]*self.n)
#        res = fmin_slsqp(Objective, xra0, ieqcons = [DiameterRatioConstraintMin,DiameterRatioConstraintMax], bounds = [(0., 1.)]*self.n,epsilon=0.1)
#        print(DiameterRatioConstraintMax())
#        res = minimize(Objective, xra0, bounds = [(0, 1)]*self.n)
        return res