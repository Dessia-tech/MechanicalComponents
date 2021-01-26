#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 17:02:33 2018

@author: jezequel
"""
import math
import numpy as npy

from dessia_common.core import DessiaObject

import matplotlib.pyplot as plt
import matplotlib.colors as colors

import volmdlr as vm
import volmdlr.primitives3d as primitives3d
import volmdlr.primitives2d as primitives2d

import pandas as pd
from pandas.plotting import scatter_matrix

from copy import copy

from scipy.optimize import minimize

import pkg_resources


class Material(DessiaObject):
    def __init__(self, volumic_mass, young_modulus, poisson_ratio, Rm, d_min = 0.12*10**-3, d_max = 12*10**-3, cost_index = 10, name = ''):
        self.volumic_mass = volumic_mass
        self.young_modulus = young_modulus
        self.poisson_ratio = poisson_ratio
        self.Rm = Rm
        self.d_min = d_min
        self.d_max = d_max
        
        self.cost_index = cost_index
        
        self.G = self.ShearModulus()
        self.tau_max = self.MaxShearStress()
        
        self.name = name
    def ShearModulus(self):
        E = self.young_modulus
        nu = self.poisson_ratio
        
        G = E/(2*(1+nu))
        
        return G
    
    def MaxShearStress(self):
        tau_max = 0.56*self.Rm
        return tau_max
    
    def Dict(self):
        return self.__dict__
        
steel1 = Material(7850, 200*10**9, 0.3, 750*10**6, 0.8*10**-3, 12*10**-3, 1, 'XC60')
steel2 = Material(7850, 200*10**9, 0.3, 1000*10**6, 3*10**-3, 12*10**-3, 2.5, 'XC95')
steel3 = Material(7850, 196*10**9, 0.3, 780*10**6, 0.8*10**-3, 12*10**-3, 3, '50CV4')
steel4 = Material(7850, 195*10**9, 0.3, 780*10**6, cost_index = 9, name = 'Z15CN17-03')
copper_tin_alloy = Material(8730, 115*10**9, 0.3, 950*10**6, cost_index = 17, name = 'CuSn6 R950')
copper_zinc_alloy = Material(8400, 110*10**9, 0.3, 700*10**6, cost_index = 18, name = 'CuZn36 R700')
#copper_beryllium_alloy = Material(8800, 120*10**9, 0.3, name = 'CuBe2')
#copper_cobalt_beryllium_alloy = Material(8800, 130*10**9, 0.3, name = 'CuCo2Be')
materials = [steel1,
             steel2,
             steel3,
             steel4,
             copper_tin_alloy,
             copper_zinc_alloy]

n_spires = npy.linspace(2.5, 9.5, 8)
diameters_mm = [0.12, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
                0.95, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.7, 2.8, 3, 3.2, 3.4, 3.5, 3.6, 3.8,
                4, 4.2, 4.5, 4.70, 4.8, 5, 5.3, 5.5, 5.6,6, 6.3, 6.5, 6.70, 7, 7.5,
                8, 8.5, 9, 9.5, 10, 11, 12] # , 13, 14]

diameters_m = [d*10**-3 for d in diameters_mm]

class Spring(DessiaObject):
    def __init__(self, D = 0.0050, d = 0.0010, n = 5, l0 = 0.010, material = steel1):
        self.D = D
        self.d = d
        self.young_modulus = material.young_modulus
        self.poisson_ratio = material.poisson_ratio
        self.n = n
        self.l0 = l0
        self.p = l0/n
        
        if d <= 10*10**-3:
            self.coil = 'cold'
            self.nt = n+2
            self.Sa = 1.5*n*(0.0015*D**2/d + 0.001*d) # Coeff : Dynamic load
        else:
            self.coil = 'hot'
            self.nt = n+1.5
            self.Sa = 2*0.002*n*(D + d) # Coeff : Dynamic load
            
        self.lc = self.nt*self.d
        self.l_min = self.lc + self.Sa
        
        self.material = material
        self.matching_products = []
        
        
        self.k = self.Stiffness()
        self.w = self.SpringIndex()
        self.i = self.InclinationAngle()
        self.lw = self.WireLength()
        self.m = self.Mass()
        self.cost = self.Cost()
        
        self.contour = self.Contour()

    def __getstate__(self):
        d=self.__dict__.copy()
        del d['contour']
        return d
    
    def __setstate__(self,d):
        self.__dict__=d
        self.contour = self.Contour()
        

    def Update(self, values):
        for key,value in values.items():
            str_split = key.split('.')
            if len(str_split) == 2:
                setattr(getattr(self, str_split[0]), str_split[1], value)
            elif len(str_split) == 1:
                setattr(self,key,value)
                
        self.p = self.l0/self.n
        
        self.k = self.Stiffness()
        self.w = self.SpringIndex()
        self.i = self.InclinationAngle()
        self.lw = self.WireLength()
        self.m = self.Mass()
        self.cost = self.Cost()
       
        self.contour = self.Contour()
#        self.volume = self.Volume()
        
    def SpringPosition(self, position):
        self.pos_x, self.pos_y = position
    
    def Stiffness(self):
        E = self.young_modulus
        nu = self.poisson_ratio
        
        G = E/(2*(1+nu))
        
        k = G*self.d**4/(8*self.n*self.D**3)
        return k
    
    def ResultingForce(self, l):
        F = self.Stiffness()*l
        
        return F
    
    def Length(self, F):
        l = self.l0 - F/self.Stiffness()
        
        return l
    
    def SpringIndex(self):
        w = self.D/self.d
        
        return w
    
    def InclinationAngle(self):
        i = math.atan(self.p/(math.pi*self.D))
        
        return i
    
    def WireLength(self):
        lw = math.pi*self.D*(2+self.n)/math.cos(self.InclinationAngle())
        return lw
    
    def Mass(self):
        m = self.WireLength()*self.material.volumic_mass*math.pi*self.d**2/4
        return m
    
    def Cost(self):
        ci = self.WireLength()*self.d*self.material.cost_index
        
        return ci
    
    def Custom2Catalog(self, accuracy_rate = 0.05):
        for cat_name, catalog in catalogs.items():
            specs = {'R' : (self.Stiffness()*(1 - accuracy_rate), self.Stiffness()*(1 + accuracy_rate)),
                     'L0' : (self.Length(0)*(1 - accuracy_rate), self.Length(0)*(1 + accuracy_rate))}
            panda_indices = catalog.FindSprings(specs)
    
            products = [Product(cat_name, ind) for ind in panda_indices.get_values()]
            self.matching_products.extend(products)
    
    def Contour(self):
        p0 = vm.Point2D(self.D, 0)
        
        l1 = vm.wires.Circle2D(p0, self.d/2)
        
        return vm.wires.Contour2D([l1])
    
    def Volume(self, F=0., position=vm.O3D, x=vm.X3D, z=vm.Z3D):
#        p_spring_coord = (self.pos_x, self.pos_y, 0)
#        p_spring = vm.Point3D(p_spring_coord)
        
#        xp_coord = (self.pos_x/(math.sqrt(self.pos_x**2 + self.pos_y**2)),
#                    self.pos_y/(math.sqrt(self.pos_x**2 + self.pos_y**2)),
#                    0)
#        p_plan = p_spring.Translation(((self.D/(2*math.sqrt(2)))*xp_coord[0],
#                                       (self.D/(2*math.sqrt(2)))*xp_coord[1], 0))
#        p0_coord = ((self.D/2 + 1), (self.D/2 + 1), 0)
#        zp_coord = (0, 0, 1)
        
#        p0 = vm.Point3D(p0_coord)
#        pc = p0.Translation((self.d/2, 0, 0))
        
#        xp = vm.Vector3D(xp_coord)
#        zp = vm.Vector3D(zp_coord)
        
        l = self.Length(F)
        p = l/self.n
        
        primitives = []
        # volume = primitives3d.HelicalExtrudedProfile(position, x, z, position,
        #                                              l*z, p, self.contour, name = 'spring')
        # primitives.append(volume)
        
        return primitives
    
    def CADExport(self,name,python_path='python',freecad_path='/usr/lib/freecad/lib/',export_types=['stl','fcstd']):
        volumes = []
        volumes.extend(self.Volume(0))
        
        model = vm.VolumeModel(volumes)
        resp = model.FreeCADExport(python_path,name,freecad_path,export_types)
        
        return resp
    
    def Dict(self):
        d=self.__dict__.copy()
        d['material']=self.material.Dict()
        matching_products = []
        for matching_product in self.matching_products:
            matching_products.append(matching_product.Dict())
        d['matching_products'] = matching_products
        
        del d['contour']
        
        return d
    
class SpringAssembly(DessiaObject):
    def __init__(self, springs, geometry):
        self.springs = springs
        self.n_springs = len(springs)
        self.geometry = geometry
        
        self.k = self.Stiffness()
        self.m = self.Mass()
        self.cost = self.Cost()
        self.l0 = self.FreeLength()
        self.PositionSprings()
        
        self.matching_product_assemblies = []
#        self.matching_product_assemblies_strict = []
        
    def Stiffness(self):
        k = sum([spring.Stiffness() for spring in self.springs])
        
        return k
    
    def Mass(self):
        m = sum([spring.Mass() for spring in self.springs])
        
        return m
    
    def Cost(self):
        c = sum([spring.Cost() for spring in self.springs])
        
        return c
    
    def FreeLength(self):
        # A modifier
        spring = self.springs[0]
        l0 = spring.l0
        
        return l0
    
    def PositionSprings(self):
        if self.geometry['pattern'] == 'circular':
            radius = self.geometry['radius']
            angle = self.geometry['angle']
            
            [spring.SpringPosition((radius*math.cos(i*angle), radius*math.sin(i*angle))) for i, spring in enumerate(self.springs)]
            
        elif self.geometry['pattern'] == 'shaft mounted':
            [spring.SpringPosition((0, 0)) for i, spring in enumerate(self.springs)]
    
    def CADExport(self,name,python_path='python',freecad_path='/usr/lib/freecad/lib/',
                  export_types=['stl','fcstd']):
        volumes = []
        for spring in self.springs:
            volumes.extend(spring.Volume(0))
            
        model = vm.VolumeModel(volumes)
        resp = model.FreeCADExport(python_path,name,freecad_path,export_types)
        return resp
    
    def Custom2Catalog(self, accuracy_rate = 0.05):
        products = []
        product_indices = []
        for spring in self.springs:
            for product in spring.matching_products:
                if product.product_index not in product_indices:
                    product_indices.append(product.product_index)
                    products.append(product)
                    
        for product in products:
            self.matching_product_assemblies.append(ProductAssembly([product]*self.n_springs, self.geometry))
            
#        for product_assembly in self.matching_product_assemblies:
#            msa = product_assembly.Instantiate(self.geometry)
#            if self.k > msa.k*(1 - accuracy_rate) and self.k < msa.k*(1 + accuracy_rate)\
#            and self.l0 > msa.l0*(1 - accuracy_rate) and self.l0 < msa.l0*(1 + accuracy_rate):
#                self.matching_product_assemblies_strict.append(product_assembly)
            
    def Dict(self):
        d=self.__dict__.copy()
        springs = []
        matching_product_assemblies = []
        for spring in self.springs:
            springs.append(spring.Dict())
        for matching_product_assembly in self.matching_product_assemblies:
            matching_product_assemblies.append(matching_product_assembly.Dict())
            
        d['springs'] = springs
        d['matching_product_assemblies'] = matching_product_assemblies
        return d
        
class SpringOptimizer(DessiaObject):
    def __init__(self, spring, specs, F2):
        self.spring = spring
        self.specs = specs
        self.F2 = F2
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
        
    def UpdateObject(self, xa):
        values = {}
        for xai, attribute, bounds in zip(xa, self.attributes, self.bounds):                    
            values[attribute] = bounds[0] + (bounds[1] - bounds[0])*xai
        self.spring.Update(values)
        
    def Optimize(self):
        def Objective(xa):
            self.UpdateObject(xa)
            return self.spring.Mass()
                
        def MinLengthConstraint(xa):
            self.UpdateObject(xa)
            return self.spring.l0 - self.F2/self.spring.Stiffness() - self.spring.l_min
                
        def MaxLengthConstraint(xa):
            self.UpdateObject(xa)
            return 0.1*math.pi*self.spring.n*self.spring.D - self.spring.l0
        
        def MinSpringIndexConstraint(xa):
            self.UpdateObject(xa)
            return self.spring.SpringIndex() - 5
        
        def MaxSpringIndexConstraint(xa):
            self.UpdateObject(xa)
            return 13 - self.spring.SpringIndex()
        
        def MinDiameterConstraint(xa):
            self.UpdateObject(xa)
            return (self.spring.D + self.spring.d)/(self.spring.D - self.spring.d) - 1.4
        
        def MaxDiameterConstraint(xa):
            self.UpdateObject(xa)
            return 2 - (self.spring.D + self.spring.d)/(self.spring.D - self.spring.d)
        
        fun_constraints = [{'type' : 'ineq', 'fun' : MinLengthConstraint},
                           {'type' : 'ineq', 'fun' : MaxLengthConstraint},
                           {'type' : 'ineq', 'fun' : MinSpringIndexConstraint},
                           {'type' : 'ineq', 'fun' : MaxSpringIndexConstraint},
                           {'type' : 'ineq', 'fun' : MinDiameterConstraint},
                           {'type' : 'ineq', 'fun' : MaxDiameterConstraint}]

#        for i in range(1000):
        xra0 = npy.random.random(self.n)
#            if DiameterRatioConstraintMin(xra0)>0:
#                if DiameterRatioConstraintMax(xra0)>0:
#                    print('valid')
#                    break
                
        res = minimize(Objective, xra0, constraints = fun_constraints, bounds = [(0., 1.)]*self.n)
        return res
        
    
class SpringDiscreteOptimizer(DessiaObject):
    def __init__(self, F1, F2, stroke, d, n, material = steel1):
        self.F1 = F1
        self.F2 = F2
        self.stroke = stroke
        self.d = d
        self.n = n
        
        self.E = material.young_modulus
        self.nu = material.poisson_ratio
        
        self.G = self.RigidityModulus()
        self.k = self.TargetStiffness()
        self.D = self.OutsideDiameter()
        self.w = self.SpringIndex()
        self.sc_factor = self.StressCorrectionFactor()
        self.tau_k = self.ShearStress()
        
        if d <= 10*10**-3:
            self.coil = 'cold'
            self.nt = n+2
            self.Sa = 1.5*n*(0.0015*self.D**2/d + 0.001*d) # Coeff : Dynamic load
        else:
            self.coil = 'hot'
            self.nt = n+1.5
            self.Sa = 2*0.002*n*(self.D + d) # Coeff : Dynamic load
            
        self.lc = self.nt*d
        self.l_min = self.lc + self.Sa
            
        self.l0 = self.MinimumFreeLength()
        self.l1 = self.Length(F1)
        
    def RigidityModulus(self):
        G = self.E/(2*(1 + self.nu))
        
        return G
        
    def TargetStiffness(self):
        k = (self.F2 - self.F1)/self.stroke
        
        return k
        
    def OutsideDiameter(self):
        D = ((self.G*(self.d)**4)/(8*self.n*self.k))**(1/3)
        
        return D
    
    def SpringIndex(self):
        w = self.D/self.d
        
        return w
    
    def StressCorrectionFactor(self):
        sc_factor = (self.w + 0.5)/(self.w - 0.75)
        
        return sc_factor
    
    def ShearStress(self):
        tau_k = self.sc_factor*(8*self.D*self.F2)/(math.pi*self.d**3)
        
        return tau_k
    
    def MinimumFreeLength(self):
        l0 = self.F2/self.k + self.l_min
        
        return l0
    
    def Length(self, F):
        l = self.l0 - F/self.k
        
        return l
    
class SpringAssemblyOptimizer(DessiaObject):
    def __init__(self, F1, F2, stroke, n_springs, r1, r2, l1_max, pattern = 'circular'):
        self.F1 = F1
        self.F2 = F2
        self.stroke = stroke
        self.n_springs = n_springs
        self.pattern = pattern
        
        self.target_k = self.TargetStiffness()
        
        self.assemblies = []
        
        if pattern == 'circular':
            for i in n_springs:
                F1eq = F1/i
                F2eq = F2/i
                angle = 2*math.pi/i
                for d in diameters_m:
                    for n in n_spires:
                        for material in materials:
                            sdo = SpringDiscreteOptimizer(F1eq, F2eq, stroke, d, n, material)
                            if sdo.tau_k < material.tau_max\
                            and d >= material.d_min and d <= material.d_max\
                            and sdo.D/d > 5 and sdo.D/d < 13\
                            and (sdo.D+d)/(sdo.D-d) > 1.4 and (sdo.D+d)/(sdo.D-d) < 2\
                            and (sdo.D + d) < r2 - r1\
                            and (sdo.D + d) < (r1 + r2)*math.sin(angle/2)\
                            and sdo.l1 < l1_max:
                                geometry = {'pattern' : pattern, 'radius' : (r1 + r2)/2, 'angle' : angle}
                                assembly = SpringAssembly([Spring(sdo.D, d, n, sdo.l0, material) for j in range(i)], geometry)
                                self.assemblies.append(assembly)
                                
        elif pattern == 'shaft mounted':
            for i in n_springs:
                F1eq = F1/i
                F2eq = F2/i
                for d in diameters_m:
                    for n in n_spires:
                        for material in materials:
                            sdo = SpringDiscreteOptimizer(F1eq, F2eq, stroke, d, n, material)
                            if sdo.tau_k < material.tau_max\
                            and d >= material.d_min and d <= material.d_max\
                            and sdo.D/d > 5 and sdo.D/d < 13\
                            and (sdo.D+d)/(sdo.D-d) > 1.4 and (sdo.D+d)/(sdo.D-d) < 2\
                            and (sdo.D - d) > r1 and (sdo.D + d) < r2\
                            and sdo.l1 < l1_max:
                                geometry = {'pattern' : pattern, 'radius' : None, 'angle' : None}
                                assembly = SpringAssembly([Spring(sdo.D, d, n, sdo.l0, material) for j in range(i)], geometry)
                                self.assemblies.append(assembly)
                                            
    def TargetStiffness(self):
        k = (self.F2 - self.F1)/self.stroke
         
        return k
    
class SpringAssemblyOptimizationResults(DessiaObject):
    def __init__(self, assemblies, input_data):
        self.type='mc_spring_assembly'
        self.assemblies = assemblies
        self.input_data = input_data
        self.l0_assemb = [assembly.l0 for assembly in assemblies]
        self.cost_assemb = [assembly.Cost() for assembly in assemblies]
        
        self.p_frontX, self.p_frontY, index = self.ParetoFrontier(self.cost_assemb, self.l0_assemb, False, False)
        
        self.results = [assemblies[i] for i in index]
        
        self.catalog_optimization_results = self.CatalogStudy()
        accuracy_rate = 0.05
        [custom_spring.Custom2Catalog(accuracy_rate) for custom_assembly in self.results
                                                     for custom_spring in custom_assembly.springs]
        [custom_assembly.Custom2Catalog(accuracy_rate) for custom_assembly in self.results]
        
    def CatalogStudy(self):
        spring_spec = self.input_data[0]
        catalog_spec = self.input_data[1]
        n_springs = [i + spring_spec['n_springs1'] for i in range(spring_spec['n_springs2'] - spring_spec['n_springs1'] + 1)]
        
        for cat_name, catalog in catalogs.items():
            co = CatalogOptimizer(catalog,
                                  spring_spec['F1'],
                                  spring_spec['F2'],
                                  spring_spec['stroke'],
                                  catalog_spec['stiffness_precision'],
                                  spring_spec['l1_max'],
                                  spring_spec['r1'],
                                  spring_spec['r2'],
                                  n_springs,
                                  spring_spec['pattern'].lower())
            
            dictionnary = {cat_name : co.opti_indices}
        cor = CatalogOptimizationResults(dictionnary, list(catalogs.keys()), self.input_data, catalog_spec['prod_volume'])
        
        return cor
        
    def PlotResults(self):
        plt.figure()
        plt.plot(self.cost_assemb, self.l0_assemb, 'b.', label = 'Assemblies')
        plt.plot(self.p_frontX, self.p_frontY, 'r', label = 'Pareto Frontier')
        plt.xlabel('Cost')
        plt.ylabel('Mass (kg)')
        plt.legend()
        
    def ParetoFrontier(self, Xs, Ys, maxX = True, maxY = True):
        # Sort the list in either ascending or descending order of X
        myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
    
        # Start the Pareto frontier with the first value in the sorted list
        p_front = [myList[0]]    
    
        # Loop through the sorted list
        for pair in myList[1:]:
            if maxY: 
                if pair[1] >= p_front[-1][1]: # Look for higher values of Y…
                    p_front.append(pair) # … and add them to the Pareto frontier
            else:
                if pair[1] <= p_front[-1][1]: # Look for lower values of Y…
                    p_front.append(pair) # … and add them to the Pareto frontier
    
        # Turn resulting pairs back into a list of Xs and Ys
        p_frontX = [pair[0] for pair in p_front]
        p_frontY = [pair[1] for pair in p_front]
        index = [i for i, input_pair in enumerate(myList) if [Xs[i], Ys[i]] in p_front]
        return p_frontX, p_frontY, index
    
    def Dict(self):
        catalog_spec = self.input_data[1]
        
        d={}
        assemblies=[]
        results = []
        product_assemblies = []
        for assembly in self.assemblies:
            assembly_d = assembly.Dict()
            assemblies.append(assembly_d)
            if assembly in self.results:
                results.append(assembly_d)
                
        for product_assembly in self.catalog_optimization_results.results:
            product_assemblies.append(product_assembly.Dict(catalog_spec['prod_volume']))
                    
        d['assemblies'] = assemblies
        d['input_data'] = self.input_data
        d['results'] = results
        d['catalog_results'] = product_assemblies
        
        return d
    
    
class Catalog(DessiaObject):
    def __init__(self, csv_file, name = ''):
        self.csv_file = csv_file
        self.products = pd.read_csv(csv_file)
        self.name = name
        
        self.CorrectionDynParameters()
        
    def CorrectionDynParameters(self):
        self.products['Fndyn'] = pd.to_numeric(self.products['Fndyn'], errors = 'coerce').fillna(0)
        self.products['shdyn'] = pd.to_numeric(self.products['shdyn'], errors = 'coerce').fillna(0)
        self.products['Lndyn'] = pd.to_numeric(self.products['Lndyn'], errors = 'coerce').fillna(0)
        
    def FindSprings(self, specs):
        products = copy(self.products)
        for key, bounds in specs.items():
            products = products[products[key] > bounds[0]]
            products = products[products[key] < bounds[1]]
            
        indices = products.index
        
        return indices
    
    def PlotStats(self):
        scatter_matrix(self.products, alpha=0.2, figsize=(6, 6), diagonal='kde')
        
        plt.figure()
        plt.hist(self.products['D']/self.products['d'])
        
    def Price(self, product_index, n_product):        
        current_dict = eval(self.products['prices'][product_index])
        
        keys = list(current_dict.keys())
        stop = False
        while not stop:
            max_key = max(keys)
            if n_product >= min(keys):
                if n_product >= max_key:
                    price = current_dict[max_key]
                    stop = True
                else:
                    keys.remove(max_key)
            else:
                 price = current_dict[min(keys)]   
                
        return price
    
    def Instantiate(self, product_indices, n_springs, pattern, r1, r2):
        for product_index in product_indices:
            spring = Spring(self.products['D'][product_index],
                            self.products['d'][product_index],
                            self.products['n'][product_index], # /!\ n_spires 0.5 & 0.25 : a vérifier
                            self.products['L0'][product_index])
                            # /!\ Voir pour matériaux
            
            springs = [spring]*n_springs
            if n_springs != 1:
                if pattern == 'circular':
                    angle = 2*math.pi/n_springs
                    geometry = {'pattern' : pattern, 'radius' : (r1 + r2)/2, 'angle' : angle}
                    springs = SpringAssembly(springs, geometry)
                
        return springs


ferroflex_file = pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                                               'mechanical_components/catalogs/ferroflex.csv')
ferroflex_catalog = Catalog(ferroflex_file, 'Ferroflex')

catalogs = {ferroflex_catalog.name : ferroflex_catalog}

class Product(DessiaObject):
    def __init__(self, catalog_name, product_index):
        self.catalog_name = catalog_name
        if isinstance(product_index, npy.generic):
            self.product_index = npy.asscalar(product_index)
        else:
            self.product_index = product_index
    
#    def __getstate__(self):
#        d=self.__dict__.copy()
#        print(type(d['poduct_index']))
#        if isinstance(d['poduct_index'], npy.generic):
#            print('ok')
#            d['poduct_index'] = npy.asscalar(d['poduct_index'])
#        return d
    
    def Instantiate(self):
        catalog=catalogs[self.catalog_name]
        spring = Spring(catalog.products['D'][self.product_index],
                        catalog.products['d'][self.product_index],
                        catalog.products['n'][self.product_index], # /!\ n_spires 0.5 & 0.25 : a vérifier
                        catalog.products['L0'][self.product_index])
                        # /!\ Voir pour matériaux
                        
        return spring
    
    def Dict(self):
        d = self.__dict__
        
        spring = self.Instantiate()
        d['k'] = spring.k
        d['l0'] = spring.l0
        d['D'] = spring.D
        d['d'] = spring.d
        
        return d

class ProductAssembly(DessiaObject):
    def __init__(self, products, geometry):
        self.products = products
        self.n_products = len(products)
        self.geometry = geometry
        
        self.l0 = self.FreeLength()
                
    def Instantiate(self, geometry):
        springs = [product.Instantiate() for product in self.products]
        
        spring_assembly = SpringAssembly(springs, geometry)
        
        return spring_assembly
    
    def FreeLength(self):
        product = self.products[0]
        catalog = catalogs[product.catalog_name]
        l0 = catalog.products['L0'][product.product_index]
        
        return l0
    
    def Price(self, prod_volume):
        for product in self.products:
            catalog = catalogs[product.catalog_name]
            assembly_price = sum([catalog.Price(product.product_index, self.n_products*prod_volume) for product in self.products])
        
        return assembly_price
    
    def Dict(self, prod_volume = None):
        d = self.__dict__.copy()
        products = []
        for product in self.products:
            products.append(product.Dict())
        d['products'] = products
        
        if prod_volume is not None:
            d['price'] = self.Price(prod_volume)
        
        spring_assembly = self.Instantiate(self.geometry)
        d['k'] = spring_assembly.k
        
        return d
    
        
class CatalogOptimizer(DessiaObject):
    def __init__(self, catalog, F1, F2, stroke, target_stiffness_percentage,
                 max_l1, r1 = 0.090, r2 = 0.120, n_springs = [1], pattern = 'shaft mounted',
                 prod_volume = 50):
        self.catalog = catalog
        self.F1 = F1
        self.F2 = F2
        self.stroke = stroke
        self.target_k_percentage = target_stiffness_percentage
        self.max_l1 = max_l1
        self.r1 = r1
        self.r2 = r2
        self.n_springs = n_springs
        self.pattern = pattern
        self.prod_volume = prod_volume
        
        self.opti_indices = self.Optimize()
        
    def TargetStiffness(self, F1eq, F2eq):
        target_k = (F2eq - F1eq)/self.stroke
        
        range_k = ((1-self.target_k_percentage)*target_k, (1+self.target_k_percentage)*target_k)
        
        return range_k
    
    def StiffnessConstraint(self, range_k):
        indices = self.catalog.FindSprings({'R' : range_k})
        
        return indices
    
    def InitialLengthConstraint(self, F1eq):
        products = copy(self.catalog.products)
        products = products[products['L0']*(F1eq/products['R']) < self.max_l1]
        
        indices = products.index
        
        return indices
        
    def MaxForceConstraint(self, F2eq):
        indices = self.catalog.FindSprings({'Fndyn' : (F2eq, math.inf)})
        
        return indices
        
    def SizeConstraint(self, angle):
        products = copy(self.catalog.products)
        if self.pattern == 'shaft mounted':
            products = products[products['D'] - products['d'] > self.r1]
            products = products[products['D'] + products['d'] < self.r2]
        elif self.pattern == 'circular':
            products = products[products['D'] + products['d'] < (self.r2 - self.r1)]
            products = products[products['D'] + products['d'] < (self.r2 + self.r1)*math.sin(angle/2)]
            
        indices = products.index
        
        return indices
            
    def Optimize(self):
        indices_dict = {}
        for ns in self.n_springs:
            F1eq = self.F1/ns
            F2eq = self.F2/ns
            angle = 2*math.pi/ns
            
            current_range_k = self.TargetStiffness(F1eq, F2eq)
            indices_StC = self.StiffnessConstraint(current_range_k)
            indices_ILC = self.InitialLengthConstraint(F1eq)
            indices_MFC = self.MaxForceConstraint(F2eq)
            indices_SiC = self.SizeConstraint(angle)
            
            indices = [indice for indice in self.catalog.products.index
                       if indice in indices_StC
                       and indice in indices_ILC
                       and indice in indices_MFC
                       and indice in indices_SiC]
            
            indices_dict[ns] = indices
            
        return indices_dict
    
class CatalogOptimizationResults(DessiaObject):
    def __init__(self, indices_dicts, catalogs_names, input_data, prod_volume = 50):
        self.indices_dicts = indices_dicts
        self.catalogs_names = catalogs_names
        spring_spec = input_data[0]
        
        self.X = []
        self.Y = []
        self.product_assemblies = []
        
        for cat_name, indices_dict in indices_dicts.items():
             for ns, indices in indices_dict.items():
                 for product_index in indices:
                     geometry = {'pattern' : spring_spec['pattern'],
                                 'radius' : (spring_spec['r1'] + spring_spec['r2'])/2,
                                 'angle' : 2*math.pi/ns}
                     
                     product = Product(cat_name, product_index)
                     product_assembly = ProductAssembly([product]*ns, geometry)
                     
                     price = product_assembly.Price(prod_volume)
                     l0 = product_assembly.l0
                     
                     self.X.append(l0)
                     self.Y.append(price)
                     self.product_assemblies.append(product_assembly)
        
        self.p_frontX, self.p_frontY, index = self.ParetoFrontier(self.X, self.Y, False, False)
        
        self.results = [self.product_assemblies[i] for i in index]
        
#    def __getstate__(self):
#        d=self.__dict__.copy()
#        del d['catalogs']
#        return d
        
    def PlotResults(self):
        plt.figure()
        plt.plot(self.X, self.Y, 'b.', label = 'Assemblies')
        plt.plot(self.p_frontX, self.p_frontY, 'r', label = 'Pareto Frontier')
        plt.xlabel('Length')
        plt.ylabel('Price')
        plt.legend()
        
    def ParetoFrontier(self, Xs, Ys, maxX = True, maxY = True):
        # Sort the list in either ascending or descending order of X
        myList = sorted([[Xs[i], Ys[i]] for i in range(len(Xs))], reverse=maxX)
    
        # Start the Pareto frontier with the first value in the sorted list
        p_front = [myList[0]]    
    
        # Loop through the sorted list
        for pair in myList[1:]:
            if maxY: 
                if pair[1] >= p_front[-1][1]: # Look for higher values of Y…
                    p_front.append(pair) # … and add them to the Pareto frontier
            else:
                if pair[1] <= p_front[-1][1]: # Look for lower values of Y…
                    p_front.append(pair) # … and add them to the Pareto frontier
    
        # Turn resulting pairs back into a list of Xs and Ys
        p_frontX = [pair[0] for pair in p_front]
        p_frontY = [pair[1] for pair in p_front]
        index = [i for i, input_pair in enumerate(myList) if [Xs[i], Ys[i]] in p_front]
        return p_frontX, p_frontY, index        
        
        
    
