#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 17:02:33 2018

@author: jezequel
"""
import math
import numpy as npy
import volmdlr as vm
import matplotlib.pyplot as plt
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

from scipy.optimize import minimize

class Material:
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

class Spring():
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
        
        self.w = self.SpringIndex()
        
        self.material = material
        self.contour = self.Contour()

    def Update(self, values):
        for key,value in values.items():
            str_split = key.split('.')
            if len(str_split) == 2:
                setattr(getattr(self, str_split[0]), str_split[1], value)
            elif len(str_split) == 1:
                setattr(self,key,value)
                
        self.p = self.l0/self.n
       
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
    
    def Contour(self):
        p0 = vm.Point2D((self.D/2, 0))
        
        l1 = vm.Circle2D(p0, self.d/2)
        
        return vm.Contour2D([l1])
    
    def Volume(self, F):
        p0_coord = (self.pos_x, self.pos_y, 0)
        zp_coord = (0, 0, 1)
        
        p0 = vm.Point3D(p0_coord)
        pc = p0.Translation((self.d/2, 0, 0))
        
        xp = vm.Vector3D((1, 0, 0))
        yp = vm.Vector3D((0, 1, 0))
        zp = vm.Vector3D(zp_coord)
        
        l = self.Length(F)
        p = l/self.n
        
        primitives = []
        volume = primitives3D.HelicalExtrudedProfile(p0, xp, zp, (self.pos_x, self.pos_y, 0), (0, 0, l), p, self.contour, name = 'spring')
        primitives.append(volume)
        
        return primitives
    
    def CADExport(self):
        volumes = []
        volumes.extend(self.volume)
        
        model = vm.VolumeModel(volumes)
        resp = model.FreeCADExport('python','spring','/usr/lib/freecad/lib/',['stl','fcstd'])
        
        return resp
    
class SpringAssembly():
    def __init__(self, springs, geometry):
        self.springs = springs
        self.n_springs = len(springs)
        self.geometry = geometry
        
        self.l0 = self.FreeLength()
        self.PositionSprings()
        
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
    
    def CADExport(self):
        volumes = []
        for spring in self.springs:
            volumes.extend(spring.Volume(0))
            
        model = vm.VolumeModel(volumes)
        resp = model.FreeCADExport('python','spring_assembly','/usr/lib/freecad/lib/',['stl','fcstd'])
        
class SpringOptimizer():    
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
        
    
class SpringDiscreteOptimizer():
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
    
class SpringAssemblyOptimizer():
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
                                            
    def TargetStiffness(self):
        k = (self.F2 - self.F1)/self.stroke
         
        return k
    
class SpringAssemblyOptimizationResults():
    def __init__(self, assemblies, bounds):
        self.assemblies = assemblies
        self.l0_assemb = [assembly.l0 for assembly in assemblies]
        self.cost_assemb = [assembly.Cost() for assembly in assemblies]
        
        self.p_frontX, self.p_frontY, index = self.ParetoFrontier(self.cost_assemb, self.l0_assemb, False, False)
        
        self.results = [assemblies[i] for i in index]
        
    def PlotResults(self):
        fig = plt.figure()
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