import math
import numpy as np
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 11:20:18 2018

@author: jezequel
"""

class Clutch:
    """
    Defines a wet clutch object
    """
    def __init__(self, plate_inner_radius, plate_outer_radius, clearance, n_friction_plates,
                 oil_dynamic_viscosity, oil_volumic_mass, input_flow, max_pressure, hydraulic_cylinder, name=''):
        self.plate_inner_radius = plate_inner_radius
        self.plate_outer_radius = plate_outer_radius
        self.clearance = clearance
        self.n_friction_plates = n_friction_plates
        self.n_separator_plates = n_friction_plates + 1
        
        self.oil_dynamic_viscosity = oil_dynamic_viscosity
        self.oil_volumic_mass = oil_volumic_mass
        self.input_flow = input_flow
        self.hydraulic_cylinder = hydraulic_cylinder
        self.name = name
        
        self.plate_contour = self.ToothContour()
        self.plate_volume = self.ToothVolume()
        
    def DragTorque(self, omega, delta_p):
        """
        Calculs drag torque when clutch discs are disengaged
        """
        N = self.n_friction_plates
        mu = self.oil_dynamic_viscosity
        rho = self.oil_volumic_mass
        h = self.clearance
        r1 = self.plate_inner_radius
        r2 = self.plate_outer_radius
        Q = self.input_flow
        
        drag_torque = []
               
        for i, val_omega in enumerate(omega):
            Qn = (((6*mu)/(math.pi*h**3))*math.log(r1/r2)\
                + math.sqrt((((6*mu)/(math.pi*h**3))*math.log(r1/r2))**2 - (81*rho**2*val_omega**2*(r2**-2-r1**-2)*(r1**2-r2**2)-540*rho*(r2**-2-r1**-2)*delta_p)/(700*math.pi**2*h**2)))\
                    /(((27*rho)/(70*math.pi**2*h**2))*(r2**-2-r1**-2))
            
            if Qn <= Q:
                r0 = r2
            else:
                r0 = math.sqrt(Q/Qn*r2**2 + (1 - Q/Qn)*r1**2)
                
            drag_torque.append(((N*val_omega*math.pi*mu)/(2*h))*(r0**4-r1**4))
        
        return drag_torque
    
#    def TransferredTorque(self):
#        """
#        Calculs the transferred torque during Inertia phase
#        """
#        tranferred_torque = self.friction_coeff*R*N*self.hydraulic_cylinder.cylinder_force
#        
#        return tranferred_torque
        
#    def ClutchAssembly(self):
#        chamber_volume = self.hydraulic_cylinder.chamber_volume
#        piston_rod_volume = self.hydraulic_cylinder.piston_rod_volume
#        piston_head_volume = self.hydraulic_cylinder.piston_head_volume
#        
#        p0 = vm.Point3D((0, 0, 0), 'origin')
        
    def ToothContour(self):
        
        n_dent = 20
        Ddent = 0.010
        alpha = math.pi/3
        theta = math.pi/n_dent
        beta12 = theta/2
        r2 = self.plate_outer_radius
#        r2 = 0.120
        r2d = 0.136
        h = r2d - r2
        
        x = vm.Vector2D((1, 0))
        
        p0 = vm.Point2D((0, 0))
        p1 = p0.Translation((Ddent/2, r2 + h))
        p2 = vm.Point2D((r2*math.sin(beta12), r2*math.cos(beta12)))
        p3 = p2.Rotation(p0, -beta12/2)
        p4 = p3.Rotation(p0, -beta12/2)
        p5 = p0.Translation((-Ddent/2, r2 + h))
        p6 = vm.Point2D((-r2*math.sin(beta12), r2*math.cos(beta12)))
        p7 = p6.Rotation(p0, beta12/2)
        p8 = p7.Rotation(p0, beta12/2)        
        
        l1 = vm.Arc2D(p8, p7, p6)
        l2 = vm.Line2D(p6, p5)
        l3 = vm.Line2D(p5, p1)
        l4 = vm.Line2D(p1, p2)
        l5 = vm.Arc2D(p2, p3, p4)
        
        primitives = [l1, l2, l3, l4, l5]
        new_primitives = primitives[:]
        for i in range(n_dent-1):
            new_primitives.extend([j.Rotation(p0, (i+1)*2*theta, True) for j in primitives])
#        primitives.extend(new_primitives)
        
        plate_contour = vm.Contour2D(new_primitives)
        return plate_contour
        
    def ToothVolume(self):
        p0 = vm.Point3D((0, 0, 0))
        
        xp = vm.Vector3D((1, 0, 0))
        yp = vm.Vector3D((0, 1, 0))
        zp = vm.Vector3D((0, 0, 1))
        
        plate_volume = primitives3D.ExtrudedProfile(p0, xp, yp, [self.plate_contour], (0, 0, 0.00180))
        
        primitives = [plate_volume]
        
        return primitives
        
    
class HydraulicCylinder:
    """
    Defines a hydraulic cylinder object
    """
    def __init__(self, inner_radius, outer_radius, chamber_width, thickness,
                 spring_constant, spring_displacement, chamber_pressure):
        # Geometry
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.chamber_width = chamber_width
        self.thickness = thickness
#        self.chamber_line, self.chamber_contour = self.ChamberContour()
        self.chamber_contour = self.ChamberContour()
#        self.area = self.chamber_line.Area()
        self.chamber_volume = self.ChamberVolume()
        
#        self.piston_rod_line, self.piston_rod_contour = self.PistonRodContour()
        self.piston_rod_contour = self.PistonRodContour()
        self.piston_volume = self.PistonVolume()
        
        # Force
        self.spring_constant = spring_constant
        self.spring_displacement = spring_displacement
        self.chamber_pressure = chamber_pressure
        
#        self.cylinder_force = self.Force(chamber_pressure)
#        
#    def Force(self, pressure):
#        spring_force = self.spring_constant*self.spring_displacement
#        F = pressure*self.area - spring_force
#        
#        return F
        
    def ChamberContour(self):
        pc = vm.Point2D((0,0))
        
        p1 = pc.Translation((0, self.inner_radius - self.thickness))
        p2 = p1.Translation((0, (self.outer_radius - self.inner_radius) + 2*self.thickness))
        p3 = p2.Translation((self.chamber_width + self.thickness, 0))
        p4 = p3.Translation((0, -self.thickness))
        p5 = p4.Translation((-self.chamber_width, 0))
        p6 = p5.Translation((0, -self.outer_radius + self.inner_radius))
        p7 = p6.Translation((self.chamber_width, 0))
        p8 = p7.Translation((0,-self.thickness))
        
        l1 = vm.Line2D(p1, p2)
        l2 = vm.Line2D(p2, p3)
        l3 = vm.Line2D(p3, p4)
        l4 = vm.Line2D(p4, p5)
        l5 = vm.Line2D(p5, p6)
        l6 = vm.Line2D(p6, p7)
        l7 = vm.Line2D(p7, p8)
        l8 = vm.Line2D(p8, p1)

        primitives = [l1, l2, l3, l4, l5, l6, l7, l8]
        
        chamber_contour = vm.Contour2D(primitives)
        
        #l = vm.Circle2D(pc, self.outer_radius)
        #c = vm.Contour2D([l])
        
#        return l, c
        return chamber_contour
          
    def ChamberVolume(self):
        p0 = vm.Point3D((0, 0, 0))
        pc = vm.Point3D((0, 0, 1))
        xp = vm.Vector3D((1, 0, 0))
        yp = vm.Vector3D((0, 1, 0))
        zp = vm.Vector3D((0, 0, 1))
        #zp = vm.Vector3D((0, 0, 1))
        
        
        primitives = []
        chamber_volume = primitives3D.RevolvedProfile(p0, zp, yp, [self.chamber_contour], p0, zp, 2*math.pi, 'cylinder_chamber')
#        chamber_cylinder = primitives3D.HollowCylinder(p0, zp, self.outer_radius, self.outer_radius + self.thickness, self.chamber_width/2)
#        end_cylinder = primitives3D.HollowCylinder((0, 0, -(self.chamber_width/2 + self.thickness)/2), zp,self.inner_radius, self.outer_radius + self.thickness, self.thickness)
        #end_hollow_cylinder = primitives3D.HollowCylinder((0, 0, (self.chamber_width/2 + self.thickness)/2), zp, self.inner_radius, self.outer_radius + self.thickness, self.thickness)
        
        primitives.append(chamber_volume)
        return primitives
   
    def PistonRodContour(self):
        pc = vm.Point2D((0,0))
        
#        primitives = []
        
        # Piston rod
        p1 = pc.Translation((0, self.inner_radius + (self.outer_radius - self.inner_radius)/2))
        p2 = p1.Translation((0, self.thickness))
        p3 = p2.Translation((0.010, 0))
        p4 = p3.Translation((0, 0.050))
        p5 = p4.Translation((0.010, 0))
        p6 = p5.Translation((0, 0.010))
        p7 = p6.Translation((self.thickness, 0))
        p8 = p7.Translation((0, -(0.010 + self.thickness)))
        p9 = p8.Translation((-0.010, 0))
        p10 = p9.Translation((0, -0.050))
        
        piston_rod_line = primitives2D.RoundedLines2D([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10], {2:0.002, 3:0.002, 4:0.002, 7:0.002, 8:0.002, 9:0.002}, True)
        piston_rod_contour = vm.Contour2D([piston_rod_line])
        
#        l1 = vm.Line2D(p1, p2)
#        l2 = vm.Line2D(p2, p3)
#        l3 = vm.Line2D(p3, p4)
#        l4 = vm.Line2D(p4, p5)
#        l5 = vm.Line2D(p5, p6)
#        l6 = vm.Line2D(p6, p7)
#        l7 = vm.Line2D(p7, p8)
#        l8 = vm.Line2D(p8, p9)
#        l9 = vm.Line2D(p9, p10)
#        l10 = vm.Line2D(p10, p1)
#
#        primitives = [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10]
        
#        piston_rod_contour = vm.Contour2D(primitives)
        
#        l = vm.Circle2D(pc, 0.005)
#        c = vm.Contour2D([l])
#        
#        return l, c
        return piston_rod_contour
    
    def PistonVolume(self): 
        p0_coord = (0, 0, 0)
        zp_coord = (0, 0, 1)       
        p0 = vm.Point3D(p0_coord)
        xp = vm.Vector3D((1, 0, 0))
        yp = vm.Vector3D((0, 1, 0))
        zp = vm.Vector3D(zp_coord)
        
        primitives = []
        
        piston_rod = primitives3D.RevolvedProfile(p0.Translation((0, 0, 0.070 + 0.050)), zp, yp, [self.piston_rod_contour], p0, zp, 2*math.pi, 'piston_rod')
        
        piston_head = primitives3D.HollowCylinder((0, 0, 0.070), zp_coord, self.inner_radius, self.outer_radius, 0.100, 'piston_head')
        
        primitives.extend([piston_rod, piston_head])
#        primitives.append(piston_head)
        
        return primitives