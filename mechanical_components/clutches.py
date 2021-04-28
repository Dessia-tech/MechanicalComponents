# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 11:20:18 2018

@author: jezequel
"""

import math
import numpy as npy

from dessia_common.core import DessiaObject

import volmdlr as vm
import volmdlr.primitives3d as primitives3d
import volmdlr.primitives2d as primitives2d

from scipy.optimize import minimize

class Clutch(DessiaObject):
    """
    Defines a wet clutch object
    """
    def __init__(self,hydraulic_cylinder,
                 plate_inner_radius = 0.090, plate_height = 0.030,
                 separator_plate_width = 0.0018, friction_plate_width = 0.0008, friction_paper_width = 0.0004,
                 separator_tooth_type = 'outer', clearance = 0.0002, n_friction_plates = 4,
                 oil_dynamic_viscosity = 0.062, oil_volumic_mass = 875,
                 input_flow = 1.5*(1/6)*10**-4,
                 max_pressure = 5000000, max_time = 0.2, max_drag_torque = 50,
                 min_torque = 30,
                 name = ''):
        
        # Plates parameters
        self.plate_inner_radius = plate_inner_radius
        self.plate_height = plate_height
        self.tooth_height = 0.005
        self.separator_plate_width = separator_plate_width
        self.friction_plate_width = friction_plate_width
        self.friction_paper_width = friction_paper_width
        self.separator_tooth_type = separator_tooth_type
        self.clearance = clearance
        self.n_friction_plates = n_friction_plates
        self.n_separator_plates = n_friction_plates + 1
        self.max_time = max_time
        self.total_clearance = self.Displacement
        
        # Oil parameters
        self.oil_dynamic_viscosity = oil_dynamic_viscosity
        self.oil_volumic_mass = oil_volumic_mass
        self.input_flow = input_flow
        self.max_pressure = max_pressure
        
        # Material parameters
        self.steel_volumic_mass = 7500
        self.paper_volumic_mass = 1000
        
        # Hydraulic cylinder part
        self.hydraulic_cylinder = hydraulic_cylinder
        self.name = name
        
        # Geometry
        self.origin = vm.Point2D((0,0))
        
        # Contours and volumes
        self.separator_plate_contours = self.SeparatorPlateContour()
        self.separator_plate_volume = self.SeparatorPlateVolume()
        self.friction_plate_contours = self.FrictionPlateContour()
        self.friction_plate_volume = self.FrictionPlateVolume()
        
        # Moments of inertia & Masses
        self.PlatesGeometryUpdate()
        
    def _get_plate_height(self):
        return self.plate_outer_radius - self.plate_inner_radius

    def _set_plate_height(self,value):
        self.plate_outer_radius = self.plate_inner_radius + value

    plate_height=property(_get_plate_height,_set_plate_height)
        
    def Update(self, values):
        for key,value in values.items():
            str_split = key.split('.')
            if len(str_split) == 2:
                setattr(getattr(self, str_split[0]), str_split[1], value)
            elif len(str_split) == 1:
                setattr(self,key,value)
        
        # Contours and volumes
        self.separator_plate_contours = self.SeparatorPlateContour()
        self.separator_plate_volume = self.SeparatorPlateVolume()
        self.friction_plate_contours = self.FrictionPlateContour()
        self.friction_plate_volume = self.FrictionPlateVolume()
        
        
        self.total_clearance = Displacement()
        self.hydraulic_cylinder.Update(self.total_clearance)
        
        # Plates geometry
        self.PlatesGeometryUpdate()
        
    def Displacement(self):
        d = 2*self.n_friction_plates*self.clearance
        
        return d
    
    def PlatesGeometryUpdate(self):
        p0 = self.origin
        
        # Temporaire (attente correction moment quadratique)
        r1 = self.plate_inner_radius
        r1eq = r1 - self.tooth_height/2
        r2 = self.plate_outer_radius
        r2eq = r2 + self.tooth_height/2
        
        self.J_paper = (self.paper_volumic_mass*math.pi*self.friction_paper_width*(r2**4-r1**4))/2
        if self.separator_tooth_type == 'outer':            
            self.J_sep = (self.n_separator_plates + 1)*(self.steel_volumic_mass*math.pi*self.separator_plate_width*(r2eq**4-r1**4))/2
            self.J_fric = self.n_friction_plates*(self.steel_volumic_mass*math.pi*self.separator_plate_width*(r2**4-r1eq**4) + 2*self.J_paper)/2
        else:
            self.J_sep = (self.n_separator_plates + 1)*(self.steel_volumic_mass*math.pi*self.separator_plate_width*(r2**4-r1eq**4))/2
            self.J_fric = self.n_friction_plates*(self.steel_volumic_mass*math.pi*self.separator_plate_width*(r2eq**4-r1**4) + 2*self.J_paper)/2
        
        # Plates are considered as cylinders and pitch radius the mean radius of outer and outer + height radii
        # Separator plate
#        sep_outer_contour = self.separator_plate_contours[0]
#        sep_inner_contour = self.separator_plate_contours[1]
#        sep_volume = (self.n_separator_plates + 1)*(sep_outer_contour.Area() - sep_inner_contour.Area())*self.separator_plate_width # Last plate is 2 times wider
        sep_volume = sum([i.Volume() for i in self.separator_plate_volume])
        
        self.sep_mass = self.steel_volumic_mass*sep_volume
#        self.sep_SMA = sep_outer_contour.SecondMomentArea(p0) - sep_inner_contour.SecondMomentArea(p0)
#        self.J_sep = (self.n_separator_plates + 1)*self.steel_volumic_mass*self.separator_plate_width*(self.sep_SMA[0, 0] + self.sep_SMA[1, 1])
        
        # Friction plate
#        fric_outer_contour = self.friction_plate_contours[0]
#        fric_inner_contour = self.friction_plate_contours[1] 
#        fric_volume = self.n_friction_plates*(fric_outer_contour.Area() - fric_inner_contour.Area())*self.friction_plate_width
#        fric_paper_volume = math.pi*self.friction_paper_width*(self.plate_outer_radius**2-self.plate_inner_radius**2)
        fric_volume = sum([i.Volume() for i in self.friction_plate_volume[:self.n_friction_plates]])
        fric_paper_volume = sum([i.Volume() for i in self.friction_plate_volume[self.n_friction_plates:]])

        
#        J_paper = (self.paper_volumic_mass*math.pi*self.friction_paper_width*(self.plate_outer_radius**4-self.plate_inner_radius**4))/2
        
        self.fric_mass = fric_volume*self.steel_volumic_mass + 2*fric_paper_volume*self.paper_volumic_mass
#        self.fric_SMA = fric_outer_contour.SecondMomentArea(p0) - fric_inner_contour.SecondMomentArea(p0)
#        self.J_fric = (self.n_friction_plates)*(self.steel_volumic_mass*self.separator_plate_width*(self.sep_SMA[0, 0] + self.sep_SMA[1, 1]) + 2*J_paper)
        
    def DragTorque(self, omega, delta_p):
        """
        Calculs drag torque when clutch discs are disengaged
        """
        # Geometry and variables
        N = self.n_friction_plates
        mu = self.oil_dynamic_viscosity
        rho = self.oil_volumic_mass
        h = self.clearance
        r1 = self.plate_inner_radius
        r2 = self.plate_outer_radius
        Q = self.input_flow
        
        drag_torque = []
               
        for i, val_omega in enumerate(omega):
            # Needed flow rate
            Qn = (((6*mu)/(math.pi*h**3))*math.log(r1/r2)\
                + math.sqrt((((6*mu)/(math.pi*h**3))*math.log(r1/r2))**2 - (81*rho**2*val_omega**2*(r2**-2-r1**-2)*(r1**2-r2**2)-540*rho*(r2**-2-r1**-2)*delta_p)/(700*math.pi**2*h**2)))\
                    /(((27*rho)/(70*math.pi**2*h**2))*(r2**-2-r1**-2))
            
            # If the needed flow rate (to have a full oil film) doesn't exceed the available flow rate
            if Qn <= Q:
                # The oil film is full and the equivalent radius is r2
                r0 = r2
            else: # I fthe needed flow rate exceeds the available flow rate
                # Oil film isn't full anymore and a equivalent radius is
                # calculated (depends on Qn & omega)
                r0 = math.sqrt(Q/Qn*r2**2 + (1 - Q/Qn)*r1**2)
                
            drag_torque.append(((N*val_omega*math.pi*mu)/(2*h))*(r0**4-r1**4))
        
        return drag_torque
    
    def SlippingTransferredTorque(self):
        """
        Calculs transferred torque when clutch is in slipping phase
        """
        return transferred_torque
        
    
    def EngagedTransferredTorque(self):
        """
        Calculs the transferred torque when clutch is engaged
        """
        n = self.n_friction_plates
        nu = 0.11 # /!\ TODO
        r2 = self.plate_outer_radius
        r1 = self.plate_inner_radius
        F = self.hydraulic_cylinder.PistonForce() - self.hydraulic_cylinder.SpringResultingForce()
        
        tranferred_torque = n*(2/3)*nu*F*(r2**3-r1**3)/(r2**2-r1**2)
        
        return tranferred_torque
    
    def PlatePressure(self):
        """
        Calculs the pressure applied on the clutch plates
        """
        plate_contact_area = 2*self.n_friction_plates*math.pi*(self.plate_outer_radius**2 - self.plate_inner_radius**2)
        
        F = self.hydraulic_cylinder.PistonForce() - self.hydraulic_cylinder.SpringResultingForce()
        
        pressure = F/plate_contact_area
        
        return pressure
    
    def EngagementTime(self, regime):
#        I = 0.1 # TODO
        I = self.J_sep + self.J_fric + self.J_paper
        Ct = self.EngagedTransferredTorque()
        delta_omega_0 = regime # TODO
        
        # Cas Ct indépendant du temps
        tf = I*delta_omega_0/Ct
        
        # Cas Ct linéaire par rapport au temps
#        tf = 2*I*delta_omega_0/Ct
        
        return tf
        
    def Mass(self):
        """ 
        Calculs the mass of the entire clutch
        """
        # Separator plates
#        sep_outer_contour = self.separator_plate_contours[0]
#        sep_inner_contour = self.separator_plate_contours[1]
#        sep_volume = (self.n_separator_plates + 1)*(sep_outer_contour.Area() - sep_inner_contour.Area())*self.separator_plate_width # Last plate is 2 times wider
#        
#        sep_mass = steel_volumic_mass*sep_volume
        
        # Friction plates
#        fric_outer_contour = self.friction_plate_contours[0]
#        fric_inner_contour = self.friction_plate_contours[1] 
#        fric_volume = self.n_friction_plates*(fric_outer_contour.Area() - fric_inner_contour.Area())*self.friction_plate_width
#        
#        fric_paper_volume = math.pi*self.friction_paper_width*(self.plate_outer_radius**2-self.plate_inner_radius**2)
#        
#        fric_mass = fric_volume*steel_volumic_mass + 2*fric_paper_volume*paper_volumic_mass
        
        mass = self.sep_mass + self.fric_mass + self.hydraulic_cylinder.Mass()
        
        return mass
        
    def SeparatorPlateContour(self):
        """
        Defines the separator plates contours (tooth profil & shaft interface)
        """
        # Geometry. /!\ Define variables in __init__ and not here
        n_dent = 20
        alpha_tirage = 1*math.pi/180
        theta = 2*math.pi/n_dent
        beta = 0.2*theta
        alpha_tot = theta-2*beta
        h = self.tooth_height
        
        r1 = self.plate_inner_radius
        r2 = self.plate_outer_radius
            
        p0 = vm.Point2D((0, 0))
        
        if self.separator_tooth_type == 'outer':
            p = p0.Translation((0, r2 + h))
            p7 = p.Rotation(p0, -theta/2)
            p67 = p7.Rotation(p0, beta/2)
            p6 = p67.Rotation(p0, beta/2)
            [x6, y6] = p6.vector
            r = y6 - x6/math.tan(alpha_tirage)
            pc = p0.Translation((0, r))
            p4 = p0.Translation((0, r2))
            p5 = p4.Rotation(pc, -alpha_tirage/2)
            p3 = p4.Rotation(pc, alpha_tirage/2)
            p2 = p6.Rotation(p0, alpha_tot)
            p12 = p2.Rotation(p0, beta/2)
            p1 = p12.Rotation(p0, beta/2)
            
            l1 = vm.Arc2D(p7, p67, p6)
            l2 = vm.Line2D(p6, p5)
            l3 = vm.Arc2D(p5, p4, p3)
            l4 = vm.Line2D(p3, p2)
            l5 = vm.Arc2D(p2, p12, p1)
            
            # Plate inner interface
            circle = vm.Circle2D(p0, r1)
                
            primitives = [l1, l2, l3, l4, l5]
            
            # Rotate tooth pattern
            new_primitives = primitives[:]
            for i in range(n_dent-1):
                new_primitives.extend([j.Rotation(p0, (i+1)*theta, True) for j in primitives])
                
            # Contours definition. /!\ Outer contour needs to be appended first 
            plate_contours = [vm.Contour2D(new_primitives), vm.Contour2D([circle])]
            
        elif self.separator_tooth_type == 'inner':
            # Points definition
            p = p0.Translation((0, r1 - h))
            p7 = p.Rotation(p0, -theta/2)
            p67 = p7.Rotation(p0, beta/2)
            p6 = p67.Rotation(p0, beta/2)
            [x6, y6] = p6.vector
            r = y6 - x6/math.tan(alpha_tirage)
            pc = p0.Translation((0, r))
            p4 = p0.Translation((0, r1))
            p5 = p4.Rotation(pc, -alpha_tirage/2)
            p3 = p4.Rotation(pc, alpha_tirage/2)
            p2 = p6.Rotation(p0, alpha_tot)
            p12 = p2.Rotation(p0, beta/2)
            p1 = p12.Rotation(p0, beta/2)
            
            l1 = vm.Arc2D(p7, p67, p6)
            l2 = vm.Line2D(p6, p5)
            l3 = vm.Arc2D(p5, p4, p3)
            l4 = vm.Line2D(p3, p2)
            l5 = vm.Arc2D(p2, p12, p1)
            
            # Plate inner interface
            circle = vm.Circle2D(p0, r2)
                
            primitives = [l1, l2, l3, l4, l5]
            
            # Rotate tooth pattern
            new_primitives = primitives[:]
            for i in range(n_dent-1):
                new_primitives.extend([j.Rotation(p0, (i+1)*theta, True) for j in primitives])
                
            # Contours definition. /!\ Outer contour needs to be appended first 
            plate_contours = [vm.Contour2D([circle]), vm.Contour2D(new_primitives)]
            
        return plate_contours
        
    def FrictionPlateContour(self):
        """
        
        """
        # Geometry. /!\ Define variables in __init__ and not here
        n_dent = 20
        alpha_tirage = 1*math.pi/180
        theta = 2*math.pi/n_dent
        beta = 0.2*theta
        alpha_tot = theta-2*beta
        h = self.tooth_height
        
        r1 = self.plate_inner_radius
        r2 = self.plate_outer_radius
            
        p0 = vm.Point2D((0, 0))
        
        if self.separator_tooth_type == 'inner':
            # => Friction plates type : outer
            # Points definition
            p = p0.Translation((0, r2 + h))
            p7 = p.Rotation(p0, -theta/2)
            p67 = p7.Rotation(p0, beta/2)
            p6 = p67.Rotation(p0, beta/2)
            [x6, y6] = p6.vector
            r = y6 - x6/math.tan(alpha_tirage)
            pc = p0.Translation((0, r))
            p4 = p0.Translation((0, r2))
            p5 = p4.Rotation(pc, -alpha_tirage/2)
            p3 = p4.Rotation(pc, alpha_tirage/2)
            p2 = p6.Rotation(p0, alpha_tot)
            p12 = p2.Rotation(p0, beta/2)
            p1 = p12.Rotation(p0, beta/2)
            
            l1 = vm.Arc2D(p7, p67, p6)
            l2 = vm.Line2D(p6, p5)
            l3 = vm.Arc2D(p5, p4, p3)
            l4 = vm.Line2D(p3, p2)
            l5 = vm.Arc2D(p2, p12, p1)
            
            circle = vm.Circle2D(p0, r1)
            
            primitives = [l1, l2, l3, l4, l5]
            
            # Rotate tooth pattern
            new_primitives = primitives[:]
            for i in range(n_dent-1):
                new_primitives.extend([j.Rotation(p0, (i+1)*theta, True) for j in primitives])
                
            # Contours definition. /!\ Outer contour needs to be appended first 
            plate_contours = [vm.Contour2D(new_primitives), vm.Contour2D([circle])]
            
        elif self.separator_tooth_type == 'outer':
            # Friction plates type : inner
            # Points definition
            p = p0.Translation((0, r1 - h))
            p7 = p.Rotation(p0, -theta/2)
            p67 = p7.Rotation(p0, beta/2)
            p6 = p67.Rotation(p0, beta/2)
            [x6, y6] = p6.vector
            r = y6 - x6/math.tan(alpha_tirage)
            pc = p0.Translation((0, r))
            p4 = p0.Translation((0, r1))
            p5 = p4.Rotation(pc, -alpha_tirage/2)
            p3 = p4.Rotation(pc, alpha_tirage/2)
            p2 = p6.Rotation(p0, alpha_tot)
            p12 = p2.Rotation(p0, beta/2)
            p1 = p12.Rotation(p0, beta/2)
            
            l1 = vm.Arc2D(p7, p67, p6)
            l2 = vm.Line2D(p6, p5)
            l3 = vm.Arc2D(p5, p4, p3)
            l4 = vm.Line2D(p3, p2)
            l5 = vm.Arc2D(p2, p12, p1)
            
            # Plate inner interface
            circle = vm.Circle2D(p0, r2)
            
            primitives = [l1, l2, l3, l4, l5]
            
            # Rotate tooth pattern
            new_primitives = primitives[:]
            for i in range(n_dent-1):
                new_primitives.extend([j.Rotation(p0, (i+1)*theta, True) for j in primitives])
                
            # Contours definition. /!\ Outer contour needs to be appended first 
            plate_contours = [vm.Contour2D([circle]), vm.Contour2D(new_primitives)]
            
        return plate_contours
        
    def SeparatorPlateVolume(self):
        """
        Defines the separator plates volume
        """
        xp = vm.Vector3D((1, 0, 0))
        yp = vm.Vector3D((0, 1, 0))
        zp = vm.Vector3D((0, 0, 1))
        
        primitives = []
        for i in range(self.n_friction_plates + 1):
            p0 = vm.Point3D((0, 0, i*(self.separator_plate_width + 2*self.clearance + 2*self.friction_paper_width + self.friction_plate_width)))
            
            if i == self.n_friction_plates:
                width = 2*self.separator_plate_width
            else:
                width = self.separator_plate_width
            
            plate_volume = primitives3d.ExtrudedProfile(p0, xp, yp, self.separator_plate_contours, (0, 0, width), 'separator_plate_{0}'.format(i))
            
            primitives.append(plate_volume)                
        
        return primitives
    
    def FrictionPlateVolume(self):
        """
        Defines the fiction plates volume
        """
        xp_coord = (1, 0, 0)
        yp_coord = (0, 1, 0)
        zp_coord = (0, 0, 1)
            
        xp = vm.Vector3D(xp_coord)
        yp = vm.Vector3D(yp_coord)
        zp = vm.Vector3D(zp_coord)
        
        primitives = []
        for i in range(self.n_friction_plates):
            p0_coord = (0, 0, i*(self.separator_plate_width + 2*self.clearance + 2*self.friction_paper_width + self.friction_plate_width) + self.separator_plate_width + self.clearance + self.friction_paper_width)
            pp1_coord = (0, 0, p0_coord[2] + self.friction_plate_width + self.friction_paper_width/2)
            pp2_coord = (0, 0, p0_coord[2] -self.friction_paper_width/2)
            
            p0 = vm.Point3D(p0_coord)
            pp1 = p0.Translation(pp1_coord)
            pp2 = p0.Translation(pp2_coord)
            
            plate_volume = primitives3d.ExtrudedProfile(p0, xp, yp, self.friction_plate_contours, (0, 0, self.friction_plate_width), 'friction_plate_{0}'.format(i))
            friction_paper_1_volume = primitives3d.HollowCylinder(pp1_coord, zp_coord, self.plate_inner_radius, self.plate_outer_radius, self.friction_paper_width, 'friction_paper_{0}1'.format(i))
            friction_paper_2_volume = primitives3d.HollowCylinder(pp2_coord, zp_coord, self.plate_inner_radius, self.plate_outer_radius, self.friction_paper_width, 'friction_paper_{0}2'.format(i))
        
            primitives.extend([plate_volume, friction_paper_1_volume, friction_paper_2_volume])
        
        return primitives

    def CADExport(self):
        volumes = []
        
        # Clutch volumes
        volumes.extend(self.separator_plate_volume)
        volumes.extend(self.friction_plate_volume)
        
        # Cylinder volumes
        volumes.extend(self.hydraulic_cylinder.chamber_volume)
        volumes.extend(self.hydraulic_cylinder.piston_volume)
        volumes.extend(self.hydraulic_cylinder.spring_volume)
        
        model = vm.VolumeModel(volumes)
        resp = model.FreeCADExport('python','clutch','/usr/lib/freecad/lib/',['stl','fcstd'])
        
        return resp
    
class HydraulicCylinder(DessiaObject):
    """
    Defines a hydraulic cylinder object
    """
    def __init__(self, inner_radius = 0.020, outer_radius = 0.050,
                 chamber_width = 0.100, thickness = 0.0005, engaged_chamber_pressure = 50000,
                 n_springs = 6,
                 spring_young_modulus = 210*10**9, spring_poisson_ratio = 0.33,
                 spring_n_windings = 5, spring_wire_diameter = 0.001, spring_outer_diameter = 0.01, 
                 spring_free_length = 0.01, spring_final_length = 0.005):
        
        # Geometry
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.chamber_width = chamber_width
        self.thickness = thickness
        self.n_springs = n_springs
        self.spring_young_modulus = spring_young_modulus
        self.spring_poisson_ratio = spring_poisson_ratio
        self.spring_n_windings = spring_n_windings
        self.spring_wire_diameter = spring_wire_diameter
        self.spring_outer_diameter = spring_outer_diameter
        self.spring_free_length = spring_free_length
        self.spring_final_length = spring_final_length
        self.displacement = spring_final_length - spring_free_length
        
        self.chamber_contour = self.ChamberContour()
        self.chamber_volume = self.ChamberVolume()
        
        self.piston_rod_contour = self.PistonRodContour()
        self.piston_volume = self.PistonVolume()
        
        self.spring_contour = self.SpringContour()
        self.spring_volume = self.SpringVolume()
        
        # Material parameter
        self.steel_volumic_mass = 7500
        
        # Force
#        self.hydraulic_force = self.HydraulicForce()
#        self.piston_force = 
        self.spring_stiffness = self.SpringStiffness()
        self.spring_resulting_force = self.SpringResultingForce()
        self.engaged_chamber_pressure = engaged_chamber_pressure
        
    def _get_height(self):
        return self.outer_radius - self.inner_radius

    def _set_height(self,value):
        self.outer_radius = self.inner_radius + value

    height=property(_get_height,_set_height)
        
    def Update(self, displacement):
        # Contour & volumes
        self.chamber_contour = self.ChamberContour()
        self.chamber_volume = self.ChamberVolume()
        
        self.piston_rod_contour = self.PistonRodContour()
        self.piston_volume = self.PistonVolume()
        
        self.spring_contour = self.SpringContour()
        self.spring_volume = self.SpringVolume()
        
        self.spring_stiffness = self.SpringStiffness()
        
        self.displacement = displacement
        
    def PressureArea(self):
        area = math.pi*(self.outer_radius**2 - self.inner_radius**2)
        return area
        
    def HydraulicForce(self):        
        F = self.engaged_chamber_pressure*self.PressureArea()
        return F     
    
    def SpringForce(self):        
        F = self.n_springs*(-self.SpringStiffness()*self.displacement)
        return F
    
    def PistonForce(self):
        force = self.HydraulicForce() - self.SpringForce()
        if force < 0:
            force = 0
        return force
    
    def SpringLengths(self):
        p0 = 10**5 # /!\ TODO lub pressure
        p1 = self.engaged_chamber_pressure
        
        htot = self.displacement
        
        l0 = htot*p0/(p1-p0)
        l1 = l0+htot
        
        return l0, l1
        
    def SpringStiffness(self):
        """
        Spring stiffness calculation
        """
        p1 = self.engaged_chamber_pressure
        S = self.PressureArea()
        
        l0, l1 = self.SpringLengths()
        
        stiffness = p1*S/l1
#        E = self.spring_young_modulus
#        nu = self.spring_poisson_ratio
#        n = self.spring_n_windings
#        d = self.spring_wire_diameter
#        D = self.spring_outer_diameter
#        
#        G = E/(2*(1+nu))
#        
#        stiffness = G*d**4/(8*n*D**3)
        
        return stiffness
    
    def Mass(self):
        """
        Calculs the mass of the hydraulinc cylinder
        """        
        # Chamber
        Vint = math.pi*self.chamber_width*self.inner_radius**2
        Vext = math.pi*(self.chamber_width + self.thickness)*self.outer_radius**2
        Vchamber = Vext - Vint
        
        # Piston
        Vpiston = sum([i.Volume() for i in self.piston_volume])
        
#        # Spring
#        Vspring = sum([i.Volume() for i in self.spring_volume])
        
        mass = Vchamber*self.steel_volumic_mass + Vpiston*self.steel_volumic_mass #+ Vspring*self.steel_volumic_mass
        
        return mass
        
    def ChamberContour(self):
        """
        Defines the contour of the cylinder chamber
        """
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
        return chamber_contour
          
    def ChamberVolume(self):
        """
        Defines the chamber volume
        """
        p0 = vm.Point3D((0, 0, 0))
        pc = vm.Point3D((0, 0, 1))
        xp = vm.Vector3D((1, 0, 0))
        yp = vm.Vector3D((0, 1, 0))
        zp = vm.Vector3D((0, 0, 1))        
        
        primitives = []
        chamber_volume = primitives3d.RevolvedProfile(p0, zp, yp, [self.chamber_contour], p0, zp, 2*math.pi, 'cylinder_chamber')
        
        primitives.append(chamber_volume)
        return primitives
   
    def PistonRodContour(self):
        """
        Defines the piston rod contour
        """
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
        
        piston_rod_line = primitives2d.RoundedLineSegments2D([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10],
                                                      {2:0.002, 3:0.002, 4:0.002, 7:0.002, 8:0.002, 9:0.002}, True)
        piston_rod_contour = vm.Contour2D([piston_rod_line])
        return piston_rod_contour
    
    def PistonVolume(self):
        """
        Defines the piston volume (piston rod & piston head)
        """     
        
        primitives = []
        
        piston_rod = primitives3d.RevolvedProfile(vm.Point3D((0, 0, 0.070 + 0.050)),
                                                  vm.z3D, vm.y3D, [self.piston_rod_contour],
                                                  vm.o3D, vm.z3D, 2*math.pi, 'piston_rod')
        
        piston_head = primitives3d.HollowCylinder((0, 0, 0.070), vm.z3D, self.inner_radius, self.outer_radius, 0.100, 'piston_head')
        
        primitives.extend([piston_rod, piston_head])
#        primitives.append(piston_head)
        
        return primitives
    
    def SpringContour(self):
        p0_coord = (0, 0)
        p0 = vm.Point2D(p0_coord)
        
        pc = p0.Translation((self.spring_outer_diameter/2, 0))
        
        circle = vm.Circle2D(pc, self.spring_wire_diameter/2)
        contour = vm.Contour2D([circle])
        
        return contour
    
    def SpringVolume(self):        
#        pc = vm.Point3D((self.spring_outer_diameter/2, 0, 0))
                
        primitives = []
        volume = primitives3d.HelicalExtrudedProfile(vm.o3D, vm.x3D, vm.z3D,
                                                     vm.o3D,
                                                     vm.Vector3D((0, self.spring_free_length, 0)),
                                                     self.spring_free_length/self.spring_n_windings,
                                                     self.spring_contour, name = 'spring')
        primitives.append(volume)
        
        return primitives
    

class ClutchOptimizer(DessiaObject):
    def __init__(self, clutch, specs):
        self.specs = specs
        self.clutch = clutch
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
        self.clutch.Update(self.fixed_values)
    
    def Optimize(self):
        def Objective(xa):
            values = {}
            for xai, attribute, bounds in zip(xa, self.attributes, self.bounds):                    
                values[attribute] = bounds[0] + (bounds[1] - bounds[0])*xai
            self.clutch.Update(values)
#            return self.clutch.DragTorque([100], 0.003)[0]
            return self.clutch.Mass()
        
#        def PressureConstraint(xa):
#            return self.clutch.max_pressure - self.clutch.PlatePressure()
        
        def TimeConstraint(xa, regime):
            return self.clutch.max_time - max(self.clutch.EngagementTime(regime))
        
#        def DragTorqueConstraint(xa, regime, delta_p):
#            return self.clutch.max_drag_torque - max(self.clutch.DragTorque(regime, delta_p))
            
        def SpringConstraint(xa):
            return self.clutch.hydraulic_cylinder.PistonForce() - self.clutch.hydraulic_cylinder.SpringResultingForce()
        
        def TorqueConstraint(xa):
            return self.clutch.EngagedTransferredTorque() - self.clutch.min_torque
        
        regime = npy.linspace(0, 2500*math.pi/30, 100)
        delta_p = 0
        
        fun_constraints = [{'type' : 'ineq', 'fun' : TimeConstraint, 'args' : [regime]},
                           {'type' : 'ineq', 'fun' : TorqueConstraint},
                           {'type' : 'ineq', 'fun' : SpringConstraint}]
#                            {'type' : 'ineq', 'fun' : PressureConstraint},
#                           {'type' : 'ineq', 'fun' : DragTorqueConstraint, 'args' : [regime, delta_p]},
        
        xra0 = npy.random.random(self.n)
        res = minimize(Objective, xra0, constraints = fun_constraints, bounds = [(0, 1)]*self.n)
        return res