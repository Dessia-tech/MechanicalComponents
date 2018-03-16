# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 11:20:18 2018

@author: jezequel
"""

import math
import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

from scipy.optimize import minimize

class Clutch:
    """
    Defines a wet clutch object
    """
    def __init__(self,hydraulic_cylinder,
                 plate_inner_radius = 0.090, plate_outer_radius = 0.120,
                 separator_plate_width = 0.0018, friction_plate_width = 0.0008, friction_paper_width = 0.0004,
                 separator_tooth_type = 'outer', clearance = 0.0002, n_friction_plates = 4,
                 oil_dynamic_viscosity = 0.062, oil_volumic_mass = 875,
                 input_flow = 1.5*(1/6)*10**-4,
                 max_pressure = 5000000, max_time = 0.2, max_drag_torque = 50,
                 name = ''):
        # Plates parameters
        self.plate_inner_radius = plate_inner_radius
        self.plate_outer_radius = plate_outer_radius
        self.separator_plate_width = separator_plate_width
        self.friction_plate_width = friction_plate_width
        self.friction_paper_width = friction_paper_width
        self.separator_tooth_type = separator_tooth_type
        self.clearance = clearance
        self.n_friction_plates = n_friction_plates
        self.n_separator_plates = n_friction_plates + 1
        self.max_time = max_time
        
        # Oil parameters
        self.oil_dynamic_viscosity = oil_dynamic_viscosity
        self.oil_volumic_mass = oil_volumic_mass
        self.input_flow = input_flow
        self.max_pressure = max_pressure
        
        # Hydraulic cylinder part
        self.hydraulic_cylinder = hydraulic_cylinder
        self.name = name
        
        # Contours and volumes
        self.separator_plate_contours = self.SeparatorPlateContour()
        self.separator_plate_volume = self.SeparatorPlateVolume()
        self.friction_plate_contours = self.FrictionPlateContour()
        self.friction_plate_volume = self.FrictionPlateVolume()
        
    def Update(self, values):
        for key,value in values.items():
            self.hydraulic_cylinder
            setattr(self,key,value)
        
        # Contours and volumes
        if math.isnan(self.plate_outer_radius):
            print(values)
        self.separator_plate_contours = self.SeparatorPlateContour()
        self.separator_plate_volume = self.SeparatorPlateVolume()
        self.friction_plate_contours = self.FrictionPlateContour()
        self.friction_plate_volume = self.FrictionPlateVolume()            
    
        
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
        
    
    def ClosedTransferredTorque(self):
        """
        Calculs the transferred torque when clutch is engaged
        """
        n = self.n_friction_plates
        nu = 1
        r2 = self.plate_outer_radius
        r1 = self.plate_inner_radius
        F = 100
        
        tranferred_torque = n*(2/3)*nu*F*(r2**3-r1**3)/(r2**2-r1**2)
        
        return tranferred_torque
    
    def PlatePressure(self):
        """
        Calculs the pressure applied on the clutch plates
        """
        plate_contact_area = 2*math.pi*(self.plate_outer_radius - self.plate_inner_radius)
        
        F = self.hydraulic_cylinder.PistonForce() - self.hydraulic_cylinder.SpringResultingForce()
        
        pressure = F/plate_contact_area
        
        return pressure
    
    def EngagingTime(self):
        I = 0.1 # A définir
        Ct = self.ClosedTransferredTorque()
        delta_omega_0 = 100 # A définir
        
        # Cas Ct indépendant du temps
        tf = I*delta_omega_0/Ct
        
        # Cas Ct linéaire par rapport au temps
#        tf = 2*I*delta_omega_0/Ct
        
        return tf
        
    def Mass(self):
        """ 
        Calculs the mass of the entire clutch
        """
        steel_volumic_mass = 7500
        paper_volumic_mass = 1000

        # Separator plates
        sep_outer_contour = self.separator_plate_contours[0]
        sep_inner_contour = self.separator_plate_contours[1]
        sep_volume = (self.n_separator_plates + 1)*(sep_outer_contour.Area() - sep_inner_contour.Area())*self.separator_plate_width # Last plate is 2 times wider
        
        sep_mass = steel_volumic_mass*sep_volume
        
        # Friction plates
        fric_outer_contour = self.friction_plate_contours[0]
        fric_inner_contour = self.friction_plate_contours[1] 
        fric_volume = self.n_friction_plates*(fric_outer_contour.Area() - fric_inner_contour.Area())*self.friction_plate_width
        
        fric_paper_volume = math.pi*self.friction_paper_width*(self.plate_outer_radius**2-self.plate_inner_radius**2)
        
        fric_mass = fric_volume*steel_volumic_mass + 2*fric_paper_volume*paper_volumic_mass
        
        mass = sep_mass + fric_mass + self.hydraulic_cylinder.Mass()
        
        return mass
        
    def SeparatorPlateContour(self):
        """
        Defines the separator plates contours (tooth profil & shaft interface)
        """
        # Geometry. /!\ Define variables in __init__ and not here
        n_dent = 20
        Ddent = 0.010
#        alpha = math.pi/3
        theta = math.pi/n_dent
        beta12 = theta/2
        
        r1 = self.plate_inner_radius
        r2 = self.plate_outer_radius
            
        p0 = vm.Point2D((0, 0))
        
        if self.separator_tooth_type == 'outer':
            h = 0.005
            
            # Points definition
            p1 = p0.Translation((Ddent/2, r2 + h))
            p2 = vm.Point2D((r2*math.sin(beta12), r2*math.cos(beta12)))
            p3 = p2.Rotation(p0, -beta12/2)
            p4 = p3.Rotation(p0, -beta12/2)
            p5 = p0.Translation((-Ddent/2, r2 + h))
            p6 = vm.Point2D((-r2*math.sin(beta12), r2*math.cos(beta12)))
            p7 = p6.Rotation(p0, beta12/2)
            p8 = p7.Rotation(p0, beta12/2)
            
            # Plate inner interface
            circle = vm.Circle2D(p0, self.plate_inner_radius)
        
            # Lines definition
            try:
                l1 = vm.Arc2D(p4, p3, p2)
                l2 = vm.Line2D(p2, p1)
                l3 = vm.Line2D(p1, p5)
                l4 = vm.Line2D(p5, p6)
                l5 = vm.Arc2D(p6, p7, p8)
            except:
                print(p0.vector)
                print('r1 = ', r1, '; r2 = ', r2, '; beta12 = ', beta12)
                print(p4.vector, p3.vector, p2.vector)
                
            
            primitives = [l1, l2, l3, l4, l5]
            
            # Rotate tooth pattern
            new_primitives = primitives[:]
            for i in range(n_dent-1):
                new_primitives.extend([j.Rotation(p0, (i+1)*2*theta, True) for j in primitives])
                
            # Contours definition. /!\ Outer contour needs to be appended first 
            plate_contours = [vm.Contour2D(new_primitives), vm.Contour2D([circle])]
            
        elif self.separator_tooth_type == 'inner':
            h = 0.005
            
            # Points definition
            p1 = p0.Translation(((Ddent/2), r1 - h))
            p2 = vm.Point2D((r1*math.sin(beta12), r1*math.cos(beta12)))
            p3 = p2.Rotation(p0, -beta12/2)
            p4 = p3.Rotation(p0, -beta12/2)
            p5 = p0.Translation((-Ddent/2, r1 - h))
            p6 = vm.Point2D((-r1*math.sin(beta12), r1*math.cos(beta12)))
            p7 = p6.Rotation(p0, beta12/2)
            p8 = p7.Rotation(p0, beta12/2)
            
            # Plate inner interface
            circle = vm.Circle2D(p0, self.plate_outer_radius)
        
            # Lines definition
            l1 = vm.Arc2D(p4, p3, p2)
            l2 = vm.Line2D(p2, p1)
            l3 = vm.Line2D(p1, p5)
            l4 = vm.Line2D(p5, p6)
            l5 = vm.Arc2D(p6, p7, p8)
            
            primitives = [l1, l2, l3, l4, l5]
            
            # Rotate tooth pattern
            new_primitives = primitives[:]
            for i in range(n_dent-1):
                new_primitives.extend([j.Rotation(p0, (i+1)*2*theta, True) for j in primitives])
                
            # Contours definition. /!\ Outer contour needs to be appended first 
            plate_contours = [vm.Contour2D([circle]), vm.Contour2D(new_primitives)]
            
        return plate_contours
        
    def FrictionPlateContour(self):
        """
        
        """
        # Geometry. /!\ Define variables in __init__ and not here
        n_dent = 20
        Ddent = 0.010
#        alpha = math.pi/3
        theta = math.pi/n_dent
        beta12 = theta/2
        
        r1 = self.plate_inner_radius
        r2 = self.plate_outer_radius
            
        p0 = vm.Point2D((0, 0))
        
        if self.separator_tooth_type == 'outer':
            # => Friction plates type : inner
            h = 0.005
            
            # Points definition
            p1 = p0.Translation(((Ddent/2), r1 - h))
            p2 = vm.Point2D((r1*math.sin(beta12), r1*math.cos(beta12)))
            p3 = p2.Rotation(p0, -beta12/2)
            p4 = p3.Rotation(p0, -beta12/2)
            p5 = p0.Translation((-Ddent/2, r1 - h))
            p6 = vm.Point2D((-r1*math.sin(beta12), r1*math.cos(beta12)))
            p7 = p6.Rotation(p0, beta12/2)
            p8 = p7.Rotation(p0, beta12/2)
            
            # Plate inner interface
            circle = vm.Circle2D(p0, self.plate_outer_radius)
        
            # Lines definition
            l1 = vm.Arc2D(p4, p3, p2)
            l2 = vm.Line2D(p2, p1)
            l3 = vm.Line2D(p1, p5)
            l4 = vm.Line2D(p5, p6)
            l5 = vm.Arc2D(p6, p7, p8)
            
            primitives = [l1, l2, l3, l4, l5]
            
            # Rotate tooth pattern
            new_primitives = primitives[:]
            for i in range(n_dent-1):
                new_primitives.extend([j.Rotation(p0, (i+1)*2*theta, True) for j in primitives])
                
            # Contours definition. /!\ Outer contour needs to be appended first 
            plate_contours = [vm.Contour2D([circle]), vm.Contour2D(new_primitives)]
            
        elif self.separator_tooth_type == 'outer':
            # Friction plates type : outer
            h = 0.015
            
            # Points definition
            p1 = p0.Translation((Ddent/2, r2 + h))
            p2 = vm.Point2D((r2*math.sin(beta12), r2*math.cos(beta12)))
            p3 = p2.Rotation(p0, -beta12/2)
            p4 = p3.Rotation(p0, -beta12/2)
            p5 = p0.Translation((-Ddent/2, r2 + h))
            p6 = vm.Point2D((-r2*math.sin(beta12), r2*math.cos(beta12)))
            p7 = p6.Rotation(p0, beta12/2)
            p8 = p7.Rotation(p0, beta12/2)
            
            # Plate inner interface
            circle = vm.Circle2D(p0, self.plate_inner_radius)
        
            # Lines definition
            l1 = vm.Arc2D(p4, p3, p2)
            l2 = vm.Line2D(p2, p1)
            l3 = vm.Line2D(p1, p5)
            l4 = vm.Line2D(p5, p6)
            l5 = vm.Arc2D(p6, p7, p8)
            
            primitives = [l1, l2, l3, l4, l5]
            
            # Rotate tooth pattern
            new_primitives = primitives[:]
            for i in range(n_dent-1):
                new_primitives.extend([j.Rotation(p0, (i+1)*2*theta, True) for j in primitives])
                
            # Contours definition. /!\ Outer contour needs to be appended first 
            plate_contours = [vm.Contour2D(new_primitives), vm.Contour2D([circle])]
            
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
            
            plate_volume = primitives3D.ExtrudedProfile(p0, xp, yp, self.separator_plate_contours, (0, 0, width))
            
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
            
            plate_volume = primitives3D.ExtrudedProfile(p0, xp, yp, self.friction_plate_contours, (0, 0, self.friction_plate_width))
            friction_paper_1_volume = primitives3D.HollowCylinder(pp1_coord, zp_coord, self.plate_inner_radius, self.plate_outer_radius, self.friction_paper_width)
            friction_paper_2_volume = primitives3D.HollowCylinder(pp2_coord, zp_coord, self.plate_inner_radius, self.plate_outer_radius, self.friction_paper_width)
        
            primitives.extend([plate_volume, friction_paper_1_volume, friction_paper_2_volume])
        
        return primitives

    def CADExport(self):
        volumes =[]
        
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
    
class HydraulicCylinder:
    """
    Defines a hydraulic cylinder object
    """
    def __init__(self, inner_radius = 0.020, outer_radius = 0.050,
                 chamber_width = 0.100, thickness = 0.0005, engaged_chamber_pressure = 500000,
                 n_springs = 6,
                 spring_young_modulus = 80000, spring_poisson_ratio = 0.33,
                 spring_n_windings = 10, spring_wire_diameter = 0.001, spring_outer_diameter = 0.01, 
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
        
        self.chamber_contour = self.ChamberContour()
        self.chamber_volume = self.ChamberVolume()
        
        self.piston_rod_contour = self.PistonRodContour()
        self.piston_volume = self.PistonVolume()
        
        self.spring_contour = self.SpringContour()
        self.spring_volume = self.SpringVolume()
        
        # Force
        self.spring_stiffness = self.SpringStiffness()
        self.spring_resulting_force = self.SpringResultingForce()
        self.engaged_chamber_pressure = engaged_chamber_pressure
        
    def PistonForce(self):
        piston_area = 2*math.pi*self.outer_radius
        
        F = self.engaged_chamber_pressure*piston_area
        return F     
    
    def SpringResultingForce(self):
        l0 = self.spring_free_length
        l = self.spring_final_length
        
        F = self.n_springs*(-self.spring_stiffness*(l - l0))
        
        return F
        
    def SpringStiffness(self):
        """
        Spring stiffness calculation
        """
        E = self.spring_young_modulus
        nu = self.spring_poisson_ratio
        n = self.spring_n_windings
        d = self.spring_wire_diameter
        D = self.spring_outer_diameter
        
        G = E/(2*(1+nu))
        
        stiffness = G*d**4/(8*n*D**3)
        
        return stiffness
    
    def Mass(self):
        """
        Calculs the mass of the hydraulinc cylinder
        """
        chamber_volumic_mass = 7500
        piston_volumic_mass = 7500
        spring_volumic_mass = 7500
        
        # Chamber
        Vint = math.pi*self.chamber_width*self.inner_radius**2
        Vext = math.pi*(self.chamber_width + self.thickness)*self.outer_radius**2
        Vchamber = Vext - Vint
        
        # Piston
        Vpiston = sum([i.Volume() for i in self.piston_volume])
        
#        # Spring
#        Vspring = sum([i.Volume() for i in self.spring_volume])
        
        mass = Vchamber*chamber_volumic_mass + Vpiston*piston_volumic_mass #+ Vspring*spring_volumic_mass
        
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
        #zp = vm.Vector3D((0, 0, 1))
        
        
        primitives = []
        chamber_volume = primitives3D.RevolvedProfile(p0, zp, yp, [self.chamber_contour], p0, zp, 2*math.pi, 'cylinder_chamber')
        
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
        
        piston_rod_line = primitives2D.RoundedLines2D([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10], {2:0.002, 3:0.002, 4:0.002, 7:0.002, 8:0.002, 9:0.002}, True)
        piston_rod_contour = vm.Contour2D([piston_rod_line])
        return piston_rod_contour
    
    def PistonVolume(self):
        """
        Defines the piston volume (piston rod & piston head)
        """
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
    
    def SpringContour(self):
        p0_coord = (0, 0)
        p0 = vm.Point2D(p0_coord)
        
        pc = p0.Translation((self.spring_outer_diameter/2, 0))
        
        circle = vm.Circle2D(pc, self.spring_wire_diameter/2)
        minimize
        contour = vm.Contour2D([circle])
        
        return contour
    
    def SpringVolume(self):
        p0_coord = (0, 0, 0)
        zp_coord = (0, 0, 1)
        
        p0 = vm.Point3D(p0_coord)
        pc = p0.Translation((self.spring_outer_diameter/2, 0, 0))
        
        xp = vm.Vector3D((1, 0, 0))
        yp = vm.Vector3D((0, 1, 0))
        zp = vm.Vector3D(zp_coord)
        
        primitives = []
        volume = primitives3D.HelicalExtrudedProfile(p0, xp, zp, (0, 0, 0), (0, self.spring_free_length, 0), self.spring_free_length/self.spring_n_windings, self.spring_contour)
        primitives.append(volume)
        
        return primitives
    
class ClutchOptimizer:
    
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
        
        def PressureConstraint(xa):
            return -self.clutch.PlatePressure() + self.clutch.max_pressure
        
        def TimeConstraint(xa):
            return -self.clutch.EngagingTime() + self.clutch.max_time
        
        def DragTorqueConstraint(xa, regime, delta_p):
            return max(self.clutch.DragTorque(regime, delta_p)) - self.clutch.max_drag_torque
        
        regime = npy.linspace(0, 2500*math.pi/30, 100)
        delta_p = 0
        
        fun_constraints = [{'type' : 'ineq', 'fun' : PressureConstraint},
                           {'type' : 'ineq', 'fun' : TimeConstraint},
                           {'type' : 'ineq', 'fun' : DragTorqueConstraint, 'args' : [regime, delta_p]}]
        
        xra0 = npy.random.random(self.n)
        res = minimize(Objective, xra0, constraints = fun_constraints, bounds = [(0, 1)]*self.n)
        return res


                


