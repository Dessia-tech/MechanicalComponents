#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 10:35:36 2019

@author: ringhausen
"""

import math
import numpy as npy
from scipy import interpolate

from typing import List
import json
import pkg_resources

from dessia_common.core import DessiaObject

import volmdlr as vm
import volmdlr.primitives2d
import volmdlr.primitives3d as vm3d


phi_fct_d1 = {'data':[[ 4.915, 0.09403],
                      [ 6.412, 0.10364],
                      [  7.91, 0.11283], 
                      [ 9.749, 0.12036],
                      [11.249, 0.12871],
                      [12.923, 0.13499],
                      [14.763, 0.1421 ], 
                      [17.108, 0.14964],
                      [ 19.79, 0.15719],
                      [22.638, 0.16598],
                      [ 25.49, 0.17311],
                      [28.342, 0.18024],
                      [31.197, 0.18612],
                      [33.883, 0.19158],
                      [36.572, 0.19621],
                      [38.923, 0.20083],
                      [41.442, 0.20587],
                      [  44.3, 0.21008],
                      [47.326, 0.21472],
                      [50.351, 0.21977],
                      [53.545, 0.2244                   ], 
                      [57.076, 0.22946],
                      [60.947, 0.23327],
                      [64.967, 0.23833],
                      [68.349, 0.24214],
                      [72.049, 0.24679],
                      [78.949, 0.25357],
                      [81.981, 0.2557 ],
                      [84.845, 0.25722],
                      [   320, 0.25722],
                          ], 'x':'Linear','y':'Linear'}

q_fct_wh = {'data':[[0.89106, 5.4519],
                    [0.90362, 5.2546],
                    [0.92116, 5.0719],
                    [ 0.9412, 4.8892],
                    [0.96623, 4.6991],
                    [0.99126, 4.5238],
                    [1.02127, 4.3484],
                    [1.05378, 4.1803],
                    [ 1.0838, 4.0122],
                    [1.12129, 3.8587],
                    [1.16378, 3.6979],
                    [1.19876, 3.5737],
                    [1.24374, 3.4202],
                    [1.28871, 3.2887],
                    [1.33868, 3.1498],
                    [1.38614, 3.0256],
                    [1.43859, 2.9013], 
                    [1.48855, 2.7771],
                    [  1.541, 2.6675],
                    [1.59843, 2.5652],
                    [1.65588, 2.4482],
                    [ 1.7183, 2.3386],
                    [1.77823, 2.2363],
                    [1.84315, 2.134 ],
                    [1.90557, 2.0463],
                    [  1.963, 1.9586],
                    [2.02043, 1.8782],
                    [2.08533, 1.8197],
                    [2.14775, 1.732 ],
                    [2.21016, 1.6809],
                    [2.26758, 1.6151],
                    [2.32748, 1.5786],
                    [2.38489, 1.5201],
                    [2.45829, 1.4543],
                    [2.52468, 1.4032],
                    [2.59956, 1.3593],
                    [2.67195, 1.3082],
                    [2.73934, 1.2643],
                    [ 2.8192, 1.2278],
                    [2.89158, 1.1912],
                    [2.95897, 1.162 ],
                    [3.02885, 1.1255],
                    [3.09873, 1.1035],
                    [ 3.1711, 1.0889],
                    [3.23848, 1.0743],
                    [3.31834, 1.0524],
                    [ 3.3957, 1.0378],
                    [3.46807, 1.0158],
                    [3.53296, 1.0012],
                    [3.60283, 1     ],
                    [      4, 1     ]
                          ], 'x':'Linear','y':'Linear'}


def AngleTrigo(u, v):
    dot = u.Dot(v)
    norm_u = u.Norm()
    norm_v = v.Norm()
    cross = u.Cross(v)
    inner_angle = math.acos(dot/(norm_u*norm_v))
    if cross < 0:
        return 2*math.pi-inner_angle
    else:
        return inner_angle



class CirclipsMaterial(DessiaObject):
    def __init__(self, E, St, t_shear, C, T, ID='', name=''):
        self.name = name
        self.ID = ID
        self.E = E              # Modulus of elasticity (Pa)
        self.St = St            # Ultimate tensile strength (minimized) (Pa)
        self.t_shear = t_shear  # Approximate shear length (pa)
        self.C = C              # Finished ring hardness (Rockwell 'C' Scale) (HRC)
        self.T = T              # Maximum recommended service temperature (°C)
        
stainless_steel = CirclipsMaterial(193000.*10**6, 1365.*10**6, 827.*10**6, 30., 205., 'ASTM-A313', 'Stainless Steel')
OT_chrome_silicon_alloy_steel = CirclipsMaterial(207000.*10**6, 1793.*10**6, 1138.*10**6, 53., 250., 'ASTM-A401', 'Oil Tempered Chrome Silicon Alloy Steel')
OT_carbon_steel = CirclipsMaterial(207000.*10**6, 1379.*10**6, 2034.*10**6, 48., 160., 'ASTM-A229', 'Oil Tempered Carbon Steel')
HD_carbon_steel = CirclipsMaterial(207000.*10**6, 1346.*10**6, 965.*10**6, 45., 150., 'ASTM-A227', 'Hard Drawn Carbon Steel')

circlips_materials = [stainless_steel, OT_chrome_silicon_alloy_steel, OT_carbon_steel, HD_carbon_steel]





class GrooveMaterial(DessiaObject):
    def __init__(self, Re, t_yield, C, name='',):
        self.Re = Re            # Ultimate tensile strength (Pa)
        self.t_yield = t_yield  # Yield stength (Pa)
        self.C = C              # Hardness (Brinell Scale) (Bhn)
        self.name = name
        
LM_carbon_steel = GrooveMaterial(462.*10**6, 310.*10**6, 150., 'Low-Mild Carbon Steel')
H_carbon_steel = GrooveMaterial(1241.*10**6, 1117.*10**6, None, 'Hardened Carbon Steel')
cast_steel = GrooveMaterial(690.*10**6, 552.*10**6, 360., 'Cast Steel')
grey_iron = GrooveMaterial(345.*10**6, None, 260., 'Grey Iron')
ductile_iron = GrooveMaterial(517.*10**6, 345.*10**6, 225., 'Ductile Iron')
cast_aluminium = GrooveMaterial(221.*10**6, 165.*10**6, 120., 'Cast Aluminium')

groove_materials = [cast_aluminium, LM_carbon_steel, H_carbon_steel, cast_steel, grey_iron, ductile_iron]


class Groove(DessiaObject):
    
    def load_capacity(self, d):
        # d the nominal diameter of the circlip ie. the shaft diameter
        Fn = self.material.Re * math.pi/4*npy.abs(d**2-self.D**2) / self.StressingCoeff(d)
        return Fn
    
    def StressingCoeff(self, d):
        h = abs(d - self.D)/2
        X = []
        Y = []
        for i in range(len(q_fct_wh['data'])):
            X.append(q_fct_wh['data'][i][0])
            Y.append(q_fct_wh['data'][i][1])
        f_q_fct_wh = interpolate.interp1d(X, Y, kind='linear', fill_value='extrapolate')
        q = f_q_fct_wh(self.w/h)
        print('q', q)
        return q


class Circlips(DessiaObject):
    _standalone_in_db = True
    
    p1 = vm.Point2D(0, 0)
    
    def RadialThickness(self, theta):
        if theta >= 0 and theta <= 2*math.pi:
            x = 1-npy.cos(theta)
            return self.b * x**(1/3)
        else : return 'The angle should be between 0 and 2*pi'
        
    
    def EllipsePoints(self, nb_points, direction):
        """
        For deporting the ellipse towards the interior of the circle, choose a negative value for 'direction'
        For deporting the ellipse towards the extorior of the circle, choose a positive value for 'direction'
        """
        x = vm.Point2D((0,-1))
        vecteurs_rot = []
        pts_ellipse = []
        for n in range(nb_points):
            vecteurs_rot.append(x.Rotation(Circlips.p1, (n+1)*2*math.pi/nb_points))
            pts_ellipse.append(vecteurs_rot[n]*(self.d/2 + npy.sign(direction)*Circlips.RadialThickness(self, (n+1)*2*math.pi/nb_points)))
        return pts_ellipse
    
#    def length(v):
#        return math.sqrt(v[0]**2+v[1]**2)
#    def dot_product(v,w):
#        return v[0]*w[0]+v[1]*w[1]
#    def determinant(v,w):
#        return v[0]*w[1]-v[1]*w[0]
#    def inner_angle(v,w):
#        cosx=Circlips.dot_product(v,w)/(Circlips.length(v)*Circlips.length(w))
#        rad=math.acos(cosx) # in radians
#        return rad*180/math.pi # returns degrees
#    def angle_counter_clockwise(A, B):
#        inner=Circlips.inner_angle(A,B)
#        det = Circlips.determinant(A,B)
#        if det<0: #this is a property of the det. If the det < 0 then B is clockwise of A
#            return 360-inner
#        else: # if the det > 0 then A is immediately clockwise of B
#            return inner
    
    
    
    def ExternalContour(self, direction):
        phaut = vm.Point2D((self.d/2, 0))
        
#        p6 = vm.Point2D((min(-6*self.d0/2,-self.d/2), -math.sqrt((self.d/2+npy.sign(direction)*self.a)**2-(6*self.d0/2)**2)))
#        p7 = vm.Point2D((min( 6*self.d0/2, self.d/2), -math.sqrt((self.d/2+npy.sign(direction)*self.a)**2-(6*self.d0/2)**2)))
#        p6 = vm.Point2D((max(-6*self.d0/2,-self.d/2), -math.sqrt((self.d/2)**2-self.d0**2/4)-npy.sign(direction)*self.a/2))
#        p7 = vm.Point2D((min( 6*self.d0/2, self.d/2), -math.sqrt((self.d/2)**2-self.d0**2/4)-npy.sign(direction)*self.a/2))
        
        
        if self.d <= 0.009:
            
            p6 = vm.Point2D((max(-6*self.d0/2,-self.d/2), -math.sqrt((self.d/2)**2-self.d0**2/4)-npy.sign(direction)*self.a/2))
            p7 = vm.Point2D((min( 6*self.d0/2, self.d/2), -math.sqrt((self.d/2)**2-self.d0**2/4)-npy.sign(direction)*self.a/2))
            
            p2 = vm.Point2D((-self.d0/2, -math.sqrt((self.d/2)**2-self.d0**2/4)))
            p3 = vm.Point2D(( self.d0/2, -math.sqrt((self.d/2)**2-self.d0**2/4)))
            p4 = vm.Point2D((-self.d0/2, -math.sqrt((self.d/2)**2-self.d0**2/4)-npy.sign(direction)*self.a))
            p5 = vm.Point2D(( self.d0/2, -math.sqrt((self.d/2)**2-self.d0**2/4)-npy.sign(direction)*self.a))
            
            pp2_arc = vm.Point2D.MiddlePoint(p2, p4)
            pp1_arc = vm.Point2D.MiddlePoint(p3, p5)
            p2_arc = vm.Vector2D(pp2_arc.vector)
            p1_arc = vm.Vector2D(pp1_arc.vector)
            p3_arc = p1_arc.Translation(vm.Vector2D((0, self.d0/2)), copy=True)
            p4_arc = p2_arc.Translation(vm.Vector2D((0, self.d0/2)), copy=True)
            p5_arc = p1_arc.Translation(vm.Vector2D((self.d0/2, 0)), copy=True)
            p6_arc = p2_arc.Translation(vm.Vector2D((-self.d0/2, 0)), copy=True)
            p7_arc = p1_arc.Translation(vm.Vector2D((0, -self.d0/2)), copy=True)
            p8_arc = p2_arc.Translation(vm.Vector2D((0, -self.d0/2)), copy=True)
            
            p8 = p6
            p9 = p7
            
            l1 = vm.LineSegment2D(p3, p3_arc)
            l2 = vm.LineSegment2D(p2, p4_arc)
            
            arc1 = vm.Arc2D(p3_arc, p5_arc, p7_arc)
            arc2 = vm.Arc2D(p4_arc, p6_arc, p8_arc)
            
            
        
        if self.d > 0.009:
            
            p6 = vm.Point2D((max(-6*self.d0/2,-self.d/2), -math.sqrt((self.d/2+npy.sign(direction)*self.a)**2-(6*self.d0/2)**2)))
            p7 = vm.Point2D((min( 6*self.d0/2, self.d/2), -math.sqrt((self.d/2+npy.sign(direction)*self.a)**2-(6*self.d0/2)**2)))
            
            p2 = vm.Point2D((-self.d0, -math.sqrt((self.d/2)**2-self.d0**2)))
            p3 = vm.Point2D(( self.d0, -math.sqrt((self.d/2)**2-self.d0**2)))
            p4 = vm.Point2D((-self.d0, -math.sqrt((self.d/2)**2-self.d0**2)-npy.sign(direction)*self.a))
            p5 = vm.Point2D(( self.d0, -math.sqrt((self.d/2)**2-self.d0**2)-npy.sign(direction)*self.a))
            
            p8 = vm.Point2D((-9*self.d0/2, -math.sqrt((self.d/2+npy.sign(direction)*self.a)**2-(9*self.d0/2)**2)))
            p9 = vm.Point2D(( 9*self.d0/2, -math.sqrt((self.d/2+npy.sign(direction)*self.a)**2-(9*self.d0/2)**2)))



        # Resolution of the external ellipse
        nb_points = 20
        
        pts_ellipse = Circlips.EllipsePoints(self, nb_points, direction)
                 
        points_to_del = []
        v0 = vm.Vector2D((0,-1))
        v9 = vm.Vector2D((p9[0],p9[1]))
        v8 = vm.Vector2D((p8[0],p8[1]))
        angle_p9 = AngleTrigo(v0,v9)
        angle_p8 = AngleTrigo(v0,v8)
        for n in range(len(pts_ellipse)):
            angle_p_ellipse = AngleTrigo(v0,vm.Vector2D((pts_ellipse[n][0],pts_ellipse[n][1])))
            if angle_p_ellipse < angle_p9 or angle_p_ellipse > angle_p8 :
                points_to_del.append(n)
                
        points_to_del.reverse()
        for point in points_to_del:
            del pts_ellipse[point]
            
        # Construction of the inner circle
        arcinterieur = vm.Arc2D(p2, phaut, p3)
                
        if self.d <= 0.009:
            raccord = []
            raccord.append(p7_arc)
            raccord.append(p5)
            raccord.append(p9)
            for i in pts_ellipse: raccord.append(i)
            raccord.append(p8)
            raccord.append(p4)
            raccord.append(p8_arc)
            
            # Radius of curvature for the Rounded Line
            radius = {}
            r = self.d/10
            radius[2] = r
            radius[3] = r
            radius[len(raccord)-4] = r
            radius[len(raccord)-3] = r
            for n in range(4, len(raccord)-4):   
                radius[n] = self.d/2
            
            
            racc = vm.primitives2d.RoundedLineSegments2D(raccord, radius, adapt_radius=True)

            contour = vm.Contour2D([arcinterieur, l1, l2, arc1, arc2, racc])
            
        if self.d > 0.009:
            raccord = []
            raccord.append(p3)
            raccord.append(p5)
            raccord.append(p7)
            raccord.append(p9)
            for i in pts_ellipse: raccord.append(i)
            raccord.append(p8)
            raccord.append(p6)
            raccord.append(p4)
            raccord.append(p2)
    
            # Radius of curvature for the Rounded Line
            radius = {}
            r = self.d/10
            radius[2] = r
            radius[3] = r
            radius[4] = r
            radius[len(raccord)-5] = r
            radius[len(raccord)-4] = r
            radius[len(raccord)-3] = r
            for n in range(5, len(raccord)-5):   
                radius[n] = self.d/2
            
            
            racc = vm.primitives2d.RoundedLineSegments2D(raccord, radius, adapt_radius=True)

            contour = vm.Contour2D([arcinterieur,racc])  

        return contour
        
    def InternalContour(self, direction):
        if self.d <= 0.009:
            return None
#        p1 = vm.Point2D(( 2*self.d0, - math.sqrt((self.d/2)**2 - npy.sign(direction)*self.d0**2/2) - npy.sign(direction)*self.a/2))
#        p2 = vm.Point2D((-2*self.d0, - math.sqrt((self.d/2)**2 - npy.sign(direction)*self.d0**2/2) - npy.sign(direction)*self.a/2))
        p2 = vm.Point2D((-self.d0, -math.sqrt((self.d/2)**2-self.d0**2)))
        p3 = vm.Point2D(( self.d0, -math.sqrt((self.d/2)**2-self.d0**2)))
        p6 = vm.Point2D((-6*self.d0/2, -math.sqrt((self.d/2+npy.sign(direction)*self.a)**2-(6*self.d0/2)**2)))
        p7 = vm.Point2D(( 6*self.d0/2, -math.sqrt((self.d/2+npy.sign(direction)*self.a)**2-(6*self.d0/2)**2)))
        
        p1 = vm.Point2D.MiddlePoint(p2, p6)
        p2 = vm.Point2D.MiddlePoint(p3, p7)
        
        c1 = vm.Circle2D(p1, self.d0/2)
        c2 = vm.Circle2D(p2, self.d0/2)
        contour1 = vm.Contour2D([c1])
        contour2 = vm.Contour2D([c2])
        return [contour1, contour2]
    
    def Plot(self):
        contour = self.ExternalContour()
        f, a = contour.MPLPlot()
        if self.InternalContour() != None:
            c1, c2 = self.InternalContour()
            c1.MPLPlot(ax=a)
            c2.MPLPlot(ax=a)
    
    def volmdlr_primitives(self, center = vm.O3D, x = vm.X3D, y = vm.Y3D):        
        z = x.Cross(y)
        z.Normalize()
        profile = vm3d.ExtrudedProfile(center, x, y, self.ExternalContour(), self.InternalContour(), self.e*z)    
        return [profile]
        
        
    def PermissibleDishingAngle(self, d1):
        X = []
        Y = []
        for i in range(len(phi_fct_d1['data'])):
            X.append(phi_fct_d1['data'][i][0])
            Y.append(phi_fct_d1['data'][i][1])
        f_phi_fct_d1 = interpolate.interp1d(X, Y, kind='linear', fill_value='extrapolate')
        phi = f_phi_fct_d1(d1)
        return phi
    
    def Groove(self, F, w_max, direction, groove_material=groove_materials[0]):
        # circlip is an external circlip object
        # direction should be a positive value for external circlips
        #           and a negative value for internal circlips
        # F the axial force in the bearing
        # w_max the maximum width available
        if direction > 0:
            D = 0.972 * self.d - 0.745/1000
            print('ext')
        elif direction < 0:
            D = 1.026 * self.d + 0.916/1000
            print('int')
        else:
            return None
        x = 1.143 * self.e
        h = abs(self.d - D)/2
        
        print('h', h)
        print('D', D)
        print(4*h,'?<=?', w_max)
        
        if 4*h <= w_max:
            w = 4*h    
        else:
            w = w_max
            
        print('w/h', w/h)
                        
        if direction > 0:
            groove = ExternalGroove(D, x, w, groove_material)
            print('load capa ext : ', groove.load_capacity(self.d))
            if groove.load_capacity(self.d) < F:
                return None # The circlip and the groove won't hold
        elif direction < 0:
            groove = InternalGroove(D, x, w, groove_material)
            print('load capa int : ', groove.load_capacity(self.d))
            if groove.load_capacity(self.d) < F:
                return None # The circlip and the groove won't hold
            return groove
        else:
            return None

        return groove
    
class ExternalCirclips(Circlips):
    """
    :param material: a material object
    """
    def __init__(self, d:float, b:float, a:float, e:float, d0:float,
                 material:CirclipsMaterial=circlips_materials[0], name:str=''):
        self.d = d                  # Inner diameter
        self.b = b                  # Radial width
        self.a = a                  # Cutoff length
        self.e = e                  # Thickness
        self.d0 = d0                # Tool holes diameter
        self.material = material 
        self.name = name
        
    def groove(self, F:float, w_max:float):
        return Circlips.Groove(self, F, w_max, 1)
    
    def ExternalContour(self):
        return Circlips.ExternalContour(self, direction=1)
              
    def InternalContour(self):        
        return Circlips.InternalContour(self, direction=1)  
    
    
    
class InternalCirclips(Circlips):  
    """
    :param material: a material object
    """
    def __init__(self, d:float, b:float, a:float, e:float, d0:float,
                 material:CirclipsMaterial=circlips_materials[0], name:str=''):
        self.d = d                  # Outter diameter
        self.b = b                  # Radial width
        self.a = a                  # Cutoff length
        self.e = e                  # Thickness
        self.d0 = d0                # Tool holes diameter
        self.name = name
        
    def Groove(self, F, w_max):
        return Circlips.Groove(self, F, w_max, -1)
    
    def ExternalContour(self):
        return Circlips.ExternalContour(self, direction=-1)
    
    def InternalContour(self):        
        return Circlips.InternalContour(self, direction=-1)
    



class ExternalCirclipsCatalog(DessiaObject):
    _standalone_in_db = True
    def __init__(self, circlips:List[ExternalCirclips], name=''):
         self.circlips = circlips    # List of circlips
         
    def ExternalShaftToCirclipGroove(self, Dshaft, F, w_max):
        # Dshaft the shaft outter diameter
        # F the axial force in the bearings
        # w_max the maximum width available of the shoulder
        circlip_groove = [] # a list of tuple (circlip, groove) possible        
        matching_circlips = [circlip for circlip in self.circlips if circlip.d == Dshaft]
        for circlip in matching_circlips:
            groove = circlip.Groove(F, w_max)
            new_pair = (circlip, groove)
            circlip_groove.append(new_pair)
        return circlip_groove
        



class InternalCirclipsCatalog(DessiaObject):
    _standalone_in_db = True
    
    def __init__(self, circlips:List[InternalCirclips], name=''):
         self.circlips = circlips    # List of circlips
         
    def InternalShaftToCirclipGroove(self, Dshaft, F, w_max):
        # Dshaft the shaft inner diameter
        # F the axial force in the bearings
        # w_max the maximum width available of the shoulder
        circlip_groove = [] # a list of tuple (circlip, groove) possible        
        matching_circlips = [circlip for circlip in self.circlips if circlip.d == Dshaft]
        for circlip in matching_circlips:
            groove = circlip.Groove(F, w_max)
            new_pair = (circlip, groove)
            circlip_groove.append(new_pair)
        return circlip_groove
        
    

class ExternalGroove(Groove):
    """
    :param material: a material object
    """
    def __init__(self, D:float, x:float, w:float,
                 material:CirclipsMaterial=groove_materials[0], name:str=''):
        self.D = D                  # Groove diameter 
        self.x = x                  # Groove width
        self.w = w                  # Shoulder width
        self.material = material    # Groove material
        self.name = name
    
    def load_capacity(self, d):
        return Groove.load_capacity(self, d)
        
    def DetachingSpeed(self):
        n_abl = 37200000*self.b/(self.d2+self.b)**2 * (self.d2-self.d3/(self.d3+self.b))**0.5
        return n_abl





class InternalGroove(Groove):
    def __init__(self, D, x, w, material=groove_materials[0]):
        self.D = D
        self.x = x
        self.w = w
        self.material = material
        
    def load_capacity(self, d):
        return Groove.load_capacity(self, d)
    
    
with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                                   'mechanical_components/catalogs/roymech_external_circlips.json') as roymech_external_circlips_json:
    EXTERNAL_CIRCLIPS_CATALOG = ExternalCirclipsCatalog.load_from_file(
        roymech_external_circlips_json)

with pkg_resources.resource_stream(pkg_resources.Requirement('mechanical_components'),
                                   'mechanical_components/catalogs/roymech_internal_circlips.json') as roymech_internal_circlips_json:
    INTERNAL_CIRCLIPS_CATALOG = InternalCirclipsCatalog.load_from_file(
        roymech_internal_circlips_json)






        
        
        
    
    
    
    

        
        
    
