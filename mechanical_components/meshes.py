#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 17 02:13:01 2018

@author: Pierre-Emmanuel Dumouchel
"""

import numpy as npy
from scipy import interpolate
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm
from scipy.optimize import minimize,fsolve
import itertools
import networkx as nx

import powertransmission.tools as tools
import matplotlib.pyplot as plt
import matplotlib.patches as patches

#data_coeff_YB_Iso
evol_coeff_yb_iso={'data':[[0.0,1.0029325508401201],
                           [4.701492563229561,0.9310850480431024],
                           [9.104477651442416,0.8782991233021732],
                           [14.5522388104227,0.8240469255458759],
                           [19.328358165913905,0.784457481990179],
                           [23.955224059269884,0.7609970656504502],
                           [28.507462609617956,0.7521994155347822],
                           [34.029850499665855,0.7507331425194141],
                           [40.0,0.7492668574805859]
                          ], 'x':'Linear','y':'Linear'}

#data_wholer_curve
wholer_hardened_alloy_steel={'data':[[4.296196199237153,1.9797762011105589],
                                     [4.824840106199563,1.9413306094362142],
                                     [5.3344338175705674,1.908892154601565],
                                     [6.115493253679078,1.8632380197445122],
                                     [6.596511629990596,1.8560294765618042],
                                     [7.144205815889171,1.8536266428508523],
                                     [7.691899918442984,1.8524252154829133],
                                     [8.010991340520903,1.8524252154829133]
                                    ], 'x':'Log','y':'Log'}

wholer_nitrided_alloy_steel={'data':[[4.104865629699472,1.9252942042661974],
                                     [4.568697315952783,1.8521640228225367],
                                     [4.887581626297173,1.8046294185503593],
                                     [5.438381821440599,1.7900033864666123],
                                     [6.402282079596832,1.7918316299646175],
                                     [7.264719174616821,1.7918316299646175],
                                     [7.989456220850952,1.793659894487549]
                                    ], 'x':'Log','y':'Log'}

wholer_through_hardened_steel={'data':[[4.172369719531124,1.895676495604088],
                                       [4.677200861168087,1.7983611100752137],
                                       [4.9677168648417585,1.741894170956562],
                                       [5.329671247836526,1.6842258044699714],
                                       [5.439210101685194,1.672211551815507],
                                       [6.091680488353632,1.6734129791834462],
                                       [7.139443246155129,1.671010124447568],
                                       [8.00146620105282,1.6758158339193243]
                                      ],'x':'Log','y':'Log'}

wholer_surface_hardened_steel={'data':[[4.281908490035029,1.7611169667937343],
                                       [4.701013626493532,1.6998443182033265],
                                       [5.015342395492649,1.6553916107142128],
                                       [5.358246582896013,1.6109389032250994],
                                       [5.620187251510196,1.5857089915731581],
                                       [6.020242109032534,1.5748961873115592],
                                       [6.567936294931109,1.5748961873115592],
                                       [7.263269725861159,1.5772990210225108],
                                       [7.996703631318779,1.5772990210225108]
                                      ], 'x':'Log','y':'Log'}

wholer_carbon_steel={'data':[[4.307791955971963,1.6419147590563592],
                             [5.242702822291173,1.535876005424268],
                             [5.938450393343521,1.4700588400224806],
                             [6.518240063668731,1.431665495290182],
                             [7.221234961844144,1.4334937598131132],
                             [7.989456220850952,1.4353220033111185]
                            ], 'x':'Log','y':'Log'}

wholer_cast_iron={'data':[[4.307791955971963,1.6419147590563592],
                          [5.242702822291173,1.535876005424268],
                          [5.938450393343521,1.4700588400224806],
                          [6.518240063668731,1.431665495290182],
                          [7.221234961844144,1.4334937598131132],
                          [7.989456220850952,1.4353220033111185]
                         ], 'x':'Log','y':'Log'}

wholer_bronze={'data':[[4.307791955971963,1.6419147590563592],
                       [5.242702822291173,1.535876005424268],
                       [5.938450393343521,1.4700588400224806],
                       [6.518240063668731,1.431665495290182],
                       [7.221234961844144,1.4334937598131132],
                       [7.989456220850952,1.4353220033111185]
                      ], 'x':'Log','y':'Log'}

wholer_grey_iron={'data':[[4.307791955971963,1.6419147590563592],
                          [5.242702822291173,1.535876005424268],
                          [5.938450393343521,1.4700588400224806],
                          [6.518240063668731,1.431665495290182],
                          [7.221234961844144,1.4334937598131132],
                          [7.989456220850952,1.4353220033111185]
                         ], 'x':'Log','y':'Log'}
#data_gear_material

sigma_hardened_alloy_steel={'data':[[1.8422104370714443,1.4645831828946267],
                                    [1.948612010770208,1.5219116983411152],
                                    [2.0605171321606295,1.5810895335609718]
                                    ,[2.141235568740199,1.6254729099758645]
                                   ], 'x':'Log','y':'Log'}

sigma_nitrided_alloy_steel={'data':[[1.8458794622934307,1.4349942652846983],
                                    [1.943108482795906,1.488624180937243],
                                    [2.0201578941534892,1.5274596179084272],
                                    [2.128393990321924,1.5866374531282839]
                                   ],'x':'Log','y':'Log'}

sigma_through_hardened_steel={'data':[[1.7798371068844516,1.292597616678765],
                                      [1.921094370898698,1.3850629693024938],
                                      [2.032999472571764,1.4627338829976548],
                                      [2.1650841833897223,1.5533499158480155]
                                     ],'x':'Log','y':'Log'}

sigma_surface_hardened_steel={'data':[[1.8312033811228403,1.115064130895591],
                                      [1.932101426847302,1.200132264055036],
                                      [2.038503000546066,1.2852003773380847]
                                     ], 'x':'Log','y':'Log'}

sigma_carbon_steel={'data':[[1.677104538690319,1.1002696720906269],
                            [1.7633265032441903,1.1723926463420797],
                            [1.8385414118494579,1.2389677010262203],
                            [1.8844041581135444,1.2796524577707729]
                           ], 'x':'Log','y':'Log'}

sigma_cast_iron={'data':[[1.4734739247717241,0.922736186307453],
                         [1.5468543306246763,0.9837633214242817],
                         [1.6073931580593532,1.0336946174064863],
                         [1.6404143456225206,1.0688314545837265]
                        ], 'x':'Log','y':'Log'}

sigma_bronze={'data':[[1.313871566195314,0.7858874572688317],
                      [1.3890864826875238,0.8487638922826322],
                      [1.4294457009773085,0.8802021097895326],
                      [1.4551288380965028,0.9097910273994609]
                     ], 'x':'Log','y':'Log'}

sigma_grey_iron={'data':[[1.354230792372041,0.7100658633470387],
                         [1.4276111785076375,0.7766409180311793],
                         [1.4936535339166166,0.84691459238566],
                         [1.5431853054026896,0.8986951882648367],
                         [1.5725374677438706,0.933832025442077]
                        ], 'x':'Log','y':'Log'}

class Material:
    """
    Gear material
    
    :param volumic_mass: A float to define the gear volumic mass
    :param data_coeff_YB_Iso: a dictionary to define the YB parameter of the ISO description
    :param data_wholer_curve: a dictionary to define the wholer slope of the ISO description
    :param data_gear_material: a dictionary to define the maximum gear stress
    
    :data_coeff_YB_Iso: - **'data'** matrix define points of the YB curve in the plane (YB, helix_angle)
        - **'x'** string define the x axis evolution ('Log' or 'Linear')
        - **'y'** string define the y axis evolution ('Log' or 'Linear') 
        
    :data_wholer_curve: - **'data'** matrix define points of the wholer slope in the plane (wholer slope, number of cycle)
        - **'x'** string define the x axis evolution ('Log' or 'Linear')
        - **'y'** string define the y axis evolution ('Log' or 'Linear')
        
    :data_gear_material: - **'data'** matrix define points of the maximum gear stress (maximum gear stress, wholer slope)
        - **'x'** string define the x axis evolution ('Log' or 'Linear')
        - **'y'** string define the y axis evolution ('Log' or 'Linear')
        
    >>> volumic_mass=7800
    >>> data_coeff_YB_Iso={'data':[[0.0,1.0029325508401201],
                           [4.701492563229561,0.9310850480431024],
                           [23.955224059269884,0.7609970656504502],
                           [40.0,0.7492668574805859]
                          ], 'x':'Linear','y':'Linear'}
    >>> data_wholer_curve={'data':[[4.307791955971963,1.6419147590563592],
                       [6.518240063668731,1.431665495290182],
                       [7.989456220850952,1.4353220033111185]
                      ], 'x':'Log','y':'Log'}
    >>> data_gear_material={'data':[[1.313871566195314,0.7858874572688317],
                      [1.4294457009773085,0.8802021097895326],
                      [1.4551288380965028,0.9097910273994609]
                     ], 'x':'Log','y':'Log'}
    >>> material1=Material(volumic_mass, data_coeff_YB_Iso, 
                           data_wholer_curve, data_gear_material)
    """
    def __init__(self,volumic_mass, data_coeff_YB_Iso, data_wholer_curve,
                 data_gear_material):
        self.volumic_mass = volumic_mass
        self.data_coeff_YB_Iso = data_coeff_YB_Iso
        self.data_wholer_curve = data_wholer_curve
        self.data_gear_material = data_gear_material
    
    def FunCoeff(self,x,data,type_x='Linear',type_y='Linear'):
        """ Interpolation of material data
        
        :param x: value of the interpolation
        :param data: dictionary of the input data
        :param type_x: type of the x axis of the data matrix ('Log' or 'Linear')
        :param type_y: type of the y axis of the data matrix ('Log' or 'Linear')
        
        :returns:  interpolation value
        
        >>> interp1=material1.FunCoeff(x = 5.2,data = data_wholer_curve,
                                       type_x = 'Log',type_y = 'Log') 
        """
        if type_x=='Log': 
            x=npy.log10(x)
        f = interpolate.interp1d(list(data[:,0]),list(data[:,1]), fill_value='extrapolate')
        sol=float(f(x))
        if type_y=='Log':
            sol=10**sol
        return sol
    
    def Dict(self):
        """Export dictionary
        
        :returns:  dictionary with all material parameters
        """
        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv==npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=float(v)
            else:
                d[k]=v
        return d
    
hardened_alloy_steel=Material(7850, evol_coeff_yb_iso,
                              wholer_hardened_alloy_steel,
                              sigma_hardened_alloy_steel)

nitrided_alloy_steel=Material(7850, evol_coeff_yb_iso, wholer_nitrided_alloy_steel,
                              sigma_nitrided_alloy_steel)

through_hardened_steel=Material(7850, evol_coeff_yb_iso,
                                wholer_through_hardened_steel,
                                sigma_through_hardened_steel)

surface_hardened_steel=Material(7850, evol_coeff_yb_iso,
                                wholer_surface_hardened_steel,
                                sigma_surface_hardened_steel)

carbon_steel=Material(7850, evol_coeff_yb_iso, wholer_carbon_steel, sigma_carbon_steel)

cast_iron=Material(7200, evol_coeff_yb_iso, wholer_cast_iron, sigma_cast_iron)

bronze=Material(8200, evol_coeff_yb_iso, wholer_bronze, sigma_bronze)

grey_iron=Material(7200, evol_coeff_yb_iso, wholer_grey_iron, sigma_grey_iron)

class Rack:
    """
    Gear rack definition
    
    :param transverse_pressure_angle: definition of the transverse pressure angle of the rack
    :type transverse_pressure_angle: radian
    
    >>> Rack1=Rack(20/180.*npy.pi) #definition of an ISO rack
    """
    def __init__(self,transverse_pressure_angle):
        self.transverse_pressure_angle=transverse_pressure_angle

    
    def RackParam(self,transverse_pressure_angle,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):
        
        self.transverse_pressure_angle=transverse_pressure_angle
        self.transverse_radial_pitch=self.module*math.pi
        self.gear_addendum=coeff_gear_addendum*self.module
        self.gear_dedendum=coeff_gear_dedendum*self.module
        self.root_radius=coeff_root_radius*self.module
        self.circular_tooth_thickness=coeff_circular_tooth_thickness*self.transverse_radial_pitch

        self.tooth_space=self.transverse_radial_pitch-self.circular_tooth_thickness
        self.whole_depth=self.gear_addendum+self.gear_dedendum
        self.clearance=self.root_radius-self.root_radius*npy.sin(self.transverse_pressure_angle)

        # trochoide parameter
        self.a=self.tooth_space/2-self.gear_dedendum*npy.tan(self.transverse_pressure_angle)-self.root_radius*npy.tan(1/2*npy.arctan(npy.cos(self.transverse_pressure_angle)/(npy.sin(self.transverse_pressure_angle))))
        self.b=self.gear_dedendum-self.root_radius

    def Update(self,module,transverse_pressure_angle=None,coeff_gear_addendum=1,
               coeff_gear_dedendum=1.25,coeff_root_radius=0.38,coeff_circular_tooth_thickness=0.5):
        """
        Update of the gear rack
        
        :param module: update of the module of the rack define on the pitch factory diameter
        :type module: m
        :param transverse_pressure_angle: update of the transverse pressure angle of the rack
        :type transverse_pressure_angle: radian
        :param coeff_gear_addendum: update of the gear addendum coefficient (gear_addendum = coeff_gear_addendum*module)
        :param coeff_gear_dedendum: update of the gear dedendum coefficient (gear_dedendum = coeff_gear_dedendum*module) (top of the rack)
        :param coeff_root_radius: update of the root radius coefficient (root_radius = coeff_root_radius*module)
        :param coeff_circular_tooth_thickness: update of the circular tooth thickness coefficient (circular_tooth_thickness = coeff_circular_tooth_thickness*transverse_radial_pitch)
        
        >>> input={'module':2*1e-3,'transverse_pressure_angle':21/180.*npy.pi}
        >>> Rack1.Update(**input) # Update of the rack definition
        """
        if transverse_pressure_angle==None:
            transverse_pressure_angle=self.transverse_pressure_angle
        self.module=module
        self.RackParam(transverse_pressure_angle,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
        
    ### Optimization Method
    
    def CheckRackViable(self):
        """ Check the viability of the rack toward the top and the root
        
        :results: boolean variable, and a list of element to be positive for the optimizer
        """
        list_ineq=[]
        list_ineq.append(self.transverse_radial_pitch-self.circular_tooth_thickness
                         -2*self.gear_dedendum*npy.tan(self.transverse_pressure_angle)
                         -2*(self.root_radius*npy.cos(self.transverse_pressure_angle)-npy.tan(self.transverse_pressure_angle)
                         *self.root_radius*(1-npy.sin(self.transverse_pressure_angle))))
        list_ineq.append(self.circular_tooth_thickness-2*(self.gear_addendum*npy.tan(self.transverse_pressure_angle)))
        check=False
        if min(list_ineq)>0:
            check=True
        return check,list_ineq
    
    def ListeIneq(self):
        """ Compilation method for inequality list used by the optimizer
        
        :results: vector of data that should be positive
        """
        check,ineq=self.CheckRackViable
        return ineq

    def Dict(self):
        """Export dictionary
        
        :returns:  dictionary with all rack parameters
        
        >>> print(Rack1.Dict())
        """
        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv==npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=float(v)
            else:
                d[k]=v
        return d
    
    def Contour(self,number_pattern):
        """ Construction of the volmdr 2D rack profile
        
        :param number_pattern: number of rack pattern to define
        """
        p1=vm.Point2D((0,0))
        p2=p1.Translation((self.gear_addendum*npy.tan(self.transverse_pressure_angle),self.gear_addendum))
        p4=p1.Translation((self.circular_tooth_thickness,0))
        p3=p4.Translation((-self.gear_addendum*npy.tan(self.transverse_pressure_angle),self.gear_addendum))
        p5=p4.Translation((self.gear_dedendum*npy.tan(self.transverse_pressure_angle),-self.gear_dedendum))
        p7=p4.Translation((self.tooth_space,0))
        p6=p7.Translation((-self.gear_dedendum*npy.tan(self.transverse_pressure_angle),-self.gear_dedendum))
        L=primitives2D.RoundedLines2D([p1,p2,p3,p4,p5,p6,p7],{4:self.root_radius,5:self.root_radius},False)
        Rack_Elem=[]
        for i in range(number_pattern):
            Rack_Elem.append(L.Translation(((i)*(p7.vector-p1.vector))))
        p10=Rack_Elem[0].points[0]
        p15=Rack_Elem[-1].points[-1]
        p11=p10.Translation((-self.circular_tooth_thickness,0))
        p12=p11.Translation((0,2*self.whole_depth))
        p14=p15.Translation((self.circular_tooth_thickness,0))
        p13=p14.Translation((0,2*self.whole_depth))
        Rack_Elem.append(primitives2D.RoundedLines2D([p10,p11,p12,p13,p14,p15],{},False))
        return Rack_Elem
    
    def Plot(self,number_pattern):
        """ Plot function of the rack
        
        :param number_pattern: number of rack pattern to draw
        """
        Rack_Elem=self.Contour(number_pattern)
        RackElem=vm.Contour2D(Rack_Elem)
        RackElem.MPLPlot()

    def CSVExport(self):
        """
        Export CSV format
        
         :returns:  list of all element in dict() function
        """
        d=self.__dict__.copy()
        return list(d.keys()),list(d.values())

class Mesh:
    """
    Gear mesh definition
    
    :param z: number of tooth
    :param db: base diameter
    :type db: m
    :param cp: coefficient profile shift of the rack
    :param transverse_pressure_angle_rack: transverse pressure angle of the rack
    :type transverse_pressure_angle_rack: radian
    :param coeff_gear_addendum: update of the gear addendum coefficient (gear_addendum = coeff_gear_addendum*module)
    :param coeff_gear_dedendum: update of the gear dedendum coefficient (gear_dedendum = coeff_gear_dedendum*module)
    :param coeff_root_radius: update of the root radius coefficient (root_radius = coeff_root_radius*module)
    :param coeff_circular_tooth_thickness: update of the circular tooth thickness coefficient (circular_tooth_thickness = coeff_circular_tooth_thickness*transverse_radial_pitch)
    :param material: class material define the gear mesh material
    :param gear_width: gear mesh width
        
    >>> input={z:13, db:40*1e-3, cp:0.3, transverse_pressure_angle_rack:20/180.*npy.pi,
                 coeff_gear_addendum:1, coeff_gear_dedendum:1, coeff_root_radius:1,
                 coeff_circular_tooth_thickness:1}
    >>> mesh1=Mesh(**input) # generation of one gear mesh
    """
    def __init__(self, z, db, cp, transverse_pressure_angle_rack,
                 coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
                 coeff_circular_tooth_thickness,material=None,gear_width=1):
        
        self.rack=Rack(transverse_pressure_angle_rack)
        self.GearParam(z,db,cp,transverse_pressure_angle_rack,
                       coeff_gear_addendum, coeff_gear_dedendum, 
                       coeff_root_radius, coeff_circular_tooth_thickness)
        
        # Definition of default parameters
        if material==None:
            self.material=hardened_alloy_steel
        self.material=material
        
        self.gear_width=gear_width
        
    def Update(self,z,db,cp,transverse_pressure_angle_rack,
               coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,
               coeff_circular_tooth_thickness,material,gear_width=1):
        """ Update of the gear mesh
        
        :param all: same parameters of this class initialisation
            
        >>> input={z:14, db:42*1e-3, cp:0.5}
        >>> mesh1.Update(**input)
        """
        self.GearParam(z,db,cp, transverse_pressure_angle_rack,
                       coeff_gear_addendum, coeff_gear_dedendum,
                       coeff_root_radius, coeff_circular_tooth_thickness)
        self.gear_width=gear_width
        
    ### geometry definition
        
    def GearParam(self,z, db, cp, transverse_pressure_angle_rack,
                  coeff_gear_addendum, coeff_gear_dedendum,
                  coeff_root_radius, coeff_circular_tooth_thickness):
        
        self.z=z
        self.db=db
        self.dff=self.db/npy.cos(transverse_pressure_angle_rack)
        module_rack=self.dff/self.z
        self.rack.Update(module_rack,transverse_pressure_angle_rack,
                         coeff_gear_addendum,coeff_gear_dedendum, 
                         coeff_root_radius,coeff_circular_tooth_thickness)
        self.coefficient_profile_shift=cp
        
        self.outside_diameter=(self.dff
                               +2*(self.rack.gear_addendum
                                   +self.rack.module*self.coefficient_profile_shift))
        self.alpha_outside_diameter = npy.arccos(self.db/self.outside_diameter)
        self.root_diameter=(self.dff
                            - 2*(self.rack.gear_dedendum
                                - self.rack.module*self.coefficient_profile_shift))
        self.root_diameter_active,self.phi_trochoide=self._RootDiameterActive()
        self.alpha_root_diameter_active=npy.arccos(self.db/self.root_diameter_active)
        
        self.alpha_pitch_diameter=npy.arccos(self.db/self.dff)
        self.circular_tooth_thickness = (self.rack.circular_tooth_thickness
                                       +self.rack.module*self.coefficient_profile_shift
                                           *npy.tan(self.rack.transverse_pressure_angle)
                                       +self.rack.module*self.coefficient_profile_shift
                                           *npy.tan(self.rack.transverse_pressure_angle))
        self.tooth_space=self.rack.transverse_radial_pitch-self.circular_tooth_thickness
        self.outside_active_angle = (2*self.circular_tooth_thickness/self.dff-2
                                     *(npy.tan(self.alpha_outside_diameter)
                                        -self.alpha_outside_diameter
                                        -npy.tan(self.alpha_pitch_diameter)
                                        +self.alpha_pitch_diameter))
        self.base_circular_tooth_thickness = (self.db/2
                                              *(2*self.circular_tooth_thickness/self.dff
                                                +2*(npy.tan(self.alpha_pitch_diameter)
                                                -self.alpha_pitch_diameter)))
        
        self.root_angle=self.tooth_space/(self.dff/2)-2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
        self.root_gear_angle=self.circular_tooth_thickness/(self.dff/2)+2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
        
    def GearSection(self,diameter):
        """ Definition of the gear section
        
        :param diameter: diameter of the gear section calculation
        :type diameter: m
        
        :results: gear section in m
            
        >>> gs=mesh1.GearSection(44*1e-3)
        """
        alpha_diameter=npy.arccos(self.db/diameter)
        theta1=(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter)-(npy.tan(alpha_diameter)-alpha_diameter)
        return diameter/2*(2*theta1+self.outside_active_angle)
        
    def _RootDiameterActive(self):
        a=self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.dff/2
        phi=-(a+b*npy.tan(npy.pi/2-self.rack.transverse_pressure_angle))/r
        root_diameter_active=2*norm(self._Trochoide(phi))
        return root_diameter_active,phi
    
    ### Optimization Method
    
    def ListeIneq(self):
        """ Compilation method for inequality list used by the optimizer
        
        :results: vector of data that should be positive
        """
        check,ineq=self.rack.CheckRackViable()
        return ineq
    
    ### Trace method
    
    def Contour(self,discret=10,list_number=[None]):
        """ Definition of the gear contour for volmdlr
        
        :param discret: number of discretization points on the gear mesh involute
        :param list_number: list of gear tooth to include on the graph
        
        :results: volmdlr profile
            
        >>> C1=mesh1.Contour(10)
        >>> G1=vm.Contour2D(C1)
        >>> G1.MPLPlot() # generate a plot with matplotlib
        """
        # Analytical tooth profil
        if list_number==[None]:
            list_number=npy.arange(int(self.z))
        L=[self._OutsideTrace(0)]
        L.append(self._InvoluteTrace(discret,0,'T'))
        L.append(self._TrochoideTrace(2*discret,0,'T'))
        L.append(self._RootCircleTrace(0))
        L.append(self._TrochoideTrace(2*discret,0,'R'))
        L.append(self._InvoluteTrace(discret,1,'R'))
        for i in list_number[1::]:
            L.append(self._OutsideTrace(i))
            L.append(self._InvoluteTrace(discret,i,'T'))
            L.append(self._TrochoideTrace(2*discret,i,'T'))
            L.append(self._RootCircleTrace(i))
            L.append(self._TrochoideTrace(2*discret,i,'R'))
            L.append(self._InvoluteTrace(discret,i+1,'R'))
        return L
    
    def _InvoluteTrace(self,discret,number,ind='T'):
        
        if ind=='T':
            drap=1
            theta=npy.linspace(npy.tan(self.alpha_outside_diameter),
                               npy.tan(self.alpha_root_diameter_active),discret)
        else:
            drap=-1
            theta=npy.linspace(npy.tan(self.alpha_root_diameter_active),
                               npy.tan(self.alpha_outside_diameter),discret)
        
        sol=self._Involute(drap*theta)
        x=sol[0]
        y=sol[1]
        p=[vm.Point2D((x[0],y[0]))]
        for i in range(1,discret):
            p.append(vm.Point2D((x[i],y[i])))
        ref=primitives2D.RoundedLines2D(p,{},False)
        if ind=='T':
            L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.z)
            self.rac=L.points[-1]
        else:
            L=ref.Rotation(vm.Point2D((0,0)),
                           self.base_circular_tooth_thickness*2/self.db)
            L=L.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.z)
            L.points[0]=self.rac
        return L
    
    def _TrochoideTrace(self, discret, number, type_flank='T'):
        # Function evolution of the trochoide
        if type_flank=='T':
            indice_flank=1
        else:
            indice_flank=-1
        
        a=indice_flank*self.rack.a # indice a in the ISO definition of the rack
        phi0=a/(self.dff/2)
        
        list_2D=[]
        if type_flank=='R':
            theta=npy.linspace(phi0,indice_flank*self.phi_trochoide,discret)
        else:
            theta=npy.linspace(indice_flank*self.phi_trochoide,phi0,discret)
        for t in theta:
            list_2D.append(vm.Point2D((self._Trochoide(t,type_flank))))
        list_2D=primitives2D.RoundedLines2D(list_2D,{},False)
        list_2D=list_2D.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        if type_flank=='T':
            export_2D=list_2D.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.z)
            export_2D.points[0]=self.rac
        else:
            export_2D=list_2D.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.z)
            self.rac=export_2D.points[-1]
        return export_2D
    
    def _RootCircleTrace(self,number):
        # 2D trace of the connection between the two trochoide
        
        # on the drive flank
        indice_flank=1
        a=indice_flank*self.rack.a
        phi0=a/(self.dff/2)
        p1=vm.Point2D((self._Trochoide(phi0,'T')))
        p1=p1.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        # on the coast flank 
        indice_flank=-1
        a=indice_flank*self.rack.a
        phi0=a/(self.dff/2)
        p2=vm.Point2D((self._Trochoide(phi0,'R')))
        p2=p2.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        list_2D=primitives2D.RoundedLines2D([p1,p2],{},False)
        export_2D=list_2D.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.z)
        return export_2D
    
    def _OutsideTrace(self,number):
        # Trace of the top of the gear mesh
        theta4=npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter
        p1=vm.Point2D((self.outside_diameter/2*npy.cos(theta4),self.outside_diameter/2*npy.sin(theta4)))
        p2=p1.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        p3=p2.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        list_2D=primitives2D.RoundedLines2D([p3,p2,p1],{},False)
        export_2D=list_2D.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.z)
        return export_2D
        
    def _Involute(self,tan_alpha):
        """ Involute function estimation
        
        :param tan_alpha: tan of the pressure angle
        :results: tuple of the involute point (x,y) with the origin of the involute position at the point (x=base_diameter/2, y=0) and pressure angle is positive in the direction counter clockwise
        """
        involute_list=(self.db/2*npy.cos(tan_alpha)+self.db/2*tan_alpha*npy.sin(tan_alpha),
             self.db/2*npy.sin(tan_alpha)-self.db/2*tan_alpha*npy.cos(tan_alpha))
        return involute_list
    
    def _Trochoide(self,phi,type_flank='T'):
        # function generation of trochoide point
        if type_flank=='T':
            indice_flank=1
        else:
            indice_flank=-1
        a=indice_flank*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.dff/2
        rho=self.rack.root_radius
        x2=rho*npy.sin(npy.arctan((a-r*phi)/b)-phi)+a*npy.cos(phi)-b*npy.sin(phi)+r*(npy.sin(phi)-phi*npy.cos(phi))
        y2=-rho*npy.cos(npy.arctan((a-r*phi)/b)-phi)-a*npy.sin(phi)-b*npy.cos(phi)+r*(npy.cos(phi)+phi*npy.sin(phi))
        export_point=(y2,x2)
        return export_point
    
    ### Method for ISO stress
    
    def _ISO_YS(self,s_thickness_iso):
        # stress concentration factor for ISO approach
        rho_f=self.rack.root_radius+self.rack.b**2/(self.dff/2+self.rack.b)
        coeff_ys_iso=1+0.15*s_thickness_iso/rho_f
        return coeff_ys_iso
    
    def GearISOSection(self,angle):
        """ Calculation of the ISO section
        
        :param angle: pressure angle of the ISO section calculation
        :type angle: radian
        
        :results: ISO section and ISO height
        """
        a=self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.dff/2
        theta0 = fsolve((lambda theta:a + b*npy.tan(theta) + r*(-angle - self.root_angle/2 - theta + npy.pi/2)) ,0)[0]
        phi0=(a-b*npy.tan(theta0))/r
        pt_iso = self._Trochoide(phi0)
        angle0 = npy.arctan(pt_iso[1]/pt_iso[0])-self.root_angle/2
        angle_iso = self.root_gear_angle-2*angle0
        diameter_iso = 2*norm(pt_iso)
        s_thickness_iso = diameter_iso*npy.sin(angle_iso/2)
        h_height_iso = (s_thickness_iso/2)/npy.tan(angle)
        return s_thickness_iso, h_height_iso
        
    ### Export and graph method
    
    def Dict(self):
        """Export dictionary
        
        :returns:  dictionary with all gear mesh parameters
        
        >>> print(mesh1.Dict())
        """
        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv==npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=float(v)
            else:
                d[k]=v
        
        d['rack']=self.rack.Dict()
        d['material']=self.material.Dict()
        del d['rac']
        return d
    
    def PlotData(self, x, heights, y, z, labels = True):
        transversal_plot_data = []
        axial_plot_data = []
        
#        points = []
#        for p in self.Contour(2):
#            for point in p.points:
#                points.extend(point.vector)
#        transversal_plot_data.append({'type' : 'line',
#                              'data' : points,
#                              'color' : [0, 0, 0],
#                              'dash' : 'none',
#                              'stroke_width' : 2,
#                              'marker' : '.',
#                              'size' : 1})
        # Outer diameter        
        transversal_plot_data.append({'type' : 'circle',
                          'cx' : y,
                          'cy' : z,
                          'r' : 0.5 * self.outside_diameter,
                          'color' : [0, 0, 0],
                          'size' : 1,
                          'group' : 3,
                          'dash' : 'none',})
        
        return transversal_plot_data, axial_plot_data


class MeshAssembly:
    """
    Gear mesh assembly definition
    
    :param Z: dictionary define the number of teeth {node1:Z1, node2:Z2, node3:Z3 ...}
    :param center_distance: list of the center distance in the order of connections definition [centerdistance1, centerdistance2 ...]
    :param connections: List of tuple define gear mesh connection [(node1,node2), (node2,node3), (node2,node4)]
    :param transverse_pressure_angle: dictionary define the transverse pressure angle of the first gear mesh in the connection list order
    :type transverse_pressure_angle: radian
    :param coefficient_profile_shift: dictionary define the profile shift of tooth {node1:cps1, node2:cps2, mesh3:cps3 ...}
    :param transverse_pressure_angle_rack: dictionary define the transverse pressure angle of tooth
    :param coeff_gear_addendum: dictionary define the gear addendum coefficient
    :param coeff_gear_dedendum: dictionary define the gear dedendum coefficient
    :param coeff_root_radius: dictionary define the root radius coefficient
    :param coeff_circular_tooth_thickness: dictionary define the circular tooth thickness coefficient
    :param material: Dictionary defining material for each gear mesh {node1:mechanical_components.meshes.Material, node2:mechanical_components.meshes.Material ...}
    :param torque: Dictionary defining all input torque, one node where the torque is not specified is define as the 'output' {node1:torque1, node2:torque2, node3:'output'}
    :param cycle: Dictionary defining the number of cycle for one node {node3: number_cycle3}
    :param safety_factor: Safety factor used for the ISO design
        
    >>> connections=[(1,2),(2,3)]
    >>> Z={1:13,2:41,3:37}
    >>> center_distance=[0.117,0.120]
    >>> transverse_pressure_angle=18/180.*npy.pi
    >>> mesh_assembly1=MeshAssembly(connections=connections,Z=Z,
                                    center_distance=center_distance,
                                    transverse_pressure_angle=transverse_pressure_angle)
    """
    def __init__(self,Z, center_distance, connections,transverse_pressure_angle,
                 coefficient_profile_shift,transverse_pressure_angle_rack,
                 coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
                 coeff_circular_tooth_thickness, material=None, torque=None, cycle=None,
                 safety_factor=1,verbose=False):
        
        self.center_distance=center_distance
        self.connections=connections
        self.transverse_pressure_angle_0=transverse_pressure_angle[0] # transverse pressure angle of the first gear mesh on the connections list order
        self.helix_angle=0
        
        # NetworkX graph construction
        list_gear=[]
        for gs in self.connections:
            for g in gs:
                if g not in list_gear:
                    list_gear.append(g)
        gear_graph=nx.Graph()
        gear_graph.add_nodes_from(list_gear)
        gear_graph.add_edges_from(self.connections)
        self.gear_graph=gear_graph
        self.list_gear=list_gear
        
        # Definition of default parameters
        minimum_gear_width=10e-3
        
        if material==None:
            material={list_gear[0]:hardened_alloy_steel}
        for ne in list_gear:
            if ne not in material.keys():
                material[ne]=hardened_alloy_steel
        
        if torque==None:
            torque={list_gear[0]:100,list_gear[1]:'output'}
            
        if cycle==None:
            cycle={list_gear[0]:1e6}
            
        self.material=material

        self.DF, DB, self.connections_dfs, self.transverse_pressure_angle = self.GearGeometryParameter(Z)
        if len(cycle.keys())<len(self.list_gear): # the gear mesh optimizer calculate this dictionary
            self.cycle = self.CycleParameter(cycle, Z)
        else:
            self.cycle=cycle
        dic_torque, self.normal_load, self.tangential_load, self.radial_load = self.GearTorque(Z, torque, DB)
        
        self.meshes={}
        for num_engr in self.list_gear:
            z=Z[num_engr]
            db=DB[num_engr]
            cp=coefficient_profile_shift[num_engr]
            ngp=self.list_gear.index(num_engr)
            tpa=transverse_pressure_angle_rack[num_engr]
            cga=coeff_gear_addendum[num_engr]
            cgd=coeff_gear_dedendum[num_engr]
            crr=coeff_root_radius[num_engr]
            cct=coeff_circular_tooth_thickness[num_engr]
            mat=self.material[num_engr]
            self.meshes[num_engr]=Mesh(z,db,cp,tpa,cga,cgd,crr,cct,mat)
                
        self.linear_backlash,self.radial_contact_ratio = \
            self.GearContactRatioParameter(Z, coefficient_profile_shift, DB,
                                           transverse_pressure_angle_rack,
                                           coeff_gear_addendum, coeff_gear_dedendum,
                                           coeff_root_radius,
                                           coeff_circular_tooth_thickness)
        self.gear_width,self.sigma_iso,self.sigma_lim=self.GearWidthDefinition(safety_factor,minimum_gear_width)
        for num_gear in self.list_gear:
            self.meshes[num_gear].gear_width=self.gear_width[num_gear]
            
    def Update(self, Z, center_distance, connections, transverse_pressure_angle,
               coefficient_profile_shift,
               transverse_pressure_angle_rack, coeff_gear_addendum,
               coeff_gear_dedendum, coeff_root_radius,coeff_circular_tooth_thickness,
               material, torque,cycle, safety_factor,verbose):
        """ Update of the gear mesh assembly
        
        :param all: same parameters of this class initialisation
            
        >>> Z={1:13,2:46,4:38}
        >>> center_distance=[0.118,0.125]
        >>> mesh_assembly1.Update(Z=Z,center_distance=center_distance)
        """
        self.center_distance=center_distance
        self.transverse_pressure_angle_0=transverse_pressure_angle[0]
        self.DF,dict_db,self.connections_dfs,self.transverse_pressure_angle=self.GearGeometryParameter(Z)
        self.linear_backlash,self.radial_contact_ratio = self.GearContactRatioParameter(
            Z, coefficient_profile_shift, dict_db, transverse_pressure_angle_rack,
            coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
            coeff_circular_tooth_thickness)
        
    ### Method gear mesh calculation
        
    def GearContactRatioParameter(self,Z,coefficient_profile_shift,DB,
                           transverse_pressure_angle_rack,coeff_gear_addendum,
                           coeff_gear_dedendum,coeff_root_radius,
                           coeff_circular_tooth_thickness):
        
        for num_engr in self.list_gear:
            z=Z[num_engr]
            db=DB[num_engr]
            cp=coefficient_profile_shift[num_engr]
            tpa=transverse_pressure_angle_rack[num_engr]
            cga=coeff_gear_addendum[num_engr]
            cgd=coeff_gear_dedendum[num_engr]
            crr=coeff_root_radius[num_engr]
            cct=coeff_circular_tooth_thickness[num_engr]
            mat=self.material[num_engr]
            self.meshes[num_engr].Update(z,db,cp,tpa,cga,cgd,crr,cct,mat)
            
        linear_backlash=[]
        radial_contact_ratio=[]
        for engr1,engr2 in self.connections_dfs:
            if (engr1,engr2) in self.connections:
                num_mesh=self.connections.index((engr1,engr2))
            elif (engr2,engr1) in self.connections:
                num_mesh=self.connections.index((engr2,engr1))
            else:
                raise RuntimeError
            circular_tooth_thickness1=self.meshes[engr1].GearSection(self.DF[num_mesh][engr1])
            circular_tooth_thickness2=self.meshes[engr2].GearSection(self.DF[num_mesh][engr2])
            transverse_radial_pitch1=npy.pi*self.DF[num_mesh][engr1]/self.meshes[engr1].z
            space_width1=transverse_radial_pitch1-circular_tooth_thickness1
            space_width2=transverse_radial_pitch1-circular_tooth_thickness2    
            linear_backlash.append(min(space_width1-circular_tooth_thickness2,space_width2-circular_tooth_thickness1))
            transverse_pressure_angle1=self.transverse_pressure_angle[num_mesh]
            center_distance1=self.center_distance[num_mesh]
            radial_contact_ratio.append((1/2*(npy.sqrt(self.meshes[engr1].outside_diameter**2
                                                       - self.meshes[engr1].db**2)
                                              + npy.sqrt(self.meshes[engr2].outside_diameter**2
                                                         - self.meshes[engr2].db**2)
                                        - 2*center_distance1*npy.sin(transverse_pressure_angle1))
                                        /(transverse_radial_pitch1*npy.cos(transverse_pressure_angle1))))
        return linear_backlash,radial_contact_ratio

    def GearGeometryParameter(self,Z):
        # Construction of pitch and base diameter
        DF={}
        db={}
        dict_transverse_pressure_angle={0:self.transverse_pressure_angle_0}
        connections_dfs=list(nx.edge_dfs(self.gear_graph, 
                            [self.connections[0][0],self.connections[0][1]]))
        for num_dfs,((engr1,engr2),cd) in enumerate(zip(connections_dfs,self.center_distance)):
            if (engr1,engr2) in self.connections:
                num_mesh=self.connections.index((engr1,engr2))
            else:
                num_mesh=self.connections.index((engr2,engr1))
            Z1=Z[engr1]
            Z2=Z[engr2]
            DF1=2*cd*Z1/Z2/(1+Z1/Z2)
            DF2=2*cd-DF1
            DF[num_mesh]={}
            DF[num_mesh][engr1]=DF1
            DF[num_mesh][engr2]=DF2
            if num_mesh==0:
                db1=float(DF1*npy.cos(self.transverse_pressure_angle_0))
                db2=float(DF2*npy.cos(self.transverse_pressure_angle_0))
            else:
                db1=db[engr1]
                dict_transverse_pressure_angle[num_mesh]=npy.arccos(db1/DF1)
                db2=DF2*npy.cos(dict_transverse_pressure_angle[num_mesh])
            db[engr1]=db1
            db[engr2]=db2
        transverse_pressure_angle=[]
        for num_mesh in sorted(dict_transverse_pressure_angle.keys()):
            tpa=dict_transverse_pressure_angle[num_mesh]
            transverse_pressure_angle.append(tpa)
        
        return DF,db,connections_dfs,transverse_pressure_angle
    
    def GearTorque(self,Z,torque,db):
        """ Calculation of the gear mesh torque
        
        :param Z: dictionary define the number of teeth {node1:Z1, node2:Z2, mesh3:Z3 ...}
        :param torque: dictionary defining all input torque, one node where the torque is not specified is define as the 'output' {node1:torque1, node2:torque2, node3:'output'}
        :param db: dictionary define the base diameter {mesh1: {node1:db1_a, node2:db2_a}, mesh2: {node2:db2_b, node3:db3_b}}
        :type db: m
        
        :results:
            * **torque1** - dictionary of the applied torque on gear mesh (torque applied by node_x on node_y) {node1:tq1, node3:tq3 ...}
            * **torque2** - dictionary of the counter drive torque on gear mesh (torque applied by node_x on node_y) {node2:-tq1, node4:-tq3 ...}
            * **normal_load** - dictionary define the normal load for each gear mesh (applied torque define the direction) {mesh1 : [fn_x1,fn_y1,fn_z1],mesh2 : [fn_x2,fn_y2,fn_z2] ...}
            * **tangential_load** - dictionary define the tangential load for each gear mesh (applied torque define the direction) {mesh1 : [ft_x1,ft_y1,ft_z1],mesh2 : [ft_x2,ft_y2,ft_z2] ...}
            * **radial_load** - dictionary define the radial load for each gear mesh (applied torque define the direction) {mesh1 : [fr_x1,fr_y1,fr_z1],mesh2 : [fr_x2,fr_y2,fr_z2] ...}
        
        be careful, due to the parameters of the gear mesh assembly (define one pressure angle for each mesh) the diameter db2_a is different to db2_b (you have to define correctly transverse_pressure_angle to have db2_a=db2_b)
        """
        if 'output' in torque.values():
            for num_gear,tq in torque.items():
                if tq=='output':
                    node_output=num_gear
            torque_graph_dfs=list(nx.dfs_edges(self.gear_graph,node_output))
            order_torque_calculation=[(eng2,eng1) for (eng1,eng2) in torque_graph_dfs[::-1]]
            # calculation torque distribution
            temp_torque={}
            for eng1 in self.list_gear:
                temp_torque[eng1]=0
            for num_mesh_tq,(eng1,eng2) in enumerate(order_torque_calculation):
                if eng1 in torque.keys():
                    temp_torque[eng1]+=torque[eng1]
                temp_torque[eng2]+=-temp_torque[eng1]*Z[eng2]/Z[eng1]
            dic_torque={}
            for num_mesh_tq,(eng1,eng2) in enumerate(order_torque_calculation):
                dic_torque[(eng1,eng2)]=temp_torque[eng1]
    
#        ggd=self.gear_graph.degree(self.list_gear) # dictionary: key num of the gear and the value the number of connection
#        liste_node_init=[]
#        for num_gear,nb_connexion in ggd:
#            if (nb_connexion==1) and (torque[num_gear]!='output'):
#                liste_node_init.append(num_gear)
#        for num_gear,tq in torque.items():
#            if tq=='output':
#                node_output=num_gear
#        if 'output' not in torque.values():
#            path_list=[nx.shortest_path(self.gear_graph,source=liste_node_init[0],target=liste_node_init[1])]
#        else:
#            path_list=[nx.shortest_path(self.gear_graph,source=liste_node_init[0],target=node_output)]
#            for nd_init in liste_node_init[1:]:
#                path_list.append(nx.shortest_path(self.gear_graph,source=nd_init,target=node_output))
#        torque1={}
#        normal_load={}
#        tangential_load={}
#        radial_load={}
#        for node_init,path in zip(liste_node_init,path_list):
#            torque_moteur_m=torque[node_init]
#            for i,eng1 in enumerate(path[0:-1]):
#                eng2=path[i+1]
#                torque_recepteur=-torque_moteur_m*Z[eng2]/Z[eng1]
#                if (eng2 in torque.keys()) and (torque[eng2]!='output'):
#                    torque_recepteur+=torque[eng2]
#                torque1[eng1]=torque_moteur_m
#                torque2={}
#                torque2[eng2]=torque_recepteur
                 
        normal_load={}
        tangential_load={}
        radial_load={}
        for num_mesh,(eng1,eng2) in enumerate(self.connections):
            if 'output' not in torque.values():
                dic_torque=torque
            try:
                tq=dic_torque[(eng1,eng2)]
            except:
                tq=dic_torque[(eng2,eng1)]
            normal_load[num_mesh]=abs(tq)*2/(db[eng1])
            tangential_load[num_mesh]=abs(tq)*2/(self.DF[num_mesh][eng1])
            radial_load[num_mesh]=npy.tan(self.transverse_pressure_angle[num_mesh])*tangential_load[num_mesh]
        return dic_torque,normal_load,tangential_load,radial_load
    
    def CycleParameter(self,cycle,Z):
        """ Calculation of the gear mesh cycle
        
        :param Z: dictionary define the number of teeth {node1:Z1, node2:Z2, node3:Z3 ...}
        :param cycle: Dictionary defining the number of cycle for one node {node3: number_cycle3}
        
        :results: dictionary define the number of cycle for each gear mesh {node1:cycle1, node2:cycle2, node3:cycle3 ...}
        """
        eng_init=list(cycle.keys())[0]
        for eng in self.list_gear:
            if eng not in cycle.keys():
                cycle[eng]=cycle[eng_init]*Z[eng_init]/Z[eng]
        return cycle
        
    def GearWidthDefinition(self,safety_factor,minimum_gear_width):
        """ Calculation of the gear width
        
        :param safety_factor: Safety factor used for the ISO design
        
        :results:
            * **gear_width** - dictionary define the gear mesh width {node1 : gw1, node2 : gw2, node3 : gw3 ...}
            * **sigma_iso** - dictionary define the ISO stress {mesh1 : {node1 sig_iso1: , node2 : sig_iso2_1}, mesh2 : {node2 : sig_iso2_2, node3 : sig_iso3} ...}
            * **sigma_lim** - dictionary define the limit material stress {mesh1 : {node1 sig_lim1: , node2 : sig_lim2}, mesh2 : {node2 : sig_lim2, node3 : sig_lim3} ...}
            
        in this function, we define the gear width for each gear mesh to respect sig_lim = sig_iso for each gear mesh
        """
        coeff_yf_iso=self._CoeffYFIso()
        coeff_ye_iso=self._CoeffYEIso()
        coeff_yb_iso=self._CoeffYBIso()
        
        sigma_lim=self.SigmaMaterialISO(safety_factor)
        gear_width={}
        for eng in self.list_gear:
            gear_width[eng]=minimum_gear_width
        for num_mesh,(eng1,eng2) in enumerate(self.connections):
            gear_width1=abs(self.tangential_load[num_mesh]
                            / (sigma_lim[num_mesh][eng1]
                                * self.meshes[eng1].rack.module)
                            *coeff_yf_iso[num_mesh][eng1]
                            *coeff_ye_iso[num_mesh]
                            *coeff_yb_iso[num_mesh][eng1])
                            
            gear_width2=abs(self.tangential_load[num_mesh]
                            /(sigma_lim[num_mesh][eng2]
                                *self.meshes[eng2].rack.module)
                            *coeff_yf_iso[num_mesh][eng2]
                            *coeff_ye_iso[num_mesh]
                            *coeff_yb_iso[num_mesh][eng2])
                            
            gear_width_set=max(gear_width1,gear_width2)
            gear_width[eng1]=max(gear_width[eng1],gear_width_set)
            gear_width[eng2]=max(gear_width[eng2],gear_width_set)
        sigma_iso=sigma_lim
        return gear_width,sigma_iso,sigma_lim
        
    def SigmaMaterialISO(self,safety_factor):
        """ Calculation of the material limit stress
        
        :param safety_factor: Safety factor used for the ISO design
        
        :results:
            * **sigma_lim** - dictionary define the limit material stress {mesh1 : {node1 sig_lim1: , node2 : sig_lim2}, mesh2 : {node2 : sig_lim2, node3 : sig_lim3} ...}
            
        in this function, we use the FunCoeff function of the Material class to interpolate the material parameters
        """
        angle=30/180*npy.pi
        sigma_lim={}
        for num_mesh,(eng1,eng2) in enumerate(self.connections):
            sigma_lim[num_mesh]={}
            
            matrice_wholer=self.material[eng1].data_wholer_curve
            matrice_material=self.material[eng1].data_gear_material
            sgla=self.material[eng1].FunCoeff(self.cycle[eng1],npy.array(matrice_wholer['data']),matrice_wholer['x'],matrice_wholer['y'])
            sgl1=self.material[eng1].FunCoeff(sgla,npy.array(matrice_material['data']),matrice_material['x'],matrice_material['y'])
            s_thickness_iso_1,h_height_iso_1=self.meshes[eng1].GearISOSection(angle)
            coeff_ys_iso=self.meshes[eng1]._ISO_YS(s_thickness_iso_1)
            sigma_lim[num_mesh][eng1]=(sgl1/(safety_factor*coeff_ys_iso))*10**7
            
            matrice_wholer=self.material[eng2].data_wholer_curve
            matrice_material=self.material[eng2].data_gear_material
            sglb=self.material[eng2].FunCoeff(self.cycle[eng2],npy.array(matrice_wholer['data']),matrice_wholer['x'],matrice_wholer['y'])
            sgl2=self.material[eng2].FunCoeff(sglb,npy.array(matrice_material['data']),matrice_material['x'],matrice_material['y'])
            s_thickness_iso_2,h_height_iso_2=self.meshes[eng2].GearISOSection(angle)
            coeff_ys_iso=self.meshes[eng2]._ISO_YS(s_thickness_iso_2)
            sigma_lim[num_mesh][eng2]=(sgl2/(safety_factor*coeff_ys_iso))*10**7
        return sigma_lim

        
    def _CoeffYFIso(self):
        # shape factor for ISO stress calculation
        angle=30/180*npy.pi
        coeff_yf_iso={}
        for num_mesh,(eng1,eng2) in enumerate(self.connections):
            coeff_yf_iso[num_mesh]={}
            s_thickness_iso_1,h_height_iso_1=self.meshes[eng1].GearISOSection(angle)
            s_thickness_iso_2,h_height_iso_2=self.meshes[eng2].GearISOSection(angle)
            coeff_yf_iso[num_mesh][eng1]=((6*(h_height_iso_1/self.meshes[eng1].rack.module)*npy.cos(self.transverse_pressure_angle[num_mesh]))
                                    /((s_thickness_iso_1/self.meshes[eng1].rack.module)**2
                                       *npy.cos(self.meshes[eng1].rack.transverse_pressure_angle)))
            coeff_yf_iso[num_mesh][eng2]=((6*(h_height_iso_2/self.meshes[eng2].rack.module)*npy.cos(self.transverse_pressure_angle[num_mesh]))
                                    /((s_thickness_iso_2/self.meshes[eng2].rack.module)**2
                                       *npy.cos(self.meshes[eng2].rack.transverse_pressure_angle)))
        return coeff_yf_iso
        
    def _CoeffYEIso(self):
        #  radial contact ratio factor for ISO stress calculation
        coeff_ye_iso=[]
        for ne,eng in enumerate(self.connections):
            coeff_ye_iso.append(1/self.radial_contact_ratio[ne])
        return coeff_ye_iso
        
    def _CoeffYBIso(self):
        # gear widht factor impact for ISO stress calculation
        coeff_yb_iso={}
        for num_mesh,(eng1,eng2) in enumerate(self.connections):
            coeff_yb_iso[num_mesh]={}
            matrice_YB=self.material[eng1].data_coeff_YB_Iso
            coeff_yb_iso[num_mesh][eng1]=self.material[eng1].FunCoeff(self.helix_angle,npy.array(matrice_YB['data']),matrice_YB['x'],matrice_YB['y'])
            matrice_YB=self.material[eng2].data_coeff_YB_Iso
            coeff_yb_iso[num_mesh][eng2]=self.material[eng2].FunCoeff(self.helix_angle,npy.array(matrice_YB['data']),matrice_YB['x'],matrice_YB['y'])
        return coeff_yb_iso
    
    ### Optimization Method
    
    def CheckMinimumBacklash(self,backlash_min=2*1e-4):
        """ Define constraint and functional for the optimizer on backlash
        
        :param backlash_min: maximum backlash available
        :results: 
            * check is a boolean (True if 0<backlash<backlash_min)
            * list_ineq a list of element that should be positive for the optimizer
            * obj is a functional on the backlash used for the optimizer
        """
        list_ineq=[] # liste of value to evaluate backlash
        obj=0
        for lb in self.linear_backlash:
            list_ineq.append(lb) # backlash > 0
            list_ineq.append(backlash_min-lb) # backlash < backlash_min so (backlash_min-backlash)>0
            obj+=10*(lb-backlash_min)**2
        check=False
        if min(list_ineq)>0:
            check=True
        return check,list_ineq,obj
    
    def CheckRadialContactRatio(self,radial_contact_ratio_min=1):
        """ Define constraint and functional for the optimizer on radial contact ratio
        
        :param radial_contact_ratio_min: minimum radial contact ratio available
        :results: 
            * check is a boolean (True if radial_contact_ratio_min<radial_contact_ratio)
            * list_ineq a list of element that should be positive for the optimizer
            * obj is a functional on the backlash used for the optimizer
        """
        list_ineq=[]
        obj=0
        for num_mesh,(eng1,eng2) in enumerate(self.connections):
            rca=self.radial_contact_ratio[num_mesh]
            list_ineq.append(rca-radial_contact_ratio_min)
            if rca>radial_contact_ratio_min:
                obj+=0.001*(rca-radial_contact_ratio_min)
            else:
                obj+=1000*(radial_contact_ratio_min-rca)
        check=False
        if min(list_ineq)>0:
            check=True
        return check,list_ineq,obj
    
    def ListeIneq(self):
        """ Compilation method for inequality list used by the optimizer
        
        :results: vector of data that should be positive
        """
        _,ineq,_=self.CheckMinimumBacklash(4*1e-4)
        _,list_ineq,_=self.CheckRadialContactRatio(1)
        ineq.extend(list_ineq)
        
        for num_gear,mesh in self.meshes.items():
            list_ineq=mesh.ListeIneq()
            ineq.extend(list_ineq)
        
        return ineq
        
    def Functional(self):
        """ Compilation method for a part of the functional used by the optimizer
        
        :results: scalar add to the global functional of the optimizer
        """
        check1,ineq1,obj1=self.CheckMinimumBacklash(4*1e-4)
        check2,ineq2,obj2=self.CheckRadialContactRatio(1)
        obj=obj1+obj2
        return obj
        
    ### Function graph and export
    
    def GearRotate(self,list_gear,list_center,list_rot):
        """ Displacement of the volmdlr gear profile (rotation and translation)
        
        :param list_gear: list of volmdlr contour [meshes.Contour, meshes.Contour ...], each contour is centered on the origin
        :param list_center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param list_rot: list of rotation for each gear mesh [node1 : rot1, node2 : rot2 ...]
        
        :results: list of volmdlr component
        """
        export=[]
        for (i,center,k) in zip(list_gear,list_center,list_rot):
            model_export=[]
            for m in i:
                model_trans=m.Translation(center)
                model_trans_rot=model_trans.Rotation(vm.Point2D(center),k)
                model_export.append(model_trans_rot)
            export.append(model_export)
        return export
    
    def InitialPosition(self,set_pos,liste_eng=()):    
        """ Calculation of the rotation for two gear mesh to initiate the contact
        
        :param list_gear: list of volmdlr contour [meshes.Contour, meshes.Contour ...], each contour is centered on the origin
        :param list_center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param list_rot: list of rotation for each gear mesh [node1 : rot1, node2 : rot2 ...]
        
        :results: list of volmdlr component
        """
        Angle1=npy.arccos(self.meshes[liste_eng[0]].db/self.DF[set_pos][liste_eng[0]])
        Angle2=npy.arccos(self.meshes[liste_eng[1]].db/self.DF[set_pos][liste_eng[1]])
        Gear1Angle=-(npy.tan(Angle1)-Angle1)
        Gear2Angle=-(npy.tan(Angle2)-Angle2)+npy.pi        
        return [Gear1Angle,Gear2Angle]
    
    def VolumeModel(self, centers = {}, axis = (1,0,0), name = ''):
        """ Generation of the 3D volume for all the gear mesh
        
        :param center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param axis: direction of gear mesh rotation
        
        :results: list of 3D volmdlr component
        """
        x = vm.Vector3D(axis)
        y = x.RandomUnitNormalVector()
        z = vm.Vector3D(npy.cross(x.vector, y.vector))  
        if len(centers)==0:
            centers = {}
            center_var = self.PosAxis({self.list_gear[0]:[0,0]})
            for engr_num in center_var.keys():
                centers[engr_num]=tuple(center_var[engr_num][0]*y.vector+center_var[engr_num][1]*z.vector)
        else:
            center_var={}
            for engr_num in centers.keys():
                center_var[engr_num]=(npy.dot(centers[engr_num],x.vector),npy.dot(centers[engr_num],y.vector),npy.dot(centers[engr_num],z.vector))
            centers=center_var
            
        Gears3D={}
        Struct=[]
        Rotation={}
        primitives=[]
        
        for set_pos_dfs,(eng1,eng2) in enumerate(self.connections_dfs):
            position1 = centers[eng1]
            position2 = centers[eng2]
            
            if (eng1,eng2) in self.connections:
                set_pos=self.connections.index((eng1,eng2))
                list_rot=self.InitialPosition(set_pos,(eng1,eng2))
            elif (eng2,eng1) in self.connections:
                set_pos=self.connections.index((eng2,eng1))
                list_rot=self.InitialPosition(set_pos,(eng2,eng1))
            Rotation[set_pos]={}
            if set_pos_dfs==0:
                Gears3D[eng1]=self.meshes[eng1].Contour(3)
            Struct.append(vm.Circle2D(vm.Point2D(position1),self.DF[set_pos][eng1]/2))
            Gears3D[eng2]=self.meshes[eng2].Contour(3)
            Struct.append(vm.Circle2D(vm.Point2D(position2),self.DF[set_pos][eng2]/2))            
            
            if position2[1]==position1[1]:
                if position2[2]-position1[2]>0:
                    angle0=npy.pi/2
                else:
                    angle0=-npy.pi/2
            else:
                angle0=-npy.arctan((position2[2]-position1[2])/(position2[1]-position1[1]))
                if (position2[2]-position1[2])<0:
                    angle0=angle0+npy.pi
            if set_pos_dfs==0:
                Rotation[set_pos][eng1]=list_rot[0]+angle0
                Rotation[set_pos][eng2]=list_rot[1]+angle0
            else:
                for k1,rot in Rotation.items():
                    if eng1 in rot.keys():
                        Rotation[set_pos][eng1]=rot[eng1]
                        delta_rot=Rotation[set_pos][eng1]-(list_rot[0]-angle0)
                Rotation[set_pos][eng2]=list_rot[1]-angle0-delta_rot*((self.meshes[eng1].z)/(self.meshes[eng2].z))
            Gears3D_Rotate=self.GearRotate([Gears3D[eng1],Gears3D[eng2]],[(position1[1::]),(position2[1::])],
                                       list_rot=[Rotation[set_pos][eng1],Rotation[set_pos][eng2]])
        
            C1=vm.Contour2D(Gears3D_Rotate[0])
            C2=vm.Contour2D(Gears3D_Rotate[1])
            
            extrusion_vector1 = (self.gear_width[eng1]*x).vector
            extrusion_vector2 = (self.gear_width[eng2]*x).vector
            
            if set_pos_dfs==0:
                vect_x=tuple(-0.5*self.gear_width[eng1]*x.vector+[npy.dot(centers[eng1],x.vector),0,0])
                t1=primitives3D.ExtrudedProfile(vm.Vector3D(vect_x),y,z,[C1],extrusion_vector1)
                primitives.append(t1)
        
            vect_x=tuple(-0.5*self.gear_width[eng2]*x.vector+[npy.dot(centers[eng2],x.vector),0,0])
            t2=primitives3D.ExtrudedProfile(vm.Vector3D(vect_x),y,z,[C2],extrusion_vector2)
            primitives.append(t2)

        model=vm.VolumeModel(primitives)
        return model
    
    def Mass(self):
        """ Estimation of gear mesh mass
        
        :results: mass of all gear mesh
        """
        DF = [0]*len(self.gear_width.keys())
        for i,(ic1, ic2) in enumerate(self.connections):
            DF[ic1] = self.DF[i][ic1]
            DF[ic2] = self.DF[i][ic2]
        
        mass = 0.
        for i,df in enumerate(DF):
            mass +=  self.gear_width[i] * self.material[i].volumic_mass* math.pi * (0.5*DF[i])**2
        return mass
    
    # Waiting for meshes to know how to plot themselves
#    def PlotData(self, x, heights, ys, zs, labels = True):
#        transversal_plot_data = []
#        axial_plot_data = []
#        # TODO remove when meshes would be a list
#        imesh = []
#        meshes = []
#        for ic, connec in enumerate(self.connections):
#            imesh.extend(connec)
#        imesh = list(set(imesh))
#        for ic, (ic1, ic2) in enumerate(self.connections):
#            if ic1 in imesh:
#                meshes.append(self.meshes[ic][ic1])
#                imesh.remove(ic1)
#            if ic2 in imesh:
#                meshes.append(self.meshes[ic][ic2])
#                imesh.remove(ic2)
#        
#        for imesh, mesh in enumerate(meshes):
#            t, a = mesh.PlotData(x, heights, ys[imesh], zs[imesh], labels)
#            transversal_plot_data.extend(t)
#            axial_plot_data.extend(a)
#            
#        # Ploting axial because mesh doesn't know its width
#        
#        return axial_plot_data, transversal_plot_data

    def FreeCADExport(self, file_path, centers = {}, axis = (1,0,0), export_types=['fcstd'], python_path = 'python',
                      freecad_path = '/usr/lib/freecad/lib'):
        """ Export 3D volume to FreeCAD
        
        :param file_path: file path for the freecad file
        :param center: list of tuple define the final position of the gear mesh center (a translation is perform, then a rotation around this axis)
        :param axis: direction of gear mesh rotation
        
        :results: export of a FreeCAD file
        """
        model = self.VolumeModel(centers, axis)
        model.FreeCADExport(python_path ,file_path, freecad_path, export_types)
        
    def PosAxis(self,position):
        # Definition of the initial center for all gear (when not given by the user)
        def fun(x):
            obj=0
            for num,it in enumerate(self.connections):
                eng1=(self.list_gear).index(it[0])
                eng2=(self.list_gear).index(it[1])
                obj+=(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-self.center_distance[num])**2
            return obj
        def eg(x):
            ine=[]
            for k,val in position.items():
                key=(self.list_gear).index(k)
                ine.append(x[2*int(key)]-val[0])
                ine.append(x[2*int(key)+1]-val[1])
            return ine
        def ineg(x):
            ine=[]
            for num,it in enumerate(self.connections):
                eng1=(self.list_gear).index(it[0])
                eng2=(self.list_gear).index(it[1])
                ine.append(((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5-0.999*self.center_distance[num])
                ine.append(1.001*self.center_distance[num]-((x[2*eng1]-x[2*eng2])**2+(x[2*eng1+1]-x[2*eng2+1])**2)**0.5)
            return ine
        cons = ({'type': 'eq','fun' : eg},{'type': 'ineq','fun' : ineg})
        drap=1
        while drap==1:
            x0=tuple(npy.random.random(2*self.gear_graph.number_of_nodes())*1)
            Bound=[[0,1]]*(self.gear_graph.number_of_nodes()*2)
            res = minimize(fun,x0, method='SLSQP', bounds=Bound,constraints=cons)
            if (min(ineg(res.x))>0) and (max(eg(res.x))<1e-7):
                drap=0
        x_opt=res.x
        centers={}
        for engr_pos,engr_num in enumerate(self.list_gear):
            centers[engr_num]=[x_opt[2*engr_pos],x_opt[2*engr_pos+1]]
        return centers
    
    # TODO change this function to make it like PlotData in PWT: output is dict of geometrical shapes
    def SVGExport(self,name,position):
        """ Export SVG graph of all gear mesh
        
        :param name: name of the svg file
        :param position: dictionary define some center position {node2 : [0,0], node4 : [0.12,0] ..}
        
        :results: SVG graph
        
        in the position dictionary, you have to be coherent with the center position
        
            * for exemple, if the center-distance of the mesh1 (node1, node2) is 0.117 m you can define position such as:
        
                * {node1 : [0,0], node2 : [0.117,0]}
                * {node1 : [0,0]}
        """
        x_opt=self.PosAxis(position)
        TG={}
        L1=[]
        Struct=[]
        Rot={}
        for num,en in enumerate(self.connections_dfs):
            position1=(x_opt[en[0]][0],x_opt[en[0]][1])   
            position2=(x_opt[en[1]][0],x_opt[en[1]][1])
            #tuple1 et 2 correspondent a la position des centres
            ne=self.connections.index(en)
            Rot[ne]={}
            if num==0:
                TG[en[0]]=self.meshes[en[0]].Contour(5)
            Struct.append(vm.Circle2D(vm.Point2D(position1),self.DF[ne][en[0]]/2))
            TG[en[1]]=self.meshes[en[1]].Contour(5)
            Struct.append(vm.Circle2D(vm.Point2D(position2),self.DF[ne][en[1]]/2))
            #Definition de la position angulaire initiale
            list_rot=self.InitialPosition(ne,en)
            if position2[0]==position1[0]:
                if position2[1]-position1[1]>0:
                    angle=npy.pi/2
                else:
                    angle=-npy.pi/2
            else:
                angle=-npy.arctan((position2[1]-position1[1])/(position2[0]-position1[0]))
            if num==0:
                Rot[ne][en[0]]=list_rot[0]-angle
                Rot[ne][en[1]]=list_rot[1]-angle
            else:
                for k1,v1 in Rot.items():
                    if en[0] in v1.keys():
                        Rot[ne][en[0]]=v1[en[0]]
                        delta_rot=Rot[ne][en[0]]-(list_rot[0]-angle)
                Rot[ne][en[1]]=list_rot[1]-angle-delta_rot*((self.meshes[en[0]].z)/(self.meshes[en[1]].z))
            sol=self.GearRotate([TG[en[0]],TG[en[1]]],[position1,position2],list_rot=[Rot[ne][en[0]],Rot[ne][en[1]]])
            if num==0:
                L1.extend(sol[0])
            L1.extend(sol[1])
        L1.extend(Struct)
        G1=vm.Contour2D(L1)
        G1.MPLPlot()
        
    def Dict(self):
        """Export dictionary
        
        :returns:  dictionary with all gear mesh assembly parameters
        """
        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv==npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=float(v)
            else:
                d[k]=v
        d['meshes']={}
        for k,v in self.meshes.items():
            d['meshes'][k]=v.Dict()
        del d['gear_graph']
        del d['material']
        del d['list_gear']
        return d

