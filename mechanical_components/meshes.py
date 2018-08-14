import numpy as npy
#import math as mt
from scipy import interpolate
#import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm
from scipy.optimize import minimize,fsolve
#from scipy.interpolate import splprep, splev
import itertools
#from jinja2 import Environment, PackageLoader, select_autoescape
import networkx as nx

#import mechanical_components.LibSvgD3 as LibSvg
import powertransmission.tools as tools

#import pyDOE

# TODO: Vérifier que les commentaires et les noms de variables soient en anglais


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
    def __init__(self,volumic_mass, data_coeff_YB_Iso, data_wholer_curve,
                 data_gear_material):
        self.volumic_mass = volumic_mass
        self.data_coeff_YB_Iso = data_coeff_YB_Iso
        self.data_wholer_curve = data_wholer_curve
        self.data_gear_material = data_gear_material
        
    def FunCoeff(self,x,data,type_x='Linear',type_y='Linear'):
        if type_x=='Log': 
            x=npy.log10(x)
        f = interpolate.interp1d(list(data[:,0]),list(data[:,1]), fill_value='extrapolate')
        sol=float(f(x))
        if type_y=='Log':
            sol=10**sol
        return sol
    
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

class Rack():
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

        #paramètre pour la trochoide
        self.a=self.tooth_space/2-self.gear_dedendum*npy.tan(self.transverse_pressure_angle)-self.root_radius*npy.tan(1/2*npy.arctan(npy.cos(self.transverse_pressure_angle)/(npy.sin(self.transverse_pressure_angle))))
        self.b=self.gear_dedendum-self.root_radius

    def Update(self,module,transverse_pressure_angle,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness):
        self.module=module
        self.RackParam(transverse_pressure_angle,coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,coeff_circular_tooth_thickness)
        

    def Dict(self):
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

    def CSVExport(self):
        d=self.__dict__.copy()
        return list(d.keys()),list(d.values())

class Mesh():
    def __init__(self, z, db, cp, transverse_pressure_angle_rack,
                 coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
                 coeff_circular_tooth_thickness):
        
        self.rack=Rack(transverse_pressure_angle_rack)
        self.GearParam(z,db,cp,transverse_pressure_angle_rack,
                       coeff_gear_addendum, coeff_gear_dedendum, 
                       coeff_root_radius, coeff_circular_tooth_thickness)
        
    def Update(self,z,db,cp,transverse_pressure_angle_rack,
               coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,
               coeff_circular_tooth_thickness):
        self.GearParam(z,db,cp, transverse_pressure_angle_rack,
                       coeff_gear_addendum, coeff_gear_dedendum,
                       coeff_root_radius, coeff_circular_tooth_thickness)
        
    def GearParam(self,z, db, cp, transverse_pressure_angle_rack,
                  coeff_gear_addendum, coeff_gear_dedendum,
                  coeff_root_radius, coeff_circular_tooth_thickness):
        
        self.Z=z
        self.DB=db
        self.DFF=self.DB/npy.cos(transverse_pressure_angle_rack)
        module_rack=self.DFF/self.Z
        self.rack.Update(module_rack,transverse_pressure_angle_rack,
                         coeff_gear_addendum,coeff_gear_dedendum, 
                         coeff_root_radius,coeff_circular_tooth_thickness)
        self.coefficient_profile_shift=cp
        
        self.outside_diameter=(self.DFF
                               +2*(self.rack.gear_addendum
                                   +self.rack.module*self.coefficient_profile_shift))
        self.alpha_outside_diameter = npy.arccos(self.DB/self.outside_diameter)
        self.root_diameter=(self.DFF
                            - 2*(self.rack.gear_dedendum
                                - self.rack.module*self.coefficient_profile_shift))
        self.root_diameter_active,self.phi_trochoide=self._RootDiameterActive()
        self.alpha_root_diameter_active=npy.arccos(self.DB/self.root_diameter_active)
        
        self.alpha_pitch_diameter=npy.arccos(self.DB/self.DFF)
        self.circular_tooth_thickness = (self.rack.circular_tooth_thickness
                                       +self.rack.module*self.coefficient_profile_shift
                                           *npy.tan(self.rack.transverse_pressure_angle)
                                       +self.rack.module*self.coefficient_profile_shift
                                           *npy.tan(self.rack.transverse_pressure_angle))
        self.tooth_space=self.rack.transverse_radial_pitch-self.circular_tooth_thickness
        self.outside_active_angle = (2*self.circular_tooth_thickness/self.DFF-2
                                     *(npy.tan(self.alpha_outside_diameter)
                                        -self.alpha_outside_diameter
                                        -npy.tan(self.alpha_pitch_diameter)
                                        +self.alpha_pitch_diameter))
        self.base_circular_tooth_thickness = (self.DB/2
                                              *(2*self.circular_tooth_thickness/self.DFF
                                                +2*(npy.tan(self.alpha_pitch_diameter)
                                                -self.alpha_pitch_diameter)))
        
        self.root_angle=self.tooth_space/(self.DFF/2)-2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
        self.root_gear_angle=self.circular_tooth_thickness/(self.DFF/2)+2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
        
    def GearSection(self,diameter):
        # TODO: traduire en englais le prochain commentaire
        #epaisseur de la dent au diameter
        alpha_diameter=npy.arccos(self.DB/diameter)
        theta1=(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter)-(npy.tan(alpha_diameter)-alpha_diameter)
        return diameter/2*(2*theta1+self.outside_active_angle)
        
    def _RootDiameterActive(self):
        # TODO: traduire en englais le prochain commentaire
        #Analyse diam pied de dent actif
        a=self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.DFF/2
#        rho=self.rack.root_radius_T
        phi=-(a+b*npy.tan(npy.pi/2-self.rack.transverse_pressure_angle))/r
        root_diameter_active=2*norm(self._Trochoide(phi))
        return root_diameter_active,phi
    
    def Contour(self,discret=10,list_number=[None]):
        #Analytical tooth profil
        if list_number==[None]:
            list_number=npy.arange(int(self.Z))
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
        
#        theta=npy.linspace(0,npy.tan(self.alpha_outside_diameter),discret)
        sol=self._Involute(drap*theta)
        x=sol[0]
        y=sol[1]
        p=[vm.Point2D((x[0],y[0]))]
        for i in range(1,discret):
            p.append(vm.Point2D((x[i],y[i])))
        ref=primitives2D.RoundedLines2D(p,{},False)
        if ind=='T':
            L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
            self.rac=L.points[-1]
        else:
            L=ref.Rotation(vm.Point2D((0,0)),
                           self.base_circular_tooth_thickness*2/self.DB)
            L=L.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
            L.points[0]=self.rac
        return L
        
    def _Involute(self,tan_alpha):
        # TODO: better variable naming than "sol"
        sol=(self.DB/2*npy.cos(tan_alpha)+self.DB/2*tan_alpha*npy.sin(tan_alpha),
             self.DB/2*npy.sin(tan_alpha)-self.DB/2*tan_alpha*npy.cos(tan_alpha))
        return sol
    
    def _TrochoideTrace(self, discret, number, ind='T'):
        # TODO: donner des noms plus explicite que ind et drap
        if ind=='T':
            drap=1
        else:
            drap=-1
        
        a=drap*self.rack.a
        # TODO: vérifier pourquoi la variable n'est pas utilisée: bug?
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        phi0=a/(self.DFF/2)
        
        # TODO: donner des noms plus explicites que ind et ref
        ref=[]
        if ind=='R':
            theta=npy.linspace(phi0,drap*self.phi_trochoide,discret)
        else:
            theta=npy.linspace(drap*self.phi_trochoide,phi0,discret)
        for t in theta:
            ref.append(vm.Point2D((self._Trochoide(t,ind))))
        ref=primitives2D.RoundedLines2D(ref,{},False)
        ref=ref.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        if ind=='T':
            L1=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
            L1.points[0]=self.rac
        else:
#            L=ref.Rotation(vm.Point2D((0,0)),self.base_circular_tooth_thickness*2/self.base_diameter)
            L1=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
            self.rac=L1.points[-1]

        return L1
    
    def _Trochoide(self,phi,ind='T'):
        # TODO: donner des noms plus explicite que ind, drap et sol
        if ind=='T':
            drap=1
        else:
            drap=-1
        a=drap*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.DFF/2
        rho=self.rack.root_radius
        x2=rho*npy.sin(npy.arctan((a-r*phi)/b)-phi)+a*npy.cos(phi)-b*npy.sin(phi)+r*(npy.sin(phi)-phi*npy.cos(phi))
        y2=-rho*npy.cos(npy.arctan((a-r*phi)/b)-phi)-a*npy.sin(phi)-b*npy.cos(phi)+r*(npy.cos(phi)+phi*npy.sin(phi))
        sol=(y2,x2)
        return sol
    
    def _RootCircleTrace(self,number):
        # TODO: A quoi sert cette fonction??
        # TODO: nommage: Trace ->Plot ?
        # TODO: donner un nom plus explicite que drap
        drap=1
        a=drap*self.rack.a
        phi0=a/(self.DFF/2)
        p1=vm.Point2D((self._Trochoide(phi0,'T')))
        p1=p1.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        drap=-1
        a=drap*self.rack.a
        phi0=a/(self.DFF/2)
        p2=vm.Point2D((self._Trochoide(phi0,'R')))
        p2=p2.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        ref=primitives2D.RoundedLines2D([p1,p2],{},False)
        L2=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
        return L2
    
    def _OutsideTrace(self,number):
        # TODO: Commentaire en anglais, Nom plus explicite pour la fonction?
        #trace du sommet des dents en arc de cercle
        theta4=npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter
        p1=vm.Point2D((self.outside_diameter/2*npy.cos(theta4),self.outside_diameter/2*npy.sin(theta4)))
        p2=p1.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        p3=p2.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        #ref=vm.Arc2D(p1,p2,p3)
        ref=primitives2D.RoundedLines2D([p3,p2,p1],{},False)
        L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.Z)
        return L

    
    def _ISO_YS(self,s_thickness_iso):
        # TODO: docstring en anglais
        #facteur de concentration de contrainte en pied de dent
        rho_f=self.rack.root_radius+self.rack.b**2/(self.DFF/2+self.rack.b)
        coeff_ys_iso=1+0.15*s_thickness_iso/rho_f
        return coeff_ys_iso
    
    def GearISOSection(self,angle):
        a=self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.DFF/2
        theta0 = fsolve((lambda theta:a + b*npy.tan(theta) + r*(-angle - self.root_angle/2 - theta + npy.pi/2)) ,0)[0]
        phi0=(a-b*npy.tan(theta0))/r
        pt_iso = self._Trochoide(phi0)
        angle0 = npy.arctan(pt_iso[1]/pt_iso[0])-self.root_angle/2
        angle_iso = self.root_gear_angle-2*angle0
        diameter_iso = 2*norm(pt_iso)
        s_thickness_iso = diameter_iso*npy.sin(angle_iso/2)
        h_height_iso = (s_thickness_iso/2)/npy.tan(angle)
        return s_thickness_iso, h_height_iso
        
    
    def Dict(self):
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


class MeshAssembly():
    def __init__(self,Z, center_distance, connections,transverse_pressure_angle,
                 coefficient_profile_shift,gear_graph, transverse_pressure_angle_rack,
                 coeff_gear_addendum,coeff_gear_dedendum,coeff_root_radius,
                 coeff_circular_tooth_thickness,list_gear,material,torque,cycle,safety_factor):
        
        self.center_distance=center_distance
        self.connections=connections
        self.transverse_pressure_angle=transverse_pressure_angle
        # TODO pourquoi l'utilisateur devrait donner un graphe? Ne peut on pas le calculer à partir des connections?
        self.gear_graph=gear_graph
        self.list_gear=list_gear
        self.material=material
        self.helix_angle=0

        # TODO: DF devrait être attaché aux objets self.meshes
        self.DF, DB, self.connections_dfs = self.GearGeometryParameter(Z)
        self.cycle = self.CycleParameter(cycle, Z)
        self.torque1, self.torque2, self.normal_load, self.tangential_load, self.radial_load = self.GearTorque(Z, torque, DB)
        
        # TODO: devrait etre une simple liste
        #instantiation des objets Gears
        self.meshes={}
        for ne,ns in enumerate(self.connections):
            self.meshes[ne]={}
            for ng in ns:
                z=Z[ng]
                db=DB[ne][ng]
                cp=coefficient_profile_shift[ng]
                ngp=self.list_gear.index(ng)
                tpa=transverse_pressure_angle_rack[ngp]
                cga=coeff_gear_addendum[ngp]
                cgd=coeff_gear_dedendum[ngp]
                crr=coeff_root_radius[ngp]
                cct=coeff_circular_tooth_thickness[ngp]
                self.meshes[ne][ng]=Mesh(z,db,cp,tpa,cga,cgd,crr,cct)
                
        self.linear_backlash,self.radial_contact_ratio = \
            self.GearContactRatioParameter(Z, coefficient_profile_shift, DB,
                                           transverse_pressure_angle_rack,
                                           coeff_gear_addendum, coeff_gear_dedendum,
                                           coeff_root_radius,
                                           coeff_circular_tooth_thickness)
        self.gear_width,self.sigma_iso,self.sigma_lim=self.GearWidthDefinition(safety_factor)
            
    def Update(self, Z, center_distance, connections, transverse_pressure_angle,
               coefficient_profile_shift, gear_graph,
               transverse_pressure_angle_rack, coeff_gear_addendum,
               coeff_gear_dedendum, coeff_root_radius,coeff_circular_tooth_thickness,
               list_gear,material, torque,cycle, safety_factor):
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.DF,DB,self.connections_dfs=self.GearGeometryParameter(Z)
        self.linear_backlash,self.radial_contact_ratio = self.GearContactRatioParameter(
            Z, coefficient_profile_shift, DB, transverse_pressure_angle_rack,
            coeff_gear_addendum, coeff_gear_dedendum, coeff_root_radius,
            coeff_circular_tooth_thickness)
        
    def GearContactRatioParameter(self,Z,coefficient_profile_shift,DB,
                           transverse_pressure_angle_rack,coeff_gear_addendum,
                           coeff_gear_dedendum,coeff_root_radius,
                           coeff_circular_tooth_thickness):
        
        for ne,ns in enumerate(self.connections):
            for ng in ns:
                z=Z[ng]
                db=DB[ne][ng]
                cp=coefficient_profile_shift[ng]
                ngp=self.list_gear.index(ng)
                tpa=transverse_pressure_angle_rack[ngp]
                cga=coeff_gear_addendum[ngp]
                cgd=coeff_gear_dedendum[ngp]
                crr=coeff_root_radius[ngp]
                cct=coeff_circular_tooth_thickness[ngp]
                self.meshes[ne][ng].Update(z,db,cp,tpa,cga,cgd,crr,cct)
            
        linear_backlash=[]
        radial_contact_ratio=[]
        for g1,g2 in self.connections_dfs:
            ne=self.connections.index((g1,g2))
            circular_tooth_thickness1=self.meshes[ne][g1].GearSection(self.DF[ne][g1])
            circular_tooth_thickness2=self.meshes[ne][g2].GearSection(self.DF[ne][g2])
            transverse_radial_pitch1=npy.pi*self.DF[ne][g1]/self.meshes[ne][g1].Z
            space_width1=transverse_radial_pitch1-circular_tooth_thickness1
            space_width2=transverse_radial_pitch1-circular_tooth_thickness2    
            linear_backlash.append(min(space_width1-circular_tooth_thickness2,space_width2-circular_tooth_thickness1))
            transverse_pressure_angle1=self.transverse_pressure_angle[ne]
            center_distance1=self.center_distance[ne]
            radial_contact_ratio.append((1/2*(npy.sqrt(self.meshes[ne][g1].outside_diameter**2
                                                       - self.meshes[ne][g1].DB**2)
                                              + npy.sqrt(self.meshes[ne][g2].outside_diameter**2
                                                         - self.meshes[ne][g2].DB**2)
                                        - 2*center_distance1*npy.sin(transverse_pressure_angle1))
                                        /(transverse_radial_pitch1*npy.cos(transverse_pressure_angle1))))
        return linear_backlash,radial_contact_ratio

    def GearGeometryParameter(self,Z):
        DF={}
        DB={}
#        transverse_pressure_angle=[self.transverse_pressure_angle_0]
        for ne,(gs,cd) in enumerate(zip(self.connections,self.center_distance)):
            Z1=Z[gs[0]]
            Z2=Z[gs[1]]
            #définition DF
            DF1=2*cd*Z1/Z2/(1+Z1/Z2)
            DF2=2*cd-DF1
            DF[ne]={}
            DF[ne][gs[0]]=DF1
            DF[ne][gs[1]]=DF2
            DB1=float(DF1*npy.cos(self.transverse_pressure_angle[ne]))
            DB2=float(DF2*npy.cos(self.transverse_pressure_angle[ne]))
            DB[ne]={}
            DB[ne][gs[0]]=DB1
            DB[ne][gs[1]]=DB2
        
        connections_dfs=list(nx.edge_dfs(self.gear_graph, [self.connections[0][0],self.connections[0][1]]))

        return DF,DB,connections_dfs
    
    def GearTorque(self,Z,torque,DB):
        ggd=self.gear_graph.degree(self.list_gear)
        liste_node_init=[]
        for ne,nb_connexion in ggd:
            if (nb_connexion==1) and (torque[ne]!='output'):
                liste_node_init.append(ne)
        for ne,tq in torque.items():
            if tq=='output':
                node_output=ne
        if 'output' not in torque.values():
            path_list=[nx.shortest_path(self.gear_graph,source=liste_node_init[0],target=liste_node_init[1])]
        else:
            path_list=[nx.shortest_path(self.gear_graph,source=liste_node_init[0],target=node_output)]
            for nd_init in liste_node_init[1:]:
                path_list.append(nx.shortest_path(self.gear_graph,source=nd_init,target=node_output))
        torque1={}
        normal_load={}
        tangential_load={}
        radial_load={}
        for node_init,path in zip(liste_node_init,path_list):
            torque_moteur_m=torque[node_init]
            for i,eng1 in enumerate(path[0:-1]):
                eng2=path[i+1]
                torque_recepteur=-torque_moteur_m*Z[eng2]/Z[eng1]
                if (eng2 in torque.keys()) and (torque[eng2]!='output'):
                    torque_recepteur+=torque[eng2]
                torque1[eng1]=torque_moteur_m
                torque2={}
                torque2[eng2]=torque_recepteur
                # TODO: Vérifier pourquoi la variable n'est pas utilisée
                torque_input_m=torque_recepteur
        for node_init,path in zip(liste_node_init,path_list):
            torque_moteur_m=torque[node_init]
            for i,eng1 in enumerate(path[0:-1]):
                eng2=path[i+1]
                if (eng1, eng2) in self.connections:
                    ne=self.connections.index((eng1, eng2))
                else:
                    ne=self.connections.index((eng2, eng1))

                tq1=torque1[eng1]
                normal_load[ne]=abs(tq1)*2/(DB[ne][eng1])
                tangential_load[ne]=abs(tq1)*2/(self.DF[ne][eng1])
                radial_load[ne]=npy.tan(self.transverse_pressure_angle[ne])*tangential_load[ne]
        return torque1,torque2,normal_load,tangential_load,radial_load
    
    def CycleParameter(self,cycle,Z):
        eng_init=list(cycle.keys())[0]
        for eng in self.list_gear:
            if eng not in cycle.keys():
                cycle[eng]=cycle[eng_init]*Z[eng_init]/Z[eng]
        return cycle
                    
    ### Stress
    
    def SigmaLewis(self):
        
        self.sigma_lewis_maximum1=6*self.tangential_load*self.gear_height_lewis1/(self.gear_width*self.Gear1.root_gear_length**2)
        self.sigma_lewis_maximum2=6*self.tangential_load*self.gear_height_lewis2/(self.gear_width*self.Gear2.root_gear_length**2)
        
    def GearWidthDefinition(self,safety_factor):
        coeff_yf_iso=self._CoeffYFIso()
        coeff_ye_iso=self._CoeffYEIso()
        coeff_yb_iso=self._CoeffYBIso()
        sigma_lim=self.SigmaMaterialISO(safety_factor)
        gear_width={}
        for eng in self.list_gear:
            gear_width[eng]=0
        for ne,(eng1,eng2) in enumerate(self.connections):
            gear_width1=abs(self.tangential_load[ne]
                            / (sigma_lim[ne][eng1]
                                * self.meshes[ne][eng1].rack.module)
                            *coeff_yf_iso[ne][eng1]
                            *coeff_ye_iso[ne]
                            *coeff_yb_iso[ne][eng1])
                            
            gear_width2=abs(self.tangential_load[ne]
                            /(sigma_lim[ne][eng2]
                                *self.meshes[ne][eng2].rack.module)
                            *coeff_yf_iso[ne][eng2]
                            *coeff_ye_iso[ne]
                            *coeff_yb_iso[ne][eng2])
                            
            gear_width_set=max(gear_width1,gear_width2)
            gear_width[eng1]=max(gear_width[eng1],gear_width_set)
            gear_width[eng2]=max(gear_width[eng2],gear_width_set)
        sigma_iso=sigma_lim
        return gear_width,sigma_iso,sigma_lim
        
    def SigmaMaterialISO(self,safety_factor):
        angle=30/180*npy.pi
        sigma_lim={}
        for ne,(eng1,eng2) in enumerate(self.connections):
            sigma_lim[ne]={}
            
            matrice_wholer=self.material[eng1].data_wholer_curve
            matrice_material=self.material[eng1].data_gear_material
            sgla=self.material[eng1].FunCoeff(self.cycle[eng1],npy.array(matrice_wholer['data']),matrice_wholer['x'],matrice_wholer['y'])
            sgl1=self.material[eng1].FunCoeff(sgla,npy.array(matrice_material['data']),matrice_material['x'],matrice_material['y'])
            s_thickness_iso_1,h_height_iso_1=self.meshes[ne][eng1].GearISOSection(angle)
            coeff_ys_iso=self.meshes[ne][eng1]._ISO_YS(s_thickness_iso_1)
            sigma_lim[ne][eng1]=(sgl1/(safety_factor*coeff_ys_iso))*10**7
            
            matrice_wholer=self.material[eng2].data_wholer_curve
            matrice_material=self.material[eng2].data_gear_material
            sglb=self.material[eng2].FunCoeff(self.cycle[eng2],npy.array(matrice_wholer['data']),matrice_wholer['x'],matrice_wholer['y'])
            sgl2=self.material[eng2].FunCoeff(sglb,npy.array(matrice_material['data']),matrice_material['x'],matrice_material['y'])
            s_thickness_iso_2,h_height_iso_2=self.meshes[ne][eng2].GearISOSection(angle)
            coeff_ys_iso=self.meshes[ne][eng2]._ISO_YS(s_thickness_iso_2)
            sigma_lim[ne][eng2]=(sgl2/(safety_factor*coeff_ys_iso))*10**7
        return sigma_lim

        
    def _CoeffYFIso(self):
        #facteur de forme pour la contrainte ISO
        angle=30/180*npy.pi
        coeff_yf_iso={}
        for ne,(eng1,eng2) in enumerate(self.connections):
            coeff_yf_iso[ne]={}
            s_thickness_iso_1,h_height_iso_1=self.meshes[ne][eng1].GearISOSection(angle)
            s_thickness_iso_2,h_height_iso_2=self.meshes[ne][eng2].GearISOSection(angle)
            coeff_yf_iso[ne][eng1]=((6*(h_height_iso_1/self.meshes[ne][eng1].rack.module)*npy.cos(self.transverse_pressure_angle[ne]))
                                    /((s_thickness_iso_1/self.meshes[ne][eng1].rack.module)**2
                                       *npy.cos(self.meshes[ne][eng1].rack.transverse_pressure_angle)))
            coeff_yf_iso[ne][eng2]=((6*(h_height_iso_2/self.meshes[ne][eng2].rack.module)*npy.cos(self.transverse_pressure_angle[ne]))
                                    /((s_thickness_iso_2/self.meshes[ne][eng2].rack.module)**2
                                       *npy.cos(self.meshes[ne][eng2].rack.transverse_pressure_angle)))
        return coeff_yf_iso
        
    def _CoeffYEIso(self):
        #facteur de conduite pour la contrainte ISO
        coeff_ye_iso=[]
        for ne,eng in enumerate(self.connections):
            coeff_ye_iso.append(1/self.radial_contact_ratio[ne])
        return coeff_ye_iso
        
    def _CoeffYBIso(self):
        #facteur de contrefort pour la contrainte ISO
        coeff_yb_iso={}
        for ne,(eng1,eng2) in enumerate(self.connections):
            coeff_yb_iso[ne]={}
            matrice_YB=self.material[eng1].data_coeff_YB_Iso
            coeff_yb_iso[ne][eng1]=self.material[eng1].FunCoeff(self.helix_angle,npy.array(matrice_YB['data']),matrice_YB['x'],matrice_YB['y'])
            matrice_YB=self.material[eng2].data_coeff_YB_Iso
            coeff_yb_iso[ne][eng2]=self.material[eng2].FunCoeff(self.helix_angle,npy.array(matrice_YB['data']),matrice_YB['x'],matrice_YB['y'])
        return coeff_yb_iso
        
    ### Fonction de trace et export
    
    def GearRotate(self,list_gear,list_center,list_rot):
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
          
        Angle1=npy.arccos(self.meshes[set_pos][liste_eng[0]].DB/self.DF[set_pos][liste_eng[0]])
        Angle2=npy.arccos(self.meshes[set_pos][liste_eng[1]].DB/self.DF[set_pos][liste_eng[1]])
        Gear1Angle=-(npy.tan(Angle1)-Angle1)
        Gear2Angle=-(npy.tan(Angle2)-Angle2)+npy.pi        
        return [Gear1Angle,Gear2Angle]
    
    def VolumeModel(self, centers = {}, axis = (1,0,0), name = ''):
        
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
            # TODO: Check if this is a bug
            engr_pos=[self.list_gear.index(eng1),self.list_gear.index(eng2)]
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
                Gears3D[eng1]=self.meshes[set_pos][eng1].Contour(3)
            Struct.append(vm.Circle2D(vm.Point2D(position1),self.DF[set_pos][eng1]/2))
            Gears3D[eng2]=self.meshes[set_pos][eng2].Contour(3)
            Struct.append(vm.Circle2D(vm.Point2D(position2),self.DF[set_pos][eng2]/2))            
            
            if position2[1]==position1[1]:
                if position2[2]-position1[2]>0:
                    angle0=npy.pi/2
                    print(1)
                else:
                    angle0=-npy.pi/2
                    print(2)
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
                Rotation[set_pos][eng2]=list_rot[1]-angle0-delta_rot*((self.meshes[set_pos][eng1].Z)/(self.meshes[set_pos][eng2].Z))
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
        # Mass function: now only DF*pi*width TODO: improve
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
        
        model = self.VolumeModel(centers, axis)
        model.FreeCADExport(python_path ,file_path, freecad_path, export_types)
        
    def PosAxis(self,position):
        # TODO: Que fait cette fonction?
        #optimisation pour le placement des axes des engrenages
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
        x_opt=self.PosAxis(position)
        print(x_opt)
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
                TG[en[0]]=self.meshes[ne][en[0]].Contour(5)
            Struct.append(vm.Circle2D(vm.Point2D(position1),self.DF[ne][en[0]]/2))
            TG[en[1]]=self.meshes[ne][en[1]].Contour(5)
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
                Rot[ne][en[1]]=list_rot[1]-angle-delta_rot*((self.meshes[ne][en[0]].Z)/(self.meshes[ne][en[1]].Z))
            sol=self.GearRotate([TG[en[0]],TG[en[1]]],[position1,position2],list_rot=[Rot[ne][en[0]],Rot[ne][en[1]]])
            if num==0:
                L1.extend(sol[0])
            L1.extend(sol[1])
        L1.extend(Struct)
        G1=vm.Contour2D(L1)
        G1.MPLPlot()
        
    def Dict(self):
        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv==npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=float(v)
            else:
                d[k]=v
        
#        gears_dicts = {}
#        for ne,gs in enumerate(self.connections):
#            d['Gear'+str(gs[0]+1)]=self.gears[ne][gs[0]].Dict()
#            d['Gear'+str(gs[1]+1)]=self.gears[ne][gs[1]].Dict()
        del d['meshes']
        del d['gear_graph']
        del d['material']
        return d

class ContinuousMeshesAssemblyOptimizer:
    def __init__(self, Z, center_distance, connections, transverse_pressure_angle,
                 coefficient_profile_shift, gear_graph,cond_init,
                 rack_list, rack_choice, list_gear, material,torque,
                 cycle, safety_factor):
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.coefficient_profile_shift=coefficient_profile_shift
        self.cond_init=cond_init
        self.rack_list=rack_list
        self.rack_choice=rack_choice
        self.list_gear=list_gear
        Bounds=list(center_distance)
        Bounds.extend(transverse_pressure_angle)
        # TODO: Check if this is a bug
        tp=len(coefficient_profile_shift.keys())
        for i in self.list_gear:
            Bounds.append(coefficient_profile_shift[i])
        self.solutions=[]
        
        self.xi={'Z': Z, 'connections': connections,'gear_graph':gear_graph,'list_gear':list_gear,
                 'material':material,'torque':torque,'cycle':cycle,'safety_factor':safety_factor}
        self.xj,self.dict_xu=self._init()
        self.xt=dict(list(self.xi.items())+list(self.xj.items()))
        
        for k,v in self.dict_xu:
            Bounds.append(self.rack_list[v][k])
        self.Bounds=npy.array(Bounds)
        
        #xk partie non optimisée du vecteur x
        #xu partie optimisée du vecteur x
        
        self.MeshAssembly = MeshAssembly(**self.xt)
        
    # TODO: Pourquoi un deuxième Init ???
    def _init(self):
        xj={'center_distance':self._init_list(self.center_distance),
                 'transverse_pressure_angle':self._init_list(self.transverse_pressure_angle),
                 'coefficient_profile_shift':self._init_item(self.coefficient_profile_shift),
                 'transverse_pressure_angle_rack':[],'coeff_gear_addendum':[],
                 'coeff_gear_dedendum':[],'coeff_root_radius':[],'coeff_circular_tooth_thickness':[]}
        dict_xu=[]
        x_list={}
        for i in list(set(list(self.rack_choice.values()))):
            x_list[i]={'transverse_pressure_angle_rack':[],
                  'coeff_gear_addendum':[],
                  'coeff_gear_dedendum':[],
                  'coeff_root_radius':[],
                  'coeff_circular_tooth_thickness':[]}
            for k,v in self.rack_list[i].items():
                if k not in ['type','name','module']:
                    if v[0]==v[1]:
                        x_list[i][k]=v[0]
                    else:
                        x_list[i][k]=(v[1]-v[0])*float(npy.random.random(1))+v[0]
                        dict_xu.append((k,i))
#        for k in sorted(list(self.rack_choice.keys())):
        for k in self.list_gear:
            v=self.rack_choice[k]
            for k2,v2 in self.rack_list[v].items():
                if k2 not in ['type','name','module']:
                    xj[k2].append(x_list[v][k2])
        return xj,dict_xu
    
    def _convert_xj2Xu(self,xj):
        sol=xj['center_distance'].copy()
        sol.extend(xj['transverse_pressure_angle'])
        # TODO: Check if this is a bug
        tp=len(xj['coefficient_profile_shift'].keys())
        for key in self.list_gear:
            sol.append(xj['coefficient_profile_shift'][key])
        for k,v in self.dict_xu:
            for k2 in sorted(list(self.rack_choice.keys())):
                v2=self.rack_choice[k2]
                if v2==v:
                    ind=k2
            sol.append(xj[k][ind])
        return npy.array(sol)

    def _convert_Xu2xj(self,X):
        X=list(X)
        tp1=len(self.center_distance)
        tp2=len(self.transverse_pressure_angle)
        tp3=len(self.coefficient_profile_shift.keys())
        xj=self.xj
        xj['center_distance']=X[0:tp1]
        xj['transverse_pressure_angle']=X[tp1:tp1+tp2]
        xj['coefficient_profile_shift']={}
        for i in self.list_gear:
            j=self.list_gear.index(i)
            xj['coefficient_profile_shift'][i]=X[tp1+tp2+j]
#        xj['transverse_pressure_angle_rack']=[0]*len(self.rack_choice.keys())
#        xj['coeff_gear_addendum']=[0]*len(self.rack_choice.keys())
#        xj['coeff_gear_dedendum']=[0]*len(self.rack_choice.keys())
#        xj['coeff_root_radius']=[0]*len(self.rack_choice.keys())
#        xj['coeff_circular_tooth_thickness']=[0]*len(self.rack_choice.keys())
        for i,(k,v) in enumerate(self.dict_xu):
            for k2 in sorted(list(self.rack_choice.keys())):
                v2=self.rack_choice[k2]
                if v2==v:
                    ind=k2
                    xj[k][ind]=X[tp1+tp2+tp3+i]
        return xj

    def _init_list(self,list):
        list1=[]
        for li in list:
            list1.append((li[1]-li[0])*float(npy.random.random(1))+li[0])
        return list1
    
    def _init_item(self,dico):
        dict1={}
        for k1,v1 in dico.items():
            dict1[k1]=(v1[1]-v1[0])*float(npy.random.random(1))+v1[0]
        return dict1
    
    def Update(self,xj):
        self.xt=dict(list(self.xi.items())+list(xj.items()))
        self.MeshAssembly.Update(**self.xt)
        return xj
    
    # TODO: implémenter les critères comme des calculs dans le MeshAssembly
    # Les appeller ensuite. Pour un exemple regarder InternalInterference du model3D PWT.
    def Fineq(self,X):
        xj=self._convert_Xu2xj(X)
        xj=self.Update(xj)
        ineq=[]
        #jeu inter-denture
        for lb in self.MeshAssembly.linear_backlash:
            ineq.append(lb)
            ineq.append((4e-4)-lb)
        #contrainte géométrique
        for ne,gs in enumerate(self.MeshAssembly.connections):
            dia1=self.MeshAssembly.meshes[ne][gs[0]].root_diameter_active
            dia2=self.MeshAssembly.meshes[ne][gs[1]].root_diameter_active
            de1=self.MeshAssembly.meshes[ne][gs[0]].outside_diameter
            de2=self.MeshAssembly.meshes[ne][gs[1]].outside_diameter
            cd=self.MeshAssembly.center_distance[ne]
            ineq.append(cd-(de1/2+dia2/2))
            ineq.append(cd-(de2/2+dia1/2))
            oaa1=self.MeshAssembly.meshes[ne][gs[0]].outside_active_angle
            oaa2=self.MeshAssembly.meshes[ne][gs[1]].outside_active_angle
            ineq.append(oaa1)
            ineq.append(oaa2)
            df1=self.MeshAssembly.DF[ne][gs[0]]
            df2=self.MeshAssembly.DF[ne][gs[1]]
            db1=self.MeshAssembly.meshes[ne][gs[0]].DB
            db2=self.MeshAssembly.meshes[ne][gs[1]].DB
            ineq.append(df1-db1)
            ineq.append(df2-db2)
        #contrainte sur le RCA
        for ne,gs in enumerate(self.MeshAssembly.connections):
            rca=self.MeshAssembly.radial_contact_ratio[ne]
            ineq.append(rca-1)
        #contrainte sur le module
        for ne,gs in enumerate(self.MeshAssembly.connections):
            for g in gs:
                mo=self.MeshAssembly.meshes[ne][g].rack.module
                list_module=self.rack_list[self.rack_choice[g]]['module']
#                if lm[0]<lm[1]:
                ineq.append(mo-list_module[0])
                ineq.append(list_module[1]-mo)
        return ineq
    
    def Feq(self,X):
        x=self._convert_Xu2xj(X)
        x=self.Update(x)
        eq=[]
#        for ng in range(self.MeshAssembly.gear_graph.number_of_nodes()):
        for ng in self.list_gear:
            nel=list(self.MeshAssembly.gear_graph.edges(ng))
            if len(nel)>1:
                list_db=[]
                for ne in nel:
                    if (ne[0],ne[1]) in self.MeshAssembly.connections:
                        nes=self.MeshAssembly.connections.index((ne[0],ne[1]))
                    elif (ne[1],ne[0]) in self.MeshAssembly.connections:
                        nes=self.MeshAssembly.connections.index((ne[1],ne[0]))
                    list_db.append(self.MeshAssembly.meshes[nes][ng].DB)
#                    print(self.MeshAssembly.meshes[nes][ng].DB,nes,ng,list_db)
                list1=itertools.combinations(list_db,2)
                for n1,n2 in list1:
                    eq.append(n1-n2)
#        for ne,gs in enumerate(self.MeshAssembly.connections):
#            for g in gs:
#                mo=self.MeshAssembly.meshes[ne][g].rack.module
#                lm=self.rack_list[self.rack_choice[g]]['module']
#                if lm[0]==lm[1]:
#                    eq.append(mo-lm[0])
        if eq==[]:
            eq=[0]
        return eq
        
    def Objective(self,X):
        x=self._convert_Xu2xj(X)
        x=self.Update(x)
        fineq=self.Fineq(X)
        feq=self.Feq(X)
        obj=0
        #Maximisation du module pour avoir des pignons avec un faible gear_width
        for ne,gs in enumerate(self.MeshAssembly.connections):
            for g in gs:
                mo=self.MeshAssembly.meshes[ne][g].rack.module
                list_module=self.rack_list[self.rack_choice[g]]['module']
                if list_module[0]<list_module[1]:
                    obj+=100*(list_module[1]-mo)**2
                
        #Minimisation des entraxes sur la borne inf
        for num_engr,list_cd in enumerate(self.center_distance):
            obj+=100*(list_cd[0]-self.MeshAssembly.center_distance[num_engr])**2
            
        for lb in self.MeshAssembly.linear_backlash:
            obj+=100*(lb)**2
            
        for i in fineq:
            if i < 0:
                obj+=1000*i**2
            else:
                obj+=0.000001*i
        for i in feq:
            if i<0:
                obj+=1000*i**2
            else:
                obj+=1000*i**2
        return obj
    
    def Optimize(self, verbose = False):
        max_iter=5
        i=0
        arret=0 
        while i<max_iter and arret==0:
            xj0,dict_xu=self._init()
            xj0=self.Update(xj0)
            X0=self._convert_xj2Xu(xj0)
            if self.Feq(X0)==[0]:
                cons = {'type': 'ineq','fun' : self.Fineq}
            else:
                cons = ({'type': 'eq','fun' : self.Feq},{'type': 'ineq','fun' : self.Fineq})
            cx = minimize(self.Objective, X0, bounds=self.Bounds,constraints=cons)
            Xsol=cx.x
            xsol=self._convert_Xu2xj(Xsol)
            xsol=self.Update(xsol)
            if verbose:
                print('Iteration n°{} with status {}, min(fineq):{}, max(eq):{}'.format(i,
                      cx.status,min(self.Fineq(Xsol)),max(npy.abs(self.Feq(Xsol)))))
            if min(self.Fineq(Xsol))>-1e-5 and max(npy.abs(self.Feq(Xsol)))<1e-5:
                self.solutions.append(xsol)
                arret=1
            i=i+1


class MeshAssemblyOptimizer:
    def __init__(self, connections, gear_speed, center_distance, Z=None,
                 transverse_pressure_angle=None, helix_angle=None,
                 gear_width=None, frequency=[[0,0]],
                 coefficient_profile_shift=None, rack_list=None,
                 rack_choice=None, material=None, torque=None, cycle=None,
                 safety_factor=1):
        # Valeur par defaut
        list_gear=[]
        for gs in connections:
            for g in gs:
                if g not in list_gear:
                    list_gear.append(g)
        # TODO: Check if this is a bug
        nb_gear=len(list_gear)
        nb_set=len(connections)
                    
        if transverse_pressure_angle==None:
            transverse_pressure_angle=[]
            for i in range(nb_set):
                transverse_pressure_angle.append([15/180*npy.pi,30/180*npy.pi])
            
        if helix_angle==None:
            helix_angle={list_gear[0]:[15/180*npy.pi,25/180*npy.pi]}
        
        if gear_width==None:
            gear_width={list_gear[0]:[15*1e-3,25*1e-3]}
        gw_min=npy.inf
        gw_max=-npy.inf
        for ne in gear_width.keys():
            if gear_width[ne][0]<gw_min:
                gw_min=gear_width[ne][0]
            if gear_width[ne][1]>gw_max:
                gw_max=gear_width[ne][1]
        for ne in list_gear:
            if ne not in gear_width.keys():
                gear_width[ne]=[gw_min,gw_max]
                
        if coefficient_profile_shift==None:
            coefficient_profile_shift={list_gear[0]:[-0.8,0.8]}
        for ne in list_gear:
            if ne not in coefficient_profile_shift.keys():
                coefficient_profile_shift[ne]=[-0.8,0.8]
                
        if rack_list==None:
            rack_list={0:{'name':'Optim_Module','module':[1*1e-3,3*1e-3],
                          'transverse_pressure_angle_rack':[20*npy.pi/180,20*npy.pi/180],
                          'coeff_gear_addendum':[1,1],
                          'coeff_gear_dedendum':[1.25,1.25],
                          'coeff_root_radius':[0.38,0.38],
                          'coeff_circular_tooth_thickness':[0.5,0.5]}}
                        
            
        if rack_choice==None:
            rack_choice={list_gear[0]:[list(rack_list.keys())[0]]}
        for ne in list_gear:
            if ne not in rack_choice.keys():
                rack_choice[ne]=[list(rack_list.keys())[0]]
                
        if material==None:
            material={list_gear[0]:hardened_alloy_steel}
        for ne in list_gear:
            if ne not in material.keys():
                material[ne]=hardened_alloy_steel
        
        if torque==None:
            torque={list_gear[0]:100,list_gear[1]:'output'}
            
        if cycle==None:
            cycle={list_gear[0]:1e6}
        
        if Z==None:
            Z={}
        self.Z=Z
        self.connections=connections
        self.gear_speed=gear_speed
        self.frequency=frequency
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.coefficient_profile_shift=coefficient_profile_shift
        self.rack_list=rack_list
        self.rack_choice=rack_choice
        self.list_gear=list_gear
        self.material=material
        self.torque=torque
        self.cycle=cycle
        self.safety_factor=safety_factor

        self.nb_gear=len(list_gear)
        gear_graph=nx.Graph()
        gear_graph.add_nodes_from(list_gear)
        gear_graph.add_edges_from(self.connections)
        self.gear_graph=gear_graph
        
        self.nb_rack=len(self.rack_list.keys())

        self.node_init=int(list(self.gear_speed.keys())[0])
        self.connections_dfs=list(nx.dfs_edges(gear_graph,self.node_init))
        
        if self.Z=={}:
            var_Z=self.AnalyseZ()
            self.Z=var_Z

        self.AnalyzeCombination()
        
        for i,plex in enumerate(self.plex_calcul):
            plex['gear_graph']=self.gear_graph
            plex['rack_list']=self.rack_list
            plex['list_gear']=self.list_gear
            plex['material']=self.material
            plex['torque']=self.torque
            plex['cycle']=self.cycle
            plex['connections']=self.connections
#            plex['center_distance']=self.center_distance
            plex['transverse_pressure_angle']=self.transverse_pressure_angle
            plex['coefficient_profile_shift']=self.coefficient_profile_shift
            plex['safety_factor']=safety_factor
            self.plex_calcul[i]=plex
            
        self.solutions=[]
        self.solutions_search=[]
        self.analyse=[]
        
    def AnalyseZ(self):
        #nombre de dents adaptatif
        Z=self.Z
        for i,(engr1,engr2) in enumerate(self.connections):
            cd_min=self.center_distance[i][0]
            cd_max=self.center_distance[i][1]
            module1_min,module1_max=(npy.inf,0)
            for rack_num in self.rack_choice[engr1]:
                mod_min,mod_max=self.rack_list[rack_num]['module']
                module1_min,module1_max=(min(module1_min,mod_min),max(module1_max,mod_max))
            module2_min,module2_max=(npy.inf,0)
            for rack_num in self.rack_choice[engr2]:
                mod_min,mod_max=self.rack_list[rack_num]['module']
                module2_min,module2_max=(min(module2_min,mod_min),max(module2_max,mod_max))
            demul_min=self.gear_speed[engr1][0]/self.gear_speed[engr2][1]
            demul_max=self.gear_speed[engr1][1]/self.gear_speed[engr2][0]
            DF1_max=2*cd_max/(1+demul_min)
            Z1_max=int(DF1_max/module1_min)+1
            DF2_max=2*cd_max*demul_max/(1+demul_max)
            Z2_max=int(DF2_max/module2_min)+1
            DF1_min=2*cd_min/(1+demul_max)
            Z1_min=int(DF1_min/module1_max)-1
            DF2_min=2*cd_min*demul_min/(1+demul_min)
            Z2_min=int(DF2_min/module2_max)-1
            
            if engr1 not in Z.keys():
                Z[engr1]=[Z1_min,Z1_max]
            else:
                Z[engr1]=[max(Z1_min,Z[engr1][0]),min(Z1_max,Z[engr1][1])]
            if engr2 not in Z.keys():
                Z[engr2]=[Z2_min,Z2_max]
            else:
                Z[engr2]=[max(Z2_min,Z[engr2][0]),min(Z2_max,Z[engr2][1])]
        return Z

    def AnalyzeCombination(self, verbose = False):
        n1=self.node_init
        list_node=[n1]
        for (n1,n2) in self.connections_dfs:
            if n2 not in list_node:
                list_node.append(n2)
                
        np=[]
        list_gear=[]
        for engr_num in list_node:
            np.append(self.Z[engr_num][1]-self.Z[engr_num][0]+1)
            list_gear.append(npy.arange(self.Z[engr_num][0],self.Z[engr_num][1]+1))
        np.extend([self.nb_rack]*self.nb_gear)

        list_rack=list(self.rack_list.keys())
#        print(list_gear)
#        print(list_node)
        
        demul_int_min=1/9.
        demul_int_max=9
#        print(np)
        dt=tools.RegularDecisionTree(np)
        
        incr=0
        self.plex_calcul=[]
        self.fonctionnel=[]
        self.fonctionnel_module=[]

        def pgcd(a,b) :
            while a%b != 0 :
                a, b = b, a%b
            return b
        while not dt.finished:
            valid=True
#            print(dt.current_node)
            Z_node = [list_gear[i][node_value] for i,node_value in enumerate(dt.current_node[:self.nb_gear])]
#            print('Znode: ', Z_node)
            if (dt.current_depth<=(self.nb_gear-1)) and (dt.current_depth>0):
                z1=list_gear[dt.current_depth-1][dt.current_node[dt.current_depth-1]]
                z2=list_gear[dt.current_depth][dt.current_node[dt.current_depth]]
                #analyse ACV engrenage 2 à 2
                if (pgcd(z1,z2)!=1):
                    valid=False
                #analyse demul interne
                if valid:
                    demul=z1/z2
#                    print(demul)
                    if (demul > demul_int_max) or (demul < demul_int_min):
                        valid=False
                #analyse ACV de l'ensemble des engrenages entre eux
#                if (valid) & (dt.current_depth<=(self.nb_gear-1)):
#                    for n in liste_node[0:dt.current_depth]:
#                        i=liste_node.index(n)
#                        z=liste_gear[dt.current_node[i]]
#                        if pgcd(z,z2)!=1:
#                            valid=False
                #analyse des vitesses du CDC
                if valid:
#                    print('@@@@@@@@@@@@@@@@@@@')
                    v0_min, v0_max = self.gear_speed[list_node[0]]
                    z0 = list_gear[0][dt.current_node[0]]
#                    print('z0', z0)
#                    print(list_node[0:dt.current_depth+1])
                    for engr_index, engr_num in enumerate(list_node[0:dt.current_depth+1]):
#                        print('###')
                        if engr_index>0:# No need of checking input
                            if engr_num in self.gear_speed.keys():
                                z=list_gear[engr_index][dt.current_node[engr_index]]
    #                            print('z', z)
                                demul=z0/z
    #                            print(demul, z, z0)
                                vsi_min=self.gear_speed[engr_num][0] # Specified speed
                                vsi_max=self.gear_speed[engr_num][1]
                                
                                vai_min = v0_min*demul # Actual speed
                                vai_max = v0_max*demul
                                
    #                            vs_min = self.gear_speed
    #                            print(demul)
    #                            print(vsi_min, vsi_max, vai_min, vai_max)
                                msi = 0.5 * (vsi_min+vsi_min)
                                mai = 0.5* (vai_min + vai_max)
                                dsi = vsi_max-vsi_min
                                dai = vai_max-vai_min
                                if abs(msi-mai) > 0.5*(dai+dsi):
                                    valid = False
    #                                print(valid)
                                    break
                #analyse frequence
                if valid:
                    for freq in self.frequency:
                        zm=list_gear[0][dt.current_node[0]]
                        for engr_num,engr_ind in enumerate(dt.current_node):
                            z=list_gear[engr_num][dt.current_node[engr_num]]
                            f_min=(2*npy.pi*v0_min*zm/z)/z
                            f_max=(2*npy.pi*v0_max*zm/z)/z
                            if (max(f_min,f_max)>freq[0]) and (min(f_min,f_max)<freq[1]):
                                valid=False
                                break

            elif (dt.current_depth>(self.nb_gear-1)):
                # Analyse de la faisabilité des cremailleres
                rack_num=list_rack[dt.current_node[-1]]
                rack_pos=list_node[dt.current_depth-self.nb_gear]
                if rack_num not in self.rack_choice[rack_pos]:
                    valid=False
                    
            if (dt.current_depth==(self.nb_gear+self.nb_gear-1)):
                 
                """
                Analyse de la viabilité module/entraxe et construction 
                d'une fonctionnelle mesurant l'écart entre 
                les entraxes estimées et les entraxes mini spécifiées
                """
                liste_DF_min={}
                module_minmax={}
                module_inf,module_sup=(0,npy.inf)
                for tree_pos,tree_val in enumerate(dt.current_node[0:self.nb_gear]):
                    engr_num=list_node[tree_pos]
                    rack_num=list_rack[dt.current_node[tree_pos+self.nb_gear]]
                    module_minmax[engr_num]=self.rack_list[rack_num]['module']
                    module_inf,module_sup=(max(module_inf,module_minmax[engr_num][0]),min(module_sup,module_minmax[engr_num][1]))
                for tree_pos,tree_val in enumerate(dt.current_node[0:self.nb_gear]):
                    z=list_gear[tree_pos][tree_val]
                    engr_num=list_node[tree_pos]
                    liste_DF_min[engr_num]=z*module_inf
                liste_pente_cd_module=[]
                for set_num,(eng1,eng2) in enumerate(self.connections):
                    cd_min=(liste_DF_min[eng1]+liste_DF_min[eng2])/2
                    liste_pente_cd_module.append(cd_min/module_inf)
                cd_minmax_nv=[]
                module_optimal=0
                for set_num,cd in enumerate(self.center_distance):
                    module_optimal=max(module_optimal,cd[0]/liste_pente_cd_module[set_num])
                if module_optimal>module_sup:
                    valid=False
                module_optimal=max(module_optimal,module_inf)
                fonctionnel=0
                for set_num,cd in enumerate(self.center_distance):
                    cd_optimal=liste_pente_cd_module[set_num]*module_optimal
                    if (cd_optimal)>(cd[1]):
                        valid=False
                        break
                    else:
                        fonctionnel+=(cd_optimal-cd[0])**2
                        cd_minmax_nv.append([cd_optimal,min(cd[1],cd_optimal*1.2)])
#            
#                # analyse coherence DB, demul et angle de pression de la cremaillere
#                for set_num,(engr1,engr2) in enumerate(self.connections):
#                    engr1_pos=liste_node.index(engr1)
#                    rack1_num=liste_rack[dt.current_node[self.nb_gear+engr1_pos]]
#                    transverse_pressure_angle=self.rack_list[rack1_num]['transverse_pressure_angle_rack']
#                    Z1=liste_gear[engr1_pos][dt.current_node[engr1_pos]]
#                    DB1_min=npy.cos(transverse_pressure_angle[0])*Z1*module_optimal
#                    DB1_max=npy.cos(transverse_pressure_angle[1])*Z1*module_sup
#                    engr2_pos=liste_node.index(engr2)
#                    rack2_num=liste_rack[dt.current_node[self.nb_gear+engr2_pos]]
#                    transverse_pressure_angle=self.rack_list[rack2_num]['transverse_pressure_angle_rack']
#                    Z2=liste_gear[engr2_pos][dt.current_node[engr2_pos]]
#                    DB2_min=npy.cos(transverse_pressure_angle[0])*Z2*module_optimal
#                    DB2_max=npy.cos(transverse_pressure_angle[1])*Z2*module_sup
#                    if ((Z1/Z2)<(DB1_min/DB2_max)) or ((Z1/Z2)>(DB1_max/DB2_min)):
#                        valid=False
#                        break
        
            if (dt.current_depth==(self.nb_gear+self.nb_gear-1)) & (valid==True):
                gear={}
                rack={}
                for n in list_node:
                    i=list_node.index(n)
                    gear[n]=list_gear[i][dt.current_node[i]]
                    rack[n]=list_rack[dt.current_node[i+self.nb_gear]]
                # TODO: nom plus explicite que Temp!!!!!
                Temp={}
                Temp['Z']=gear
#                Temp['cond_init']=res.x
                Temp['cond_init']=0
                Temp['rack_choice']=rack
                Temp['center_distance']=cd_minmax_nv
                self.fonctionnel.append(fonctionnel)
                self.fonctionnel_module.append(module_optimal)
                self.plex_calcul.append(Temp)
                incr+=1
#            print('valid sent to dt', valid)
            dt.NextNode(valid)
        if incr>1:
            if verbose:
                print('Number of combination found: {}'.format(incr))
        else:
            if verbose:
                print('No teeth combination found: increase center distances')


    def Optimize(self,nb_sol=1, list_sol=None, verbose=False):
        
        compt_nb_sol=0
        if list_sol==None:
            liste_plex=self.plex_calcul
        else:
            liste_plex=[]
            for ind_plex in list_sol:
                liste_plex.append(self.plex_calcul[ind_plex])
            nb_sol=len(liste_plex)
        for plex in liste_plex:
            ga=ContinuousMeshesAssemblyOptimizer(**plex)
            try:
                ga.Optimize()
            # BUG: Définir l'erreur à attraper!
            except ValueError:
                print('Convergence problem')
            if len(ga.solutions)>0:
                xsol=ga.solutions[-1]
                xt=dict(list(ga.xi.items())+list(xsol.items()))
                self.solutions.append(MeshAssembly(**xt))
                compt_nb_sol+=1
                if verbose:
                    print('Mesh sections: {}'.format(self.solutions[-1].gear_width))
                    print('Numbers of teeth: {}'.format(plex['Z']))
                    print('Center distances: {}'.format(self.solutions[-1].center_distance))
                if compt_nb_sol==nb_sol:
                    break

    # TODO: Rename to optimize?
    def SearchOptimumCD(self, nb_sol = 1, verbose = False,
                        progress_callback = lambda x:x):
        #Optimisation des nb_sol meilleures solutions vis à vis de la fonctionnelle
        list_fonctionnel=npy.array(self.fonctionnel)
        
        #En cours de construction
        list_fonctionnel_module=npy.array(self.fonctionnel_module)
        sort_fonct_module=npy.argsort(list_fonctionnel_module)
        compt_fonct_mod=0
        fonct_entraxe=[]
        for ind_plex in sort_fonct_module[::-1]:
            fonct_entraxe.append(list_fonctionnel[ind_plex])
            if compt_fonct_mod==5*nb_sol:
                break
            compt_fonct_mod+=1
        list_fonct_entraxe=npy.array(fonct_entraxe)
        sort_list_fonct_entraxe=npy.argsort(list_fonct_entraxe)
        plex_analyse=[]
        for ind_plex in range(nb_sol):
            plex_analyse.append(sort_fonct_module[::-1][sort_list_fonct_entraxe[ind_plex]])
        
        compt_nb_sol=0
        liste_solutions=[]
#        for num_plex in npy.argsort(list_fonctionnel):
        for num_plex in plex_analyse:
            plex=self.plex_calcul[num_plex]
            ga=ContinuousMeshesAssemblyOptimizer(**plex)
            try:
                ga.Optimize(verbose = verbose)
            # BUG: Définir l'erreur
            except ValueError:
                print('Convergence Problem')
            if len(ga.solutions)>0:
                xsol=ga.solutions[-1]
                xt=dict(list(ga.xi.items())+list(xsol.items()))
                solutions=MeshAssembly(**xt)
                valid_cd=True
                for engr_num,cd in enumerate(self.center_distance):
                    if (solutions.center_distance[engr_num])<(cd[0]):
                        valid_cd=False
                    elif (solutions.center_distance[engr_num])>(cd[1]):
                        valid_cd=False
                if valid_cd:
                    liste_solutions.append(solutions)
                    compt_nb_sol+=1
                    if verbose:
                        print('valid solution n°{}'.format(compt_nb_sol))
                        print('Module: {}'.format(solutions.meshes[0][list(solutions.meshes[0].keys())[0]].rack.module))
                    if compt_nb_sol==nb_sol:
                        break
                else:
                    if verbose: 
                        print('unvalid solution')
            else:
                if verbose:
                    print('Solution non convergée')
        # TODO: commentaire en anglais
        #Tri des solutions convergées en fonction de la fonctionnelle actualisée
        list_fonctionnel_nv=[]
        for solutions in liste_solutions:
            fonctionnel=0
            for engr_num,cd in enumerate(self.center_distance):
                fonctionnel+=(solutions.center_distance[engr_num]-cd[0])**2
            list_fonctionnel_nv.append(fonctionnel)
        for indice_plex,num_plex in enumerate(npy.argsort(npy.array(list_fonctionnel_nv))):
            self.solutions.append(liste_solutions[num_plex])
            if verbose and indice_plex==0:
                print('Mesh sections: {}'.format(self.solutions[-1].gear_width))
                print('Center distances: {}'.format(self.solutions[-1].center_distance))
                print('Module: {}'.format(self.solutions[-1].meshes[0][list(self.solutions[-1].meshes[0].keys())[0]].rack.module))

