import numpy as npy
import math as mt
from scipy import interpolate
import os
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm
from scipy.optimize import minimize,fsolve
from scipy.interpolate import splprep, splev
import cma
#from sympy import *
import itertools
from jinja2 import Environment, PackageLoader, select_autoescape

import mechanical_components.LibSvgD3 as LibSvg

import persistent
#from dessia_common import ResultsDBClient
#import pyDOE

class Material(persistent.Persistent):
    
    def __init__(self):
        
        data_array,dico_axis,dico_nom,type_axis=self.DataCoeffYBIso()
        self.data_coeff_YB_Iso=self.GenereCoeff(data_array,dico_axis,dico_nom,type_axis)
    
        data_array,dico_axis,dico_nom,type_axis=self.DataWholerCurve()
        self.data_wholer_curve=self.GenereCoeff(data_array,dico_axis,dico_nom,type_axis)
        
        data_array,dico_axis,dico_nom,type_axis=self.DataGearMaterial()
        self.data_gear_material=self.GenereCoeff(data_array,dico_axis,dico_nom,type_axis)
        
    def GenereCoeff(self,data_array,dico_axis,dico_nom,type_axis):
        structure={}
        vect_x=npy.linspace(data_array[dico_nom['x']][0,0],data_array[dico_nom['x']][-1,0],len(data_array[dico_nom['x']][:,0]))
        vect_y=npy.linspace(data_array[dico_nom['y']][0,1],data_array[dico_nom['y']][-1,1],len(data_array[dico_nom['y']][:,1]))
        
        for i,j in type_axis.items():
            if j=='Log':
                axe_reel=[mt.log10(dico_axis[i][0]),mt.log10(dico_axis[i][-1])]
            if j=='Linear':
                axe_reel=dico_axis[i]
            if i=='x':
                ax,bx=self.AxisLinear(axe_reel,vect_x)
            if i=='y':
                ay,by=self.AxisLinear(axe_reel,vect_y)

        for key,data in dico_nom.items():
            if not key in ['x','y']:
                export=[]
                for item in data_array[data]:
                    data_x=item[0]*ax+bx
                    data_y=item[1]*ay+by
                    export.append([data_x,data_y])
                export=npy.array(export)
                structure[key]=export
        structure['x']=type_axis['x']
        structure['y']=type_axis['y']
        return structure
    
    def DataCoeffYBIso(self):
        data_svg=['56.097561,477.97196 81.707319,-1.21951 81.70732,1.21951 82.92682,0 79.2683,-1.21951 82.92683,0 81.70731,-1.21951 81.70732,0 81.70732,1.21951',
                               '57.317073,476.75245 -1.219512,-82.92683 0,-82.92683 1.219512,-84.14634 0,-82.92683 1.219512,-82.92683',
                               '56.097561,58.459766 76.829269,59.756094 71.95122,43.90244 89.02439,45.12195 78.04878,32.92683 75.60976,19.5122 74.39024,7.31707 90.2439,1.21951 97.56098,1.21952']
        dico_nom={'evol_coeff_yb_iso':2,'x':0,'y':1}
        dico_axis={'x':[0,5,10,15,20,25,30,35,40],'y':[0.5,0.6,0.7,0.8,0.9,1]}
        type_axis={'x':'Linear','y':'Linear'}
        data_array=self.SVG2Array(data_svg)
        return data_array,dico_axis,dico_nom,type_axis
    
    def DataWholerCurve(self):
        data_svg=['199.46895,423.19149 0,-42.67247 -0.99238,-41.68007 0.99238,-31.75625 -0.99238,-28.77911 0,-20.84004 0,-19.84765 0,-20.84004 0,-18.85528 0.99238,-17.86289 0,-14.88574 0,-16.87051 0,-12.90098 0,-9.92382 0,-11.9086 0,-11.90859 0,-10.916213','200.46133,423.19149 34.7334,-0.99239 20.84004,0.99239 14.88574,-0.99239 11.9086,0 36.29445,0.74708 36.12719,-0.80283 21.67632,0.80283 14.45088,0 11.23957,0 38.93709,-0.40141 36.1272,0.40141 21.67631,0 16.85936,-0.40141 11.23957,0 36.1272,-0.40142 37.33143,-0.40141 19.26784,0.40141 14.04947,0.80283 9.2325,0 36.93002,-0.40141','236,100.3622 63.42857,18.28572 61.14286,15.42857 93.71428,21.71429 57.71429,3.42857 65.71429,1.14285 65.71428,0.57143 38.28572,0','213.04348,126.27525 55.65217,34.78261 38.26087,22.60869 66.08696,6.95652 115.65217,-0.86956 103.47826,0 86.95652,-0.86957','221.14286,140.3622 60.57143,46.28572 34.85714,26.85714 43.42857,27.42857 13.14286,5.71429 78.28571,-0.57143 125.71429,1.14286 103.42857,-2.28572','234.28571,204.3622 50.28572,29.14286 37.71428,21.14286 41.14286,21.14286 31.42857,12 48,5.14285 65.71429,0 83.42857,-1.14285 88,0','237.3913,261.05786 112.17392,50.43478 83.47826,31.30435 69.56522,18.26087 84.34782,-0.86957 92.17391,-0.86956']
        dico_nom={'hardened_alloy_steel':2,'nitrided_alloy_steel':3,'through_hardened_steel':4,'surface_hardened_steel':5,'x':1,'y':0,'carbon_steel':6,'cast_iron':6,'bronze':6,'grey_iron':6}
        dico_axis={'x':[10000,100000000],'y':[20,100]}
        type_axis={'x':'Log','y':'Log'}
        data_array=self.SVG2Array(data_svg)
        return data_array,dico_axis,dico_nom,type_axis
    
    def DataGearMaterial(self):
        data_svg=['76.2931,579.71715 91.17956,0.9304 59.54583,0.9304 47.45059,-0.9304 134.90853,0 110.71804,-0.9304','76.2931,579.71715 0,-148.86459 -0.930404,-90.24915 0,-240.974553','350.76218,194.53003 53.96341,-28.84251 56.75463,-29.77292 40.93776,-22.32969','352.62299,209.41649 49.31139,-26.98171 39.07695,-19.53847 54.89382,-29.77292','319.12846,281.05757 71.64108,-46.52018 56.75462,-39.07696 66.98906,-45.58977','345.17976,370.37632 51.1722,-42.79857 53.96341,-42.79856','267.02585,377.81955 43.72897,-36.28574 38.14655,-33.49453 23.26009,-20.46888','163.75104,467.1383 37.21615,-30.70332 30.70332,-25.1209 16.74727,-17.67767','82.805926,535.98817 38.146554,-31.63372 20.46888,-15.81686 13.02565,-14.88646','103.27481,574.13472 37.21614,-33.49453 33.49453,-35.35534 25.1209,-26.0513 14.88646,-17.67767']
        dico_nom={'hardened_alloy_steel':2,'nitrided_alloy_steel':3,'through_hardened_steel':4,'surface_hardened_steel':5,'carbon_steel':6,'cast_iron':7,'bronze':8,'grey_iron':9,'x':0,'y':1}
        dico_axis={'x':[20,150],'y':[5,45]}
        type_axis={'x':'Log','y':'Log'}
        data_array=self.SVG2Array(data_svg)
        return data_array,dico_axis,dico_nom,type_axis
    
    def FunCoeff(self,x,data,type_x='Linear',type_y='Linear'):
        if type_x=='Log': 
            x=npy.log10(x)
        f = interpolate.interp1d(list(data[:,0]),list(data[:,1]), fill_value='extrapolate')
        sol=float(f(x))
        if type_y=='Log':
            sol=10**sol
        return sol
    
    def AxisLinear(self,axe,vect):
        a=(axe[0]-axe[-1])/(vect[0]-vect[-1])
        b=axe[0]-a*vect[0]
        return a,b
    
    def SVG2Array(self,data):
        # en entrée une liste de data SVG positionnee en relatif et 
        # en sortie liste Array avec positionnement en absolu
        export={}
        for i,dat in enumerate(data):
            data_temp=[]
            sol=dat.split(' ')
            for item in sol:
                temp=item.split(',')
                data_temp.append([float(temp[0]),float(temp[1])])
            export_temp=[data_temp[0]]
            for item in data_temp[1::]:
                export_temp.append([item[0]+export_temp[-1][0],item[1]+export_temp[-1][1]])
            export[i]=npy.array(export_temp)
        return export

class Rack(persistent.Persistent):
    def __init__(self,transverse_radial_pitch,transverse_pressure_angle_T=None,circular_tooth_thickness=None,gear_addendum=None,
                 gear_dedendum=None,root_radius_T=None,root_radius_R=None,transverse_pressure_angle_R=None):
        
        self.RackParam(transverse_radial_pitch,transverse_pressure_angle_T,circular_tooth_thickness,gear_addendum,gear_dedendum,
                       root_radius_T,root_radius_R,transverse_pressure_angle_R)
        
    def RackParam(self,transverse_radial_pitch,transverse_pressure_angle_T,circular_tooth_thickness,gear_addendum,gear_dedendum,
                  root_radius_T,root_radius_R,transverse_pressure_angle_R):
    
        self.transverse_pressure_angle_T=transverse_pressure_angle_T
        self.transverse_pressure_angle_R=transverse_pressure_angle_R
        self.transverse_radial_pitch=transverse_radial_pitch
        self.circular_tooth_thickness=circular_tooth_thickness
        self.gear_addendum=gear_addendum
        self.gear_dedendum=gear_dedendum
        self.root_radius_T=root_radius_T
        self.root_radius_R=root_radius_R
        self.transverse_pressure_angle_R=transverse_pressure_angle_R
        self.module=self.transverse_radial_pitch/npy.pi
        
        #if nothing define, we put the ISO definition for the rack
        if gear_addendum==None:
            self.gear_addendum=self.module
        if gear_dedendum==None:
            self.gear_dedendum=1.25*self.module
        if root_radius_T==None:
            self.root_radius_T=0.38*self.module
        if circular_tooth_thickness==None:
            self.circular_tooth_thickness=transverse_radial_pitch/2
        if transverse_pressure_angle_T==None:
            self.transverse_pressure_angle_T=20/180*npy.pi
        if transverse_pressure_angle_R==None:
            self.transverse_pressure_angle_R=self.transverse_pressure_angle_T
        if root_radius_R==None:
            self.root_radius_R=self.root_radius_T
            
        self.tooth_space=self.transverse_radial_pitch-self.circular_tooth_thickness
        self.whole_depth=self.gear_addendum+self.gear_dedendum
        self.clearance_T=self.root_radius_T-self.root_radius_T*npy.sin(self.transverse_pressure_angle_T)
        self.clearance_R=self.root_radius_R-self.root_radius_R*npy.sin(self.transverse_pressure_angle_R)
        self.clearance=npy.max([self.clearance_T,self.clearance_R])
        
        #paramètre pour la trochoide
        self.a=self.tooth_space/2-self.gear_dedendum*npy.tan(self.transverse_pressure_angle_T)-self.root_radius_T*npy.tan(1/2*npy.arctan(npy.cos(self.transverse_pressure_angle_T)/(npy.sin(self.transverse_pressure_angle_T))))
        self.b=self.gear_dedendum-self.root_radius_T
        
    def Update(self,transverse_radial_pitch,transverse_pressure_angle_T,circular_tooth_thickness=None,gear_addendum=None,
               gear_dedendum=None,root_radius_T=None,root_radius_R=None,transverse_pressure_angle_R=None):
        
        self.RackParam(transverse_radial_pitch,transverse_pressure_angle_T,circular_tooth_thickness,gear_addendum,gear_dedendum,
                  root_radius_T,root_radius_R,transverse_pressure_angle_R)
        
    def RackContours(self,list_number):
        p1=vm.Point2D((0,0))
        p2=p1.Translation((self.gear_addendum*npy.tan(self.transverse_pressure_angle_R),self.gear_addendum))
        p4=p1.Translation((self.circular_tooth_thickness,0))
        p3=p4.Translation((-self.gear_addendum*npy.tan(self.transverse_pressure_angle_T),self.gear_addendum))
        p5=p4.Translation((self.gear_dedendum*npy.tan(self.transverse_pressure_angle_T),-self.gear_dedendum))
        p7=p4.Translation((self.tooth_space,0))
        p6=p7.Translation((-self.gear_dedendum*npy.tan(self.transverse_pressure_angle_R),-self.gear_dedendum))
        L=primitives2D.RoundedLines2D([p1,p2,p3,p4,p5,p6,p7],{4:self.root_radius_T,5:self.root_radius_R},False)
        Rack_Elem=[]
        for i in list_number:
            Rack_Elem.append(L.Translation(((i)*(p7.vector-p1.vector))))
        return Rack_Elem
    
    def RackComplete(self,list_number):
        
        L=self.RackContours(list_number)
        p1=L[0].points[0]
        p6=L[-1].points[-1]
        p2=p1.Translation((-self.circular_tooth_thickness,0))
        p3=p2.Translation((0,2*self.whole_depth))
        p5=p6.Translation((self.circular_tooth_thickness,0))
        p4=p5.Translation((0,2*self.whole_depth))
        Lbis=primitives2D.RoundedLines2D([p1,p2,p3,p4,p5,p6],{},False)
        L.append(Lbis)
        return L
        
    def CriteriaEq(self):
        crit=[0]
        return crit
        
    def CriteriaIneq(self):
        crit=[self.transverse_radial_pitch-self.circular_tooth_thickness-self.gear_dedendum*npy.tan(self.transverse_pressure_angle_T)
        -self.gear_dedendum*npy.tan(self.transverse_pressure_angle_R)-(self.root_radius_T*npy.cos(self.transverse_pressure_angle_T)
        -npy.tan(self.transverse_pressure_angle_T)*self.root_radius_T*(1-npy.sin(self.transverse_pressure_angle_T)))
        -(self.root_radius_R*npy.cos(self.transverse_pressure_angle_R)-npy.tan(self.transverse_pressure_angle_R)
        *self.root_radius_R*(1-npy.sin(self.transverse_pressure_angle_R)))]
        crit.append(self.circular_tooth_thickness-(self.gear_addendum*npy.tan(self.transverse_pressure_angle_T)+self.gear_addendum*npy.tan(self.transverse_pressure_angle_R)))
        return crit
        
    def Mass(self):
        pass
    
    def VolumeModel(self):
        pass
    
    def SVGExport(self,number,name):
        
        R1=self.RackComplete(npy.arange(number))
        Ref=[vm.Line2D(vm.Point2D((-self.transverse_radial_pitch,0)),vm.Point2D(((number+1)*self.transverse_radial_pitch,0)))]
        Ref.append(vm.Line2D(vm.Point2D((-self.transverse_radial_pitch,self.gear_addendum)),vm.Point2D(((number+1)*self.transverse_radial_pitch,self.gear_addendum))))
        Ref.append(vm.Line2D(vm.Point2D((-self.transverse_radial_pitch,-self.gear_dedendum)),vm.Point2D(((number+1)*self.transverse_radial_pitch,-self.gear_dedendum))))
        SVG1=LibSvg.SVGTrace(700)
        SVG1.Convert(R1,'R1','black',0.0002,0)
        SVG1.Convert(Ref,'Ref','black',0.0001,1,'0.001px, 0.008px')
        SVG1.Show(name)
        
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
    
class Gear(persistent.Persistent):
    def __init__(self,tooth_number,rack_data,coefficient_profile_shift,gear_width,material,rim,alpha_rim,boring_diameter,thickness_rim):
        
        self.rack=rack_data
        
        save=[self.rack.transverse_radial_pitch,tooth_number,coefficient_profile_shift]
        self.save=save[:]
        self.GearParam(tooth_number,rack_data,coefficient_profile_shift,gear_width,material,rim,alpha_rim,boring_diameter,thickness_rim)
        self._RootDiameterActive()
        
    ### Data geometry
        
    def GearParam(self,tooth_number,rack_data,coefficient_profile_shift,gear_width,material,rim,alpha_rim,boring_diameter,thickness_rim):
        
        self.material=material
        self.gear_width=gear_width
        self.tooth_number=tooth_number
        self.coefficient_profile_shift=coefficient_profile_shift
        self.pitch_diameter_factory=rack_data.transverse_radial_pitch*tooth_number/npy.pi
        
        self.outside_diameter=self.pitch_diameter_factory+2*(rack_data.gear_addendum+rack_data.module*self.coefficient_profile_shift)
        self.root_diameter=self.pitch_diameter_factory-2*(rack_data.gear_dedendum-rack_data.module*self.coefficient_profile_shift)
        
        #data sur le cercle primitif ()
        self.circular_tooth_thickness=rack_data.circular_tooth_thickness+rack_data.module*self.coefficient_profile_shift*npy.tan(rack_data.transverse_pressure_angle_T)+rack_data.module*self.coefficient_profile_shift*npy.tan(rack_data.transverse_pressure_angle_R)
        self.tooth_space=rack_data.transverse_radial_pitch-self.circular_tooth_thickness
        
        self.transverse_pressure_angle_T_factory=rack_data.transverse_pressure_angle_T
        circular_tooth_thickness_angle=self.circular_tooth_thickness/(self.pitch_diameter_factory/2)
        tooth_space_angle=self.tooth_space/(self.pitch_diameter_factory/2)
        self.base_diameter=self.pitch_diameter_factory*npy.cos(self.transverse_pressure_angle_T_factory)
        temp=self.pitch_diameter_factory-2*rack_data.gear_dedendum+2*rack_data.root_radius_T-2*rack_data.root_radius_T*npy.sin(rack_data.transverse_pressure_angle_T)-2*rack_data.module*self.coefficient_profile_shift
        
        #prise de la valeur absolue pour gérer les cas de changement de signe lié au déport
        theta2=abs(fsolve((lambda theta:self.base_diameter/2*npy.cos(rack_data.transverse_pressure_angle_T-theta)/npy.cos(theta)-temp/2) ,0)[0])
        
        #print(theta2)
        self.theta0=npy.tan(theta2)-rack_data.transverse_pressure_angle_T
        self.alpha_root_diameter_active=theta2
        self.tan_alpha_root_diameter_active=npy.tan(self.alpha_root_diameter_active)
        self.alpha_outside_diameter=npy.arccos(self.base_diameter/self.outside_diameter)
        self.alpha_pitch_diameter=npy.arccos(self.base_diameter/self.pitch_diameter_factory)
        self.inv_alpha_pitch_diameter=npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter
        self.base_circular_tooth_thickness=self.base_diameter/2*(2*self.circular_tooth_thickness/self.pitch_diameter_factory+2*self.inv_alpha_pitch_diameter)
        self.transverse_base_pitch=rack_data.transverse_radial_pitch*npy.cos(rack_data.transverse_pressure_angle_T)
        
        self.root_active_angle=self.tooth_space/(self.pitch_diameter_factory/2)-2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter-npy.tan(self.alpha_root_diameter_active)+self.alpha_root_diameter_active)
        self.root_angle=self.tooth_space/(self.pitch_diameter_factory/2)-2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
        self.root_gear_angle=self.circular_tooth_thickness/(self.pitch_diameter_factory/2)+2*(npy.tan(self.alpha_pitch_diameter)-self.alpha_pitch_diameter)
        
        self.outside_active_angle=2*self.circular_tooth_thickness/self.pitch_diameter_factory-2*(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter-npy.tan(self.alpha_pitch_diameter)+self.alpha_pitch_diameter)
        
        
#        self.root_diameter_active=self.base_diameter/npy.cos(theta2)
#        self.alt_internal_contact=self.root_diameter_active/2*npy.sin(rack_data.transverse_pressure_angle_T-theta2)
#        self._RootDiameterActive()
        
        save=[self.rack.transverse_radial_pitch,self.tooth_number,
              self.coefficient_profile_shift]
        if not self.save==[self.rack.transverse_radial_pitch,self.tooth_number,
                           self.coefficient_profile_shift]:
            self._RootDiameterActive()
#            print(save,self.root_diameter_active)
            
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.pitch_diameter_factory/2
        rho=self.rack.root_radius_T
        self.root_diameter_active_bis=2*(rho+b**2/(r+b))
        
        #initialisation rim
        self.rim=rim
        self.alpha_rim=alpha_rim
        self.boring_diameter=boring_diameter
        self.thickness_rim=thickness_rim
        self.RimDef()
        
        self.save=save[:]
        
        
    def Update(self,tooth_number,rack_data,coefficient_profile_shift,gear_width,material,rim,alpha_rim,boring_diameter,thickness_rim):
        
        self.GearParam(tooth_number,rack_data,coefficient_profile_shift,gear_width,material,rim,alpha_rim,boring_diameter,thickness_rim)
        self.RimDef()
        save=[self.rack.transverse_radial_pitch,self.tooth_number,self.coefficient_profile_shift]
        if not self.save==[self.rack.transverse_radial_pitch,self.tooth_number,self.coefficient_profile_shift]:
            self._RootDiameterActive()
#            print(save,self.root_diameter_active)
        self.save=save[:]
        
    def RimDef(self):
        if str(self.rim) in 'rim_gear':
            self.boring_diameter_out=self.boring_diameter+4e-3
            self.rim_diameter_int=self.root_diameter-self.alpha_rim*(self.outside_diameter-self.root_diameter)
            rim_diam_max=(self.rim_diameter_int-self.boring_diameter_out)/2
            rim_diam_thickness=(self.gear_width-self.thickness_rim)
            self.rim_diam=0.8*min(rim_diam_max,rim_diam_thickness)
            self.rim_diam_x=self.thickness_rim/2+self.rim_diam/2
            self.rim_diam_r1=self.boring_diameter_out/2+self.rim_diam/2
            self.rim_diam_r2=self.rim_diameter_int/2-self.rim_diam/2
    
    def RimContour(self):
        hg=self.alpha_rim*(self.outside_diameter-self.root_diameter)/2
        if str(self.rim) in 'rim_gear':
            p=[vm.Point2D((0,self.boring_diameter/2))]
            p.append(vm.Point2D((-self.gear_width/2,self.boring_diameter/2)))
            p.append(vm.Point2D((-self.gear_width/2,self.boring_diameter_out/2)))
            p.append(vm.Point2D((-self.thickness_rim/2,self.boring_diameter_out/2)))
            p.append(vm.Point2D((-self.thickness_rim/2,self.rim_diameter_int/2)))
            p.append(vm.Point2D((-self.gear_width/2,self.rim_diameter_int/2)))
            p.append(vm.Point2D((-self.gear_width/2,self.root_diameter/2-hg/2)))
            p.append(vm.Point2D((self.gear_width/2,self.root_diameter/2-hg/2)))
            p.append(vm.Point2D((self.gear_width/2,self.rim_diameter_int/2)))
            p.append(vm.Point2D((self.thickness_rim/2,self.rim_diameter_int/2)))
            p.append(vm.Point2D((self.thickness_rim/2,self.boring_diameter_out/2)))
            p.append(vm.Point2D((self.gear_width/2,self.boring_diameter_out/2)))
            p.append(vm.Point2D((self.gear_width/2,self.boring_diameter/2)))
            p.append(p[0])
            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{3:self.rim_diam/2,4:self.rim_diam/2,9:self.rim_diam/2,10:self.rim_diam/2},False).primitives)
        elif str(self.rim) in 'shaft_gear':
            p=[vm.Point2D((0,0))]
            p.append(vm.Point2D((-self.gear_width/2,0)))
            p.append(vm.Point2D((-self.gear_width/2,self.root_diameter/2-hg/2)))
            p.append(vm.Point2D((self.gear_width/2,self.root_diameter/2-hg/2)))
            p.append(vm.Point2D((self.gear_width/2,0)))
            p.append(vm.Point2D((0,0)))
            ref=vm.Contour2D(primitives2D.RoundedLines2D(p,{},False).primitives)
#        print(str(self.rim) in 'shaft_gear',str(self.rim) in 'rim_gear',self.rim)
        return ref
        
    
    def GearSection(self,diameter):
        #epaisseur de la dent au diameter
        alpha_diameter=npy.arccos(self.base_diameter/diameter)
        theta1=(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter)-(npy.tan(alpha_diameter)-alpha_diameter)
        return diameter/2*(2*theta1+self.outside_active_angle)
    
    def GearSectionISO(self,angle):
        
        drap=1
        a=drap*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.pitch_diameter_factory/2
        rho=self.rack.root_radius_T
        
        theta0=fsolve((lambda theta:a + b*npy.tan(theta) + r*(-angle - self.root_angle/2 - theta + npy.pi/2)) ,0)[0]
        phi0=(a-b*npy.tan(theta0))/r
        pt_iso=self._Trochoide(phi0)
        angle0=npy.arctan(pt_iso[1]/pt_iso[0])-self.root_angle/2
        self.angle_iso=self.root_gear_angle-2*angle0
        self.diameter_iso=2*norm(pt_iso)
        self.s_thickness_iso=self.diameter_iso*npy.sin(self.angle_iso/2)
        self.h_height_iso=(self.s_thickness_iso/2)/npy.tan(angle)
        self.angle_iso=angle
        
    def _RootDiameterActive(self):
        
        #Analyse diam pied de dent actif
        drap=1
        a=drap*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.pitch_diameter_factory/2
        rho=self.rack.root_radius_T
        phi=-(a+b*npy.tan(npy.pi/2-self.transverse_pressure_angle_T_factory))/r
        data=2*norm(self._Trochoide(phi))
        self.root_diameter_active=data
        self.alpha_root_diameter_active=npy.arccos(self.base_diameter/self.root_diameter_active)
        self.tan_alpha_root_diameter_active=npy.tan(self.alpha_root_diameter_active)
        self.phi_trochoide=phi
        #corde de la dent au diam de pied actif
        self.root_gear_length=npy.sin((self.root_gear_angle-2*(npy.tan(self.alpha_root_diameter_active)-self.alpha_root_diameter_active))/2)*(self.root_diameter_active/2)*2
        
    ### Contrainte
        
    def CriteriaIneq(self):
        crit=[self.root_diameter_active-self.base_diameter-0.3*self.transverse_base_pitch]
        crit.append(self.outside_diameter-self.root_diameter-0.3*self.transverse_base_pitch)
        #Distance top of the rack
        #DistRI=self.DistRackInvolute()
        #crit.append(DistRI)
        #Top of the gear must be positive
        crit.append(self.base_circular_tooth_thickness*2/self.base_diameter-2*(npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter))
        return crit
        
    def CriteriaEq(self):
        crit=[0]
        return crit
    
    ### Trace
        
    def GearContours(self,discret=10,list_number=[None]):
        #Analytical tooth profil
        self._RootDiameterActive()
        if list_number==[None]:
            list_number=npy.arange(int(self.tooth_number))
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
                               self.tan_alpha_root_diameter_active,discret)
        else:
            drap=-1
            theta=npy.linspace(self.tan_alpha_root_diameter_active,
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
            L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
            self.rac=L.points[-1]
        else:
            L=ref.Rotation(vm.Point2D((0,0)),
                           self.base_circular_tooth_thickness*2/self.base_diameter)
            L=L.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
            L.points[0]=self.rac
        return L
        
    def _Involute(self,tan_alpha):
        
        sol=(self.base_diameter/2*npy.cos(tan_alpha)+self.base_diameter/2*tan_alpha*npy.sin(tan_alpha),
             self.base_diameter/2*npy.sin(tan_alpha)-self.base_diameter/2*tan_alpha*npy.cos(tan_alpha))
        return sol
    
    def _TrochoideTrace(self,discret,number,ind='T'):
        if ind=='T':
            drap=1
        else:
            drap=-1
        
        a=drap*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.pitch_diameter_factory/2
        rho=self.rack.root_radius_T
        self.phi0=npy.arctan((a)/(self.pitch_diameter_factory-b))
        self.phi0=a/(self.pitch_diameter_factory/2)
        
        ref=[]
        self._RootDiameterActive()
        
        if ind=='R':
            theta=npy.linspace(self.phi0,drap*self.phi_trochoide,discret)
        else:
            theta=npy.linspace(drap*self.phi_trochoide,self.phi0,discret)
        for t in theta:
            ref.append(vm.Point2D((self._Trochoide(t,ind))))
        ref=primitives2D.RoundedLines2D(ref,{},False)
        ref=ref.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        if ind=='T':
            L1=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
            L1.points[0]=self.rac
        else:
#            L=ref.Rotation(vm.Point2D((0,0)),self.base_circular_tooth_thickness*2/self.base_diameter)
            L1=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
            self.rac=L1.points[-1]

        return L1
    
    def _Trochoide(self,phi,ind='T'):
        
        if ind=='T':
            drap=1
        else:
            drap=-1
        a=drap*self.rack.a
        b=self.rack.b-self.rack.module*self.coefficient_profile_shift
        r=self.pitch_diameter_factory/2
        rho=self.rack.root_radius_T
        
        x2=rho*npy.sin(npy.arctan((a-r*phi)/b)-phi)+a*npy.cos(phi)-b*npy.sin(phi)+r*(npy.sin(phi)-phi*npy.cos(phi))
        y2=-rho*npy.cos(npy.arctan((a-r*phi)/b)-phi)-a*npy.sin(phi)-b*npy.cos(phi)+r*(npy.cos(phi)+phi*npy.sin(phi))
        sol=(y2,x2)
        return sol
    
    def _RootCircleTrace(self,number):
        
        theta4=-self.root_angle/2
        
        drap=1
        a=drap*self.rack.a
        self.phi0=a/(self.pitch_diameter_factory/2)
        p1=vm.Point2D((self._Trochoide(self.phi0,'T')))
        p1=p1.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        drap=-1
        a=drap*self.rack.a
        self.phi0=a/(self.pitch_diameter_factory/2)
        p2=vm.Point2D((self._Trochoide(self.phi0,'R')))
        p2=p2.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        
        ref=primitives2D.RoundedLines2D([p1,p2],{},False)
        L2=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        return L2
    
    def _OutsideTrace(self,number):
        #trace du sommet des dents en arc de cercle
        theta4=npy.tan(self.alpha_outside_diameter)-self.alpha_outside_diameter
        p1=vm.Point2D((self.outside_diameter/2*npy.cos(theta4),self.outside_diameter/2*npy.sin(theta4)))
        p2=p1.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        p3=p2.Rotation(vm.Point2D((0,0)),self.outside_active_angle/2)
        #ref=vm.Arc2D(p1,p2,p3)
        ref=primitives2D.RoundedLines2D([p3,p2,p1],{},False)
        L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.tooth_number)
        return L
    
    def _TrochoideSecondary(self,list_number,ind,angle,number=10):
        
        ref=[]
        for t in npy.linspace(-angle,angle,number):
            ref.append(vm.Point2D((self._Trochoide(t,ind))))
        ref=primitives2D.RoundedLines2D(ref,{},False)
        ref=ref.Rotation(vm.Point2D((0,0)),-self.root_angle/2)
        L=[]
        for i in list_number:
            L.append(ref.Rotation(vm.Point2D((0,0)),-i*2*npy.pi/self.tooth_number))
        return L
        
    def _PosRack(self,alpha,list_number=[1]):
        
        angle=npy.tan(alpha)-self.rack.transverse_pressure_angle_T
        pt0=vm.Point2D((self.pitch_diameter_factory/2,0))
        repere=[pt0,pt0.Translation((1,0)),pt0.Translation((0,1))]
        for i in range(3):
            repere[i]=repere[i].Rotation(vm.Point2D((0,0)),angle)
        ptT=vm.Point2D(self._Involute(npy.tan(alpha)))
        DistY=npy.dot(ptT.vector-repere[0].vector,repere[2].vector-repere[0].vector)
        DistX=npy.dot(ptT.vector-repere[0].vector,repere[1].vector-repere[0].vector)
        decal=DistY+self.rack.circular_tooth_thickness-npy.tan(self.rack.transverse_pressure_angle_T)*DistX+npy.tan(self.rack.transverse_pressure_angle_T)*(self.rack.module*self.coefficient_profile_shift)
        temp=self.rack.RackComplete(list_number)
        for (i,j) in enumerate(temp):
            j=j.Rotation(vm.Point2D((0,0)),-npy.pi/2)
            j=j.Translation((self.pitch_diameter_factory/2+self.rack.module*self.coefficient_profile_shift,decal))
            j=j.Rotation(vm.Point2D((0,0)),angle)
            temp[i]=j
        rack_profil=temp
        return repere,rack_profil
    
    ### Stress
    
    def _CoeffYSIso(self):
        #facteur de concentration de contrainte en pied de dent
        rho_f=self.rack.root_radius_T+self.rack.b**2/(self.pitch_diameter_factory/2+self.rack.b)
        self.coeff_ys_iso=1+0.15*self.s_thickness_iso/rho_f
    
    ### Generique

    def _BSpline(self,pt):
        x=[]
        y=[]
        for i in pt:
            x.append(i.vector[0])
            y.append(i.vector[1])
        tck, u=splprep([x, y], s=0)
        pas=0.01
        NewPT=splev(npy.arange(0,1+pas,pas), tck)
        pt=[]
        for i,j in zip(NewPT[0],NewPT[1]):
            pt.append(vm.Point2D((i,j)))
        return tck,pt
    
    def _Trace(self,x,y):
        p=[vm.Point2D((x[0],y[0]))]
        for i in range(1,len(x)):
            p.append(vm.Point2D((x[i],y[i])))
        ref=primitives2D.RoundedLines2D(p,{},False)
        return ref
        
    def Mass(self):
        rho=7800
        masse=(self.root_diameter/2)**2*npy.pi*self.gear_width*rho
        masse+=npy.sin((self.root_gear_angle-self.outside_active_angle)/2)*self.root_diameter/2*(self.outside_diameter/2-self.root_diameter/2)*self.gear_width*rho
        masse+=((self.outside_diameter/2)**2*npy.pi/(2*npy.pi)*self.outside_active_angle-(self.root_diameter/2)**2*npy.pi/(2*npy.pi)*self.outside_active_angle)*self.gear_width*rho
        return masse
    
    def VolumeModel(self):
        p=vm.Point3D((0,0,0))
        x=vm.Vector3D((1,0,0))
        y=vm.Vector3D((0,1,0))
        z=vm.Vector3D((0,0,self.width))
        wheel3D=vm.primitives3D.ExtrudedProfile(p,x,y,self.WheelContours(),z)
        lever3D=vm.primitives3D.ExtrudedProfile(p,x,y,self.WheelContours(),z)
        return vm.VolumeModel([wheel3D,lever3D])
    
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

        d['mass']=self.Mass()
        
        try:
            del d['rac']
        except KeyError:
            pass
        
        del d['save']
        d['rack']=self.rack.Dict()
        return d
    
    ### Export
    
    def GearGenerationSVGExport(self,name):

        L1=self.GearContours(20,[int(self.tooth_number)-1,0,1])
        L2=[]
        for i in npy.linspace(self.alpha_root_diameter_active,self.alpha_outside_diameter,5):
            L=self._PosRack(i,[-1,0,1])
            L2.extend(L[1])
        L3=self._TrochoideSecondary([0,self.tooth_number-1],'T',0.8,100)
        Temp=self._TrochoideSecondary([0,self.tooth_number-1],'R',0.8,100)
        L3.extend(Temp)
        #G1=vm.Contour2D(L1)
        #G1.MPLPlot()
        SVG1=LibSvg.SVGTrace(700)
        SVG1.Convert(L1,'Rack','black',0.005,0)
        SVG1.Convert(L2,'Rack','black',0.005,0,'0.01px, 0.08px')
        SVG1.Convert(L3,'Rack','blue',0.005,0)
        SVG1.Show(name)
    
    def CSVExport(self):
        d=self.__dict__.copy()
        (d1,d2)=self.rack.CSVExport()
        del d['rack']
        del d['save']
        return list(d.keys())+d1,list(d.values())+d2
        
class GearAssembly(persistent.Persistent):
    def __init__(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
             coefficient_profile_shift2,gear_width,maximum_torque,
             material1,material2,nb_cycle1,
             rim1,rim2,alpha_rim1,alpha_rim2,boring_diameter1,boring_diameter2,thickness_rim1,thickness_rim2,
             transverse_pressure_angle_rack_T1=None,
             transverse_pressure_angle_rack_T2=None,circular_tooth_thickness_rack1=None,
             circular_tooth_thickness_rack2=None,
             gear_addendum_rack1=None,gear_addendum_rack2=None,gear_dedendum_rack1=None,gear_dedendum_rack2=None,
             root_radius_T1=None,root_radius_T2=None,root_radius_R1=None,
             root_radius_R2=None,transverse_pressure_angle_rack_R1=None,transverse_pressure_angle_rack_R2=None):
        
        self.GearAssemblyParam(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
               coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
               transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
               gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
               root_radius_T1,root_radius_T2,root_radius_R1,
               root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)
        
        self.Rack1=Rack(self.transverse_radial_pitch_rack1,self.transverse_pressure_angle_rack_T1,circular_tooth_thickness_rack1,
                        gear_addendum_rack1,gear_dedendum_rack1,root_radius_T1,root_radius_R1,transverse_pressure_angle_rack_R1)
        self.Rack2=Rack(self.transverse_radial_pitch_rack2,self.transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack2,
                        gear_addendum_rack2,gear_dedendum_rack2,root_radius_T2,root_radius_R2,transverse_pressure_angle_rack_R2)
        
        self.Gear1=Gear(self.Z1,self.Rack1,self.coefficient_profile_shift1,self.gear_width,material1,rim1,alpha_rim1,boring_diameter1,thickness_rim1)
        self.Gear2=Gear(self.Z2,self.Rack2,self.coefficient_profile_shift2,self.gear_width,material2,rim2,alpha_rim2,boring_diameter2,thickness_rim2)
        
        self.GearAssemblyParam2(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
               coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
               transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
               gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
               root_radius_T1,root_radius_T2,root_radius_R1,
               root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)
        
        #chargement des abaques
        self.material1=material1
        self.material2=material2
        self.nb_cycle1=nb_cycle1
        self.material=Material()
        
        #initialisation contrainte ISO
        self.SigmaISO()
        self.SigmaMaterialISO()
        
        self.rim1=rim1
        self.alpha_rim1=alpha_rim1
        self.boring_diameter1=boring_diameter1
        self.thickness_rim1=thickness_rim1
        self.rim2=rim2
        self.alpha_rim2=alpha_rim2
        self.boring_diameter2=boring_diameter2
        self.thickness_rim2=thickness_rim2
         
    ### Data Geometry
        
    def GearAssemblyParam(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2):
        
        self.Z1=Z1
        self.Z2=Z2
        self.center_distance=center_distance
        self.DF1=2*self.center_distance*self.Z1/self.Z2/(1+self.Z1/self.Z2)
        self.DF2=2*self.center_distance-self.DF1
        self.gear_width=gear_width
        
        self.transverse_pressure_angle=transverse_pressure_angle
        self.transverse_radial_pitch=npy.pi*self.DF1/self.Z1
        self.helix_angle=helix_angle
        if npy.round(helix_angle,1)==0:
            self.normal_circular_pitch=self.transverse_radial_pitch
            self.axial_pitch=None
            self.normal_pressure_angle=self.transverse_pressure_angle
        else:
            self.normal_circular_pitch=self.transverse_radial_pitch*npy.cos(self.helix_angle)
            self.axial_pitch=self.normal_circular_pitch/npy.sin(helix_angle)
            self.normal_pressure_angle=npy.arctan(npy.tan(self.transverse_pressure_angle)*npy.cos(self.helix_angle))
            
        self.DB1=self.DF1*npy.cos(self.transverse_pressure_angle)
        self.DB2=self.DF2*npy.cos(self.transverse_pressure_angle)
        self.coefficient_profile_shift1=coefficient_profile_shift1
        self.coefficient_profile_shift2=coefficient_profile_shift2

        self.transverse_pressure_angle_rack_T1=transverse_pressure_angle_rack_T1
        self.transverse_pressure_angle_rack_T2=transverse_pressure_angle_rack_T2
        if transverse_pressure_angle_rack_T1==None:
            self.transverse_pressure_angle_rack_T1=20/180*npy.pi
        if transverse_pressure_angle_rack_T2==None:
            self.transverse_pressure_angle_rack_T2=20/180*npy.pi
            
        self.pitch_diameter_factory1=self.DB1/npy.cos(self.transverse_pressure_angle_rack_T1)
        self.pitch_diameter_factory2=self.DB2/npy.cos(self.transverse_pressure_angle_rack_T2)
        self.transverse_radial_pitch_rack1=npy.pi*self.pitch_diameter_factory1/self.Z1
        self.transverse_radial_pitch_rack2=npy.pi*self.pitch_diameter_factory2/self.Z2
        
        self.maximum_torque=maximum_torque
        self.normal_load=self.maximum_torque*2/(self.DB1)
        self.tangential_load=self.maximum_torque*2/(self.DF1)
        self.radial_load=npy.tan(self.transverse_pressure_angle)*self.tangential_load
        
    def GearAssemblyParam2(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
                       coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
                       transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
                       gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
                       root_radius_T1,root_radius_T2,root_radius_R1,
                       root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2):
        
        self.transverse_base_pitch=self.transverse_radial_pitch*npy.cos(self.transverse_pressure_angle)
        
        #self.transverse_radial_pitch_rack1=self.transverse_base_pitch/npy.cos(self.transverse_pressure_angle_rack_T1)
        #self.transverse_radial_pitch_rack2=self.transverse_base_pitch/npy.cos(self.transverse_pressure_angle_rack_T2)
        self.circular_tooth_thickness1=self.Rack1.circular_tooth_thickness*npy.cos(self.transverse_pressure_angle_rack_T1)/npy.cos(self.transverse_pressure_angle)
        self.circular_tooth_thickness1=self.Gear1.GearSection(self.DF1)
        self.circular_tooth_thickness2=self.Rack2.circular_tooth_thickness*npy.cos(self.transverse_pressure_angle_rack_T2)/npy.cos(self.transverse_pressure_angle)
        self.circular_tooth_thickness2=self.Gear2.GearSection(self.DF2)
        self.space_width1=self.transverse_radial_pitch-self.circular_tooth_thickness1
        self.space_width2=self.transverse_radial_pitch-self.circular_tooth_thickness2
        
        self.linear_backlash=self.space_width1-self.circular_tooth_thickness2
        self.linear_backlash2=self.space_width2-self.circular_tooth_thickness1
        
        #formule spur gear
        self.radial_contact_ratio=1/2*(npy.sqrt(self.Gear1.outside_diameter**2-self.DB1**2)+npy.sqrt(self.Gear2.outside_diameter**2-self.DB2**2)-2*self.center_distance*npy.sin(self.transverse_pressure_angle))/(self.transverse_radial_pitch*npy.cos(self.transverse_pressure_angle))
        self.axial_contact_ratio=self.gear_width*npy.sin(self.helix_angle)/self.normal_circular_pitch
        self.total_contact_ratio=self.radial_contact_ratio+self.axial_contact_ratio
        
        #calcul de la hauteur pour la contrainte de Lewis
        self.gear_height_lewis1=self.Gear1.outside_diameter/2-self.Gear1.root_diameter_active/2-(self.Gear1.outside_active_angle/2*self.Gear1.outside_diameter/2*npy.tan(self.transverse_pressure_angle))
        self.gear_height_lewis2=self.Gear2.outside_diameter/2-self.Gear2.root_diameter_active/2-(self.Gear2.outside_active_angle/2*self.Gear2.outside_diameter/2*npy.tan(self.transverse_pressure_angle))
        
    def Update(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
           coefficient_profile_shift2,gear_width,maximum_torque,
           material1,material2,nb_cycle1,
           rim1,rim2,alpha_rim1,alpha_rim2,boring_diameter1,boring_diameter2,thickness_rim1,thickness_rim2,
           transverse_pressure_angle_rack_T1,
           transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
           gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
           root_radius_T1,root_radius_T2,root_radius_R1,
           root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2):
        
        self.GearAssemblyParam(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
           coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
           transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
           gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
           root_radius_T1,root_radius_T2,root_radius_R1,
           root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)
        
        self.Rack1.Update(self.transverse_radial_pitch_rack1,transverse_pressure_angle_rack_T1,circular_tooth_thickness_rack1,
                        gear_addendum_rack1,gear_dedendum_rack1,root_radius_T1,root_radius_R1,transverse_pressure_angle_rack_R1)
        self.Rack2.Update(self.transverse_radial_pitch_rack2,transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack2,
                        gear_addendum_rack2,gear_dedendum_rack2,root_radius_T2,root_radius_R2,transverse_pressure_angle_rack_R2)
        self.Gear1.Update(self.Z1,self.Rack1,coefficient_profile_shift1,self.gear_width,self.material1,rim1,alpha_rim1,boring_diameter1,thickness_rim1)
        self.Gear2.Update(self.Z2,self.Rack2,coefficient_profile_shift2,self.gear_width,self.material2,rim2,alpha_rim2,boring_diameter2,thickness_rim2)
        
        self.GearAssemblyParam2(Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,
           coefficient_profile_shift2,gear_width,maximum_torque,transverse_pressure_angle_rack_T1,
           transverse_pressure_angle_rack_T2,circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
           gear_addendum_rack1,gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
           root_radius_T1,root_radius_T2,root_radius_R1,
           root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)
        
        self.SigmaISO()
        self.SigmaMaterialISO()
        
    ### Contraintes
        
    def CriteriaIneq(self):
#        print(self.radial_contact_ratio,self.transverse_radial_pitch,self.transverse_pressure_angle)
        crit=[
              (2*self.center_distance-self.Gear1.root_diameter_active-self.Gear2.outside_diameter),
              2*self.center_distance-self.Gear2.root_diameter_active-self.Gear1.outside_diameter,
              self.space_width1-self.circular_tooth_thickness2,
              self.space_width2-self.circular_tooth_thickness1,
              self.Gear1.outside_diameter-self.DF1,
              self.Gear2.outside_diameter-self.DF2,
              self.DF1-self.DB1,
              self.DF2-self.DB2,
              self.radial_contact_ratio-1,
              self.linear_backlash-0.05*10**(-3),
#              (0.5*10**(-3)-self.linear_backlash),
              (self.sigma_lim1-self.sigma_iso_1)/self.sigma_lim1,
              (self.sigma_lim2-self.sigma_iso_2)/self.sigma_lim2,
              ]
        return crit
        
    def CriteriaEq(self):
        
        crit=[
              self.Gear1.tooth_number-int(self.Gear1.tooth_number),
              self.Gear2.tooth_number-int(self.Gear2.tooth_number),
              ]
        return crit
    
    ### Trace
    
    def GearAssemblyTrace(self,list_gear,list_center,list_rot):
        
        TG=[]
        for (i,j,k) in zip(list_gear,list_center,list_rot):
            temp=[]
            for m in i:
                temp1=m.Translation(j)
                temp2=temp1.Rotation(vm.Point2D(j),k)
                temp.append(temp2)
            TG.append(temp)
        return TG
        
    def InitialPosition(self):
        
        fun = (lambda tan_alpha : (norm(self.Gear1._Involute(tan_alpha))-(self.center_distance-self.DF2/2))**2)
#        bnds = (0,1)
        sol=minimize(fun,[0.1], method='SLSQP', tol=1e-20)
        xsol=sol.x
        Angle1=xsol[0]-npy.arctan(xsol[0])
        fun = (lambda tan_alpha : (norm(self.Gear2._Involute(tan_alpha))-(self.DF2/2))**2)
#        bnds = (0,1)
        sol=minimize(fun,[0.1], method='SLSQP', tol=1e-20)
        xsol=sol.x
        Angle2=xsol[0]-npy.arctan(xsol[0])
        
        Angle1=npy.arccos(self.Gear1.base_diameter/self.DF1)
        Angle2=npy.arccos(self.Gear2.base_diameter/self.DF2)
        
        #Gear2Angle=-Angle2+npy.pi-Angle1*self.Z1/self.Z2
        Gear1Angle=-(npy.tan(Angle1)-Angle1)
        Gear2Angle=-(npy.tan(Angle2)-Angle2)+npy.pi
        #Gear1Angle=0
        
        return [Gear1Angle,Gear2Angle]
        
    ### Stress
    
    def SigmaLewis(self):
        
        self.sigma_lewis_maximum1=6*self.tangential_load*self.gear_height_lewis1/(self.gear_width*self.Gear1.root_gear_length**2)
        self.sigma_lewis_maximum2=6*self.tangential_load*self.gear_height_lewis2/(self.gear_width*self.Gear2.root_gear_length**2)
        
    def SigmaISO(self):
        
        self._CoeffYFIso()
        self._CoeffYEIso()
        self._CoeffYBIso()
        self.sigma_iso_1=abs(self.tangential_load/(self.gear_width*self.Gear1.rack.module)*self.coeff_yf_iso1*self.coeff_ye_iso*self.coeff_yb_iso)
        self.sigma_iso_2=abs(self.tangential_load/(self.gear_width*self.Gear2.rack.module)*self.coeff_yf_iso2*self.coeff_ye_iso*self.coeff_yb_iso)
        
    def SigmaMaterialISO(self):
        
        safety_factor=1
        data=self.material.data_wholer_curve
        sgl1=self.material.FunCoeff(self.nb_cycle1,data[self.material1],data['x'],data['y'])
        data=self.material.data_gear_material
        sgl1=self.material.FunCoeff(sgl1,data[self.material1],data['x'],data['y'])
        self.Gear1._CoeffYSIso()
        self.sigma_lim1=(sgl1/(safety_factor*self.Gear1.coeff_ys_iso))*10**7
        
        data=self.material.data_wholer_curve
        sgl2=self.material.FunCoeff(self.nb_cycle1*self.Z1/self.Z2,data[self.material2],data['x'],data['y'])
        data=self.material.data_gear_material
        sgl2=self.material.FunCoeff(sgl2,data[self.material2],data['x'],data['y'])
        self.Gear2._CoeffYSIso()
        self.sigma_lim2=(sgl2/(safety_factor*self.Gear2.coeff_ys_iso))*10**7
        
    def _CoeffYFIso(self):
        #facteur de forme pour la contrainte ISO
        angle=30/180*npy.pi
        self.Gear1.GearSectionISO(angle)
        self.Gear2.GearSectionISO(angle)
        self.coeff_yf_iso1=(6*(self.Gear1.h_height_iso/self.Gear1.rack.module)*npy.cos(self.transverse_pressure_angle))/((self.Gear1.s_thickness_iso/self.Gear1.rack.module)**2*npy.cos(self.Gear1.rack.transverse_pressure_angle_T))
        self.coeff_yf_iso2=(6*(self.Gear2.h_height_iso/self.Gear2.rack.module)*npy.cos(self.transverse_pressure_angle))/((self.Gear2.s_thickness_iso/self.Gear2.rack.module)**2*npy.cos(self.Gear2.rack.transverse_pressure_angle_T))
        
    def _CoeffYEIso(self):
        #facteur de conduite pour la contrainte ISO
        self.coeff_ye_iso=1/self.radial_contact_ratio
        
    def _CoeffYBIso(self):
        #facteur de contrefort pour la contrainte ISO
        data=self.material.data_coeff_YB_Iso
        self.coeff_yb_iso=self.material.FunCoeff(self.helix_angle,data['evol_coeff_yb_iso'],data['x'],data['y'])
#        self.coeff_yb_iso=float(self.fun_coeff_yb_iso(self.helix_angle)[1])

    ### Generique
    
    def Mass(self):
        return self.Gear1.Mass()+self.Gear2.Mass()
        
    def Dict(self):
        self.SigmaLewis()

        d={}
        for k,v in self.__dict__.items():
            tv=type(v)
            if tv==npy.int64:
                d[k]=int(v)
            elif tv==npy.float64:
                d[k]=float(v)
            else:
                d[k]=v

        del d['Rack1']
        del d['Rack2']
        
        del d['material']
        
        d['mass']=self.Mass()
        
        d['Gear1']=self.Gear1.Dict()
        d['Gear2']=self.Gear2.Dict()

        return d
    
    def RoundData(self,prec):
        self.Update(self.Z1,self.Z2,
           round(float(self.center_distance),prec),
           round(float(self.transverse_pressure_angle),prec),
           round(float(self.helix_angle),prec),
           round(float(self.coefficient_profile_shift1),prec),
           round(float(self.coefficient_profile_shift2),prec),
           round(float(self.gear_width),prec),
           self.maximum_torque,self.material1,self.material2,self.nb_cycle1,
           self.rim1,self.rim2,self.alpha_rim1,self.alpha_rim2,self.boring_diameter1,
           self.boring_diameter2,self.thickness_rim1,self.thickness_rim2,
           transverse_pressure_angle_rack_T1=None,
           transverse_pressure_angle_rack_T2=None,circular_tooth_thickness_rack1=None,
           circular_tooth_thickness_rack2=None,gear_addendum_rack1=None,gear_addendum_rack2=None,
           gear_dedendum_rack1=None,gear_dedendum_rack2=None,root_radius_T1=None,root_radius_T2=None,
           root_radius_R1=None,root_radius_R2=None,transverse_pressure_angle_rack_R1=None,
           transverse_pressure_angle_rack_R2=None)
        for key,ob in (self.Dict()).items():
            if 'float' in str(ob.__class__):
                data=round(float(ob),prec)
                setattr(self,key,data)
    
    ### Export
    
    def FreeCADExport(self,file_path,position1,position2,python_path,freecad_lib_path,export_types):
        RIM1=self.Gear1.RimContour()
        RIM2=self.Gear2.RimContour()
#        RIM1.MPLPlot()
#        RIM2.MPLPlot()
        gear1=primitives3D.RevolvedProfile(vm.Point3D((position1[0],position1[1],0.5*self.gear_width)),vm.Vector3D((0,0,1)),
                                           vm.Vector3D((0,1,0)),[RIM1],vm.Vector3D((position1[0],position1[1],0)),
                                           vm.Vector3D((0,0,1)),angle=2*math.pi,name='Rim1')
        gear2=primitives3D.RevolvedProfile(vm.Point3D((position2[0],position2[1],0.5*self.gear_width)),vm.Vector3D((0,0,1)),
                                           vm.Vector3D((0,1,0)),[RIM2],vm.Vector3D((position2[0],position2[1],0)),
                                           vm.Vector3D((0,0,1)),angle=2*math.pi,name='Rim2')
        
        # Teeth
        TG1=self.Gear1.GearContours(5)
        TG2=self.Gear2.GearContours(5)
        list_rot=self.InitialPosition()
        L1=self.GearAssemblyTrace([TG1,TG2],[position1,position2],list_rot)
        C1=vm.Contour2D(L1[0])
        C2=vm.Contour2D(L1[1])
#        C1.MPLPlot()
        h1=0.5*self.Gear1.alpha_rim*(self.Gear1.outside_diameter-self.Gear1.root_diameter)
        h2=0.5*self.Gear2.alpha_rim*(self.Gear2.outside_diameter-self.Gear2.root_diameter)
        r1=self.Gear1.root_diameter/2-h1/2
        r2=self.Gear2.root_diameter/2-h2/2
#        print(h1,h2,r1,r2)
        c1=primitives3D.Cylinder((position1[0],position1[1],0.48*self.gear_width),
                                 (0,0,1),r1,1.05*self.gear_width)
        c2=primitives3D.Cylinder((position2[0],position2[1],0.48*self.gear_width),
                                 (0,0,1),r2,1.05*self.gear_width)
        
#        C1int=vm.Contour2D([vm.Circle2D(vm.Point2D(position1),r1)])
#        C2int=vm.Contour2D([vm.Circle2D(vm.Point2D(position2),r2)])
        if self.helix_angle==0.:            
            t1=primitives3D.ExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
                                            vm.Vector3D((0,1,0)),[C1],(0,0,self.gear_width))
            t2=primitives3D.ExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
                                            vm.Vector3D((0,1,0)),[C2],(0,0,self.gear_width))
        else:
            t1=primitives3D.HelicalExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
                                                   vm.Vector3D((0,1,0)),(position1[0],position1[1],0),
                                                   (0,0,self.gear_width),self.DF1*mt.pi/mt.tan(self.helix_angle),
                                                   C1)
            t2=primitives3D.HelicalExtrudedProfile(vm.Point3D((0,0,0)),vm.Vector3D((1,0,0)),
                                                   vm.Vector3D((0,1,0)),(position2[0],position2[1],0),
                                                   (0,0,self.gear_width),-self.DF2*mt.pi/mt.tan(self.helix_angle),
                                                   C2)

        # Creating holes in teeths
        t1=primitives3D.Cut(t1,c1,name='Teeth1')
        t2=primitives3D.Cut(t2,c2,name='Teeth2')

        model=vm.VolumeModel([gear1,t1,gear2,t2])
#        model=vm.VolumeModel([gear1,t1,gear2])
        model.FreeCADExport('python',file_path,'/usr/lib/freecad/lib',export_types)
    
    def CSVExport(self):
        self.SigmaLewis()
        d=self.__dict__.copy()
        (d1,d2)=self.Gear1.CSVExport()
        d1=['gear1_'+i for i in d1]
        (d3,d4)=self.Gear2.CSVExport()
        d3=['gear2_'+i for i in d3]
        return list(d.keys())+d1+d3,list(d.values())+d2+d4
    
        
    def SVGExport(self,name,position1,position2):
        #tuple1 et 2 correspondent a la position des centres
        TG1=self.Gear1.GearContours(5)
        TG2=self.Gear2.GearContours(5)
        list_rot=self.InitialPosition()
        L1=self.GearAssemblyTrace([TG1,TG2],[(0,0),(0,0)],list_rot)
        L2=[]
        L2.append(vm.Circle2D(vm.Point2D(position1),self.DF1/2))
        L2.append(vm.Circle2D(vm.Point2D((self.center_distance,0)),self.DF2/2))
        L2.append(vm.Circle2D(vm.Point2D(position1),self.Gear1.base_diameter/2))
        L2.append(vm.Circle2D(vm.Point2D((self.center_distance,0)),self.Gear2.base_diameter/2))
        L2.append(vm.Circle2D(vm.Point2D(position1),self.pitch_diameter_factory1/2))
        L2.append(vm.Circle2D(vm.Point2D((self.center_distance,0)),self.pitch_diameter_factory2/2))
        L3=[]
        L3.append(vm.Circle2D(vm.Point2D(position1),self.Gear1.root_diameter_active/2))
        L3.append(vm.Circle2D(vm.Point2D((self.center_distance,0)),self.Gear2.root_diameter_active/2))
        #G1=vm.Contour2D(LR)
        #G1.MPLPlot()
        SVG1=LibSvg.SVGTrace(1000)
        SVG1.Convert(L1[0],'gear1','black',1/50,0)
        SVG1.Convert(L1[1],'gear2','red',1/50,0)
        SVG1.Convert(L2,'Construction','blue',1/50,0,'0.1px, 0.1px')
        SVG1.Convert(L3,'Construction','red',1/50,0,'0.1px, 0.1px')
        SVG1.Show(name,{'gear1':{'R':[2*npy.pi/self.Gear1.tooth_number,0,0]},'gear2':{'R':[-2*npy.pi/self.Gear2.tooth_number,self.center_distance,0]}})
        
    def SVGGearSet(self,name,position1,position2,local=0):
        #tuple1 et 2 correspondent a la position des centres
        TG1=self.Gear1.GearContours(5)
        TG2=self.Gear2.GearContours(5)
        
        #Definition de la position angulaire initiale
        list_rot=self.InitialPosition()
        if position2[0]==position1[0]:
            if position2[1]-position1[1]>0:
                angle=npy.pi/2
            else:
                angle=-npy.pi/2
        else:
            angle=-npy.arctan((position2[1]-position1[1])/(position2[0]-position1[0]))
        list_rot[0]=list_rot[0]-angle
        list_rot[1]=list_rot[1]-angle
        
        L1=self.GearAssemblyTrace([TG1,TG2],[(0,0),(0,0)],list_rot)
        L2=[]
        L2.append(vm.Circle2D(vm.Point2D(position1),self.DF1/2))
        L2.append(vm.Circle2D(vm.Point2D(position2),self.DF2/2))
        L2.append(vm.Circle2D(vm.Point2D(position1),self.Gear1.base_diameter/2))
        L2.append(vm.Circle2D(vm.Point2D(position2),self.Gear2.base_diameter/2))
        L2.append(vm.Circle2D(vm.Point2D(position1),self.pitch_diameter_factory1/2))
        L2.append(vm.Circle2D(vm.Point2D(position2),self.pitch_diameter_factory2/2))
        L2.append(vm.Circle2D(vm.Point2D(position1),self.Gear1.root_diameter_active/2))
        L2.append(vm.Circle2D(vm.Point2D(position2),self.Gear2.root_diameter_active/2))
        
        L3=[]
        #quote
        alpha=npy.pi/4
        L3.append(vm.Line2D(vm.Point2D((position2[0]-npy.cos(alpha)*self.DF2/2,position2[1]-npy.sin(alpha)*self.DF2/2)),vm.Point2D((position2[0]+npy.cos(alpha)*self.DF2/2,position2[1]+npy.sin(alpha)*self.DF2/2))))
        alpha=npy.pi/4+0.4
        L3.append(vm.Line2D(vm.Point2D((position2[0]-npy.cos(alpha)*self.Gear2.base_diameter/2,position2[1]-npy.sin(alpha)*self.Gear2.base_diameter/2)),vm.Point2D((position2[0]+npy.cos(alpha)*self.Gear2.base_diameter/2,position2[1]+npy.sin(alpha)*self.Gear2.base_diameter/2))))
        alpha=npy.pi/4+0.8
        L3.append(vm.Line2D(vm.Point2D((position2[0]-npy.cos(alpha)*self.Gear2.root_diameter_active/2,position2[1]-npy.sin(alpha)*self.Gear2.root_diameter_active/2)),vm.Point2D((position2[0]+npy.cos(alpha)*self.Gear2.root_diameter_active/2,position2[1]+npy.sin(alpha)*self.Gear2.root_diameter_active/2))))
        
        alpha=npy.pi/4
        L3.append(vm.Line2D(vm.Point2D((position1[0]-npy.cos(alpha)*self.DF1/2,position1[1]-npy.sin(alpha)*self.DF1/2)),vm.Point2D((position1[0]+npy.cos(alpha)*self.DF1/2,position1[1]+npy.sin(alpha)*self.DF1/2))))
        alpha=npy.pi/4+0.4
        L3.append(vm.Line2D(vm.Point2D((position1[0]-npy.cos(alpha)*self.Gear1.base_diameter/2,position1[1]-npy.sin(alpha)*self.Gear1.base_diameter/2)),vm.Point2D((position1[0]+npy.cos(alpha)*self.Gear1.base_diameter/2,position1[1]+npy.sin(alpha)*self.Gear1.base_diameter/2))))
        alpha=npy.pi/4+0.8
        L3.append(vm.Line2D(vm.Point2D((position1[0]-npy.cos(alpha)*self.Gear1.root_diameter_active/2,position1[1]-npy.sin(alpha)*self.Gear1.root_diameter_active/2)),vm.Point2D((position1[0]+npy.cos(alpha)*self.Gear1.root_diameter_active/2,position1[1]+npy.sin(alpha)*self.Gear1.root_diameter_active/2))))
        
        #G1=vm.Contour2D(LR)
        #G1.MPLPlot()
        
        #definition de la viewbox
#        boxX_min,boxX_max,boxY_min,boxY_max=npy.inf,-npy.inf,npy.inf,-npy.inf
#        for li in Temp:
#            if li['inbox']==0:
#                boxX_min=min(boxX_min,float(min(npy.array(li['data'])[:,0])))
#                boxY_min=min(boxY_min,float(min(npy.array(li['data'])[:,1])))
#                boxX_max=max(boxX_max,float(max(npy.array(li['data'])[:,0])))
#                boxY_max=max(boxY_max,float(max(npy.array(li['data'])[:,1])))
        boxX_min=min(position1[0]-self.Gear1.outside_diameter/2,position2[0]-self.Gear2.outside_diameter/2)*1000
        boxX_max=max(position1[0]+self.Gear1.outside_diameter/2,position2[0]+self.Gear2.outside_diameter/2)*1000
        boxY_min=min(position1[1]-self.Gear1.outside_diameter/2,position2[1]-self.Gear2.outside_diameter/2)*1000
        boxY_max=max(position1[1]+self.Gear1.outside_diameter/2,position2[1]+self.Gear2.outside_diameter/2)*1000
        
        view_x=(boxX_max-boxX_min)
        view_y=(boxY_max-boxY_min)
        width=700
        scale=width/view_x
        height=scale*view_y
        vb1=boxX_min
        vb2=boxY_min
        vb3=boxX_max-boxX_min
        vb4=boxY_max-boxY_min
        
        Temp,ListGear1=self.ConvertBspline(L1[0],0,1,2,1/scale,0)
        Temp2,ListGear2=self.ConvertBspline(L1[1],0,2,2,1/scale,0)
        Temp.extend(Temp2)
        Temp2=self.ConvertGeom(L2,3,2,0.5/scale,0,0)
        Temp.extend(Temp2)
        Temp2=self.ConvertGeom(L3,3,2,0.5/scale,1,0)
        Temp.extend(Temp2)
        
        data=str(Temp)
        data=data.replace(chr(39)+'data'+chr(39),'data')
        data=data.replace(chr(39)+'curve'+chr(39),'curve')
        data=data.replace(chr(39)+'group'+chr(39),'group')
        data=data.replace(chr(39)+'color'+chr(39),'color')
        data=data.replace(chr(39)+'size'+chr(39),'size')
        data=data.replace(chr(39)+'inbox'+chr(39),'inbox')
        data=data.replace(chr(39)+'cx'+chr(39),'cx')
        data=data.replace(chr(39)+'cy'+chr(39),'cy')
        data=data.replace(chr(39)+'r'+chr(39),'r')
        data=data.replace(chr(39)+'dash'+chr(39),'dash')
        
        rot1=2*npy.pi/self.Gear1.tooth_number*180/npy.pi
        pos1_x=position1[0]*1000
        pos1_y=position1[1]*1000
        rot2=-2*npy.pi/self.Gear2.tooth_number*180/npy.pi
        pos2_x=position2[0]*1000
        pos2_y=position2[1]*1000
        
        if local==1:
            with open(name,'w') as file:
                file.write(self.ExportSVGGearSet(str(ListGear1),str(ListGear2),data,width,height,1/scale,0.5/scale,vb1,vb2,vb3,vb4,name,rot1,pos1_x,pos1_y,rot2,pos2_x,pos2_y))
        
#        SVG1=LibSvg.SVGTrace(1000)
#        SVG1.Convert(L1[0],'gear1','black',1/50,0)
#        SVG1.Convert(L1[1],'gear2','red',1/50,0)
#        SVG1.Convert(L2,'Construction','blue',1/50,0,'0.1px, 0.1px')
#        SVG1.Convert(L3,'Construction','red',1/50,0,'0.1px, 0.1px')
#        SVG1.Show(name,{'gear1':{'R':[2*npy.pi/self.Gear1.tooth_number,0,0]},'gear2':{'R':[-2*npy.pi/self.Gear2.tooth_number,self.center_distance,0]}})
        
    def ExportSVGGearSet(self,ListGear1,ListGear2,data,width,height,trait_ep,trait_ep3,vb1,vb2,vb3,vb4,
                         name,rot1,pos1_x,pos1_y,rot2,pos2_x,pos2_y):
        
        env = Environment(loader=PackageLoader('mechanical_components', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))
        
        template = env.get_template('template_animate2.html')
        
        return template.render(ListGear1=ListGear1,ListGear2=ListGear2,list_name=name,data=data,width=width,height=height,vb1=vb1,vb2=vb2,vb3=vb3,vb4=vb4,trait_ep=trait_ep,trait_ep3=trait_ep3,rot1=rot1,pos1_x=pos1_x,pos1_y=pos1_y,rot2=rot2,pos2_x=pos2_x,pos2_y=pos2_y)
        
    def ConvertBspline(self,liste,curve,group,color,size,inbox):
        #Export ensemble de ligne uniquement
        List=[]
        ListTotal=[]
        for i,obj in enumerate(liste):
            dico={}
            dico['data']=[]
            for j,objL in enumerate(obj.primitives):
                dico['data'].append([1000*objL.points[0].vector[0],1000*objL.points[0].vector[1]]) 
            dico['data'].append([1000*objL.points[-1].vector[0],1000*objL.points[-1].vector[1]])
            dico['curve']=curve
            dico['group']=group
            dico['color']=color
            dico['size']=size
            dico['inbox']=inbox
            List.append(dico)
            ListTotal.extend(dico['data'])
#        ListTotal.extend(List[0]['data'])
        return List,ListTotal
    
    def ConvertGeom(self,liste,group,color,size,dash,inbox):
        List=[]
        for i,obj in enumerate(liste):
            dico={}
            if 'Circle2D' in str(obj.__class__):
                dico['cx']=float(obj.center.vector[0])*1000
                dico['cy']=float(obj.center.vector[1])*1000
                dico['r']=float(obj.radius)*1000
                dico['curve']=1
            if 'Line2D' in str(obj.__class__):
                dico['data']=[[float(obj.points[0].vector[0])*1000,float(obj.points[0].vector[1])*1000],[float(obj.points[1].vector[0])*1000,float(obj.points[1].vector[1])*1000]]
                dico['curve']=2
            dico['group']=group
            dico['color']=color
            dico['size']=size
            dico['inbox']=inbox
            dico['dash']=dash
            List.append(dico)
#        ListTotal.extend(List[0]['data'])
        return List
    
    def MeshingSVGExport(self,name,gear):
        if gear=='Z1':
            dent=self.Gear1
            diam_fonct=self.DF1
        elif gear=='Z2':
            dent=self.Gear2
            diam_fonct=self.DF2
        Ldev=[dent._InvoluteTrace(10,0,'T')]
        Ltroc=[dent._TrochoideTrace(20,0,'T')]
        Lpied=[dent._RootCircleTrace(0)]
        Lout=[dent._OutsideTrace(0)]
        Lcomplet=dent.GearContours(3,[-1,0,1])
        sol=dent._Involute(npy.linspace(0,npy.tan(dent.alpha_outside_diameter),20))
        Lconst=[dent._Trace(sol[0],sol[1])]
        Lconst.extend(dent._TrochoideSecondary([0],'T',0.8,40))
        L2=[vm.Circle2D(vm.Point2D((0,0)),dent.base_diameter/2)]
        L2.append(vm.Circle2D(vm.Point2D((0,0)),diam_fonct/2))
        L2.append(vm.Circle2D(vm.Point2D((0,0)),dent.root_diameter_active/2))
        L2.append(vm.Circle2D(vm.Point2D((0,0)),dent.root_diameter/2))
        L2.append(vm.Circle2D(vm.Point2D((0,0)),dent.outside_diameter/2))
        L3=[vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.outside_diameter/2,0)))]
        L4=[vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.root_diameter/2*npy.cos(-dent.root_angle/2),dent.root_diameter/2*npy.sin(-dent.root_angle/2))))]
        L2.append(vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.root_diameter/2*npy.cos(-dent.root_angle/2-dent.phi0),dent.root_diameter/2*npy.sin(-dent.root_angle/2-dent.phi0)))))
        L2.append(vm.Line2D(vm.Point2D((0,0)),vm.Point2D((dent.outside_diameter/2*npy.cos(-dent.root_angle/2-dent.phi_trochoide),dent.outside_diameter/2*npy.sin(-dent.root_angle/2-dent.phi_trochoide)))))
        SVG1=LibSvg.SVGTrace(700)
        SVG1.Convert(Ldev,'Ldev','black',0.005,0)
        SVG1.Convert(Ltroc,'Ltroc','blue',0.005,0)
        SVG1.Convert(Lpied,'Lpied','red',0.005,0)
        SVG1.Convert(Lout,'Lout','green',0.005,0)
        SVG1.Convert(Lcomplet,'Lcomplet','black',0.005,1,'0.01px, 0.08px')
        SVG1.Convert(Lconst,'Gc','blue',0.005,1,'0.01px, 0.08px')
        SVG1.Convert(L2,'L2','black',0.005,1,'0.01px, 0.08px')
        SVG1.Convert(L3,'L3','black',0.005,1)
        SVG1.Convert(L4,'L4','black',0.005,1)
        SVG1.Show(name)



class ContinuousGearAssemblyOptimizer:
    def __init__(self,Z1,Z2,center_distance,transverse_pressure_angle,
           helix_angle,coefficient_profile_shift1,
           coefficient_profile_shift2,gear_width,maximum_torque,
           material1,material2,nb_cycle1,
           rim1,rim2,alpha_rim1,alpha_rim2,boring_diameter1,boring_diameter2,thickness_rim1,thickness_rim2,
           transverse_pressure_angle_rack_T1=None,
           transverse_pressure_angle_rack_T2=None,
           circular_tooth_thickness_rack1=None,circular_tooth_thickness_rack2=None,
           gear_addendum_rack1=None,
           gear_addendum_rack2=None,gear_dedendum_rack1=None,gear_dedendum_rack2=None,
           root_radius_T1=None,root_radius_T2=None,root_radius_R1=None,
           root_radius_R2=None,transverse_pressure_angle_rack_R1=None,transverse_pressure_angle_rack_R2=None):
        
        self.GearAssemblyParam(Z1,Z2,center_distance,transverse_pressure_angle,
           helix_angle,coefficient_profile_shift1,
           coefficient_profile_shift2,gear_width,maximum_torque,
           material1,material2,nb_cycle1,
           rim1,rim2,alpha_rim1,alpha_rim2,boring_diameter1,boring_diameter2,thickness_rim1,thickness_rim2,
           transverse_pressure_angle_rack_T1,
           transverse_pressure_angle_rack_T2,
           circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
           gear_addendum_rack1,
           gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
           root_radius_T1,root_radius_T2,root_radius_R1,
           root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2)

        self.bounds={'center_distance':center_distance,
                 'transverse_pressure_angle':transverse_pressure_angle,
                 'helix_angle':helix_angle,
                 'coefficient_profile_shift1':coefficient_profile_shift1,
                 'coefficient_profile_shift2':coefficient_profile_shift2,
                 'gear_width':gear_width
                 }
        self.Bounds=npy.array([center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,coefficient_profile_shift2,gear_width])
        
        self.xk={'Z1':self.Z1,'Z2':self.Z2,'maximum_torque':self.maximum_torque,
                 'material1':self.material1,'material2':self.material2,'nb_cycle1':self.nb_cycle1,
                 'rim1':rim1,
                 'rim2':rim2,
                 'alpha_rim1':alpha_rim1,
                 'alpha_rim2':alpha_rim2,
                 'boring_diameter1':boring_diameter1,
                 'boring_diameter2':boring_diameter2,
                 'thickness_rim1':thickness_rim1,
                 'thickness_rim2':thickness_rim2,
                 }
        self.xo={'transverse_pressure_angle_rack_T1':transverse_pressure_angle_rack_T1,
                 'transverse_pressure_angle_rack_T2':transverse_pressure_angle_rack_T2,
                 'circular_tooth_thickness_rack1':circular_tooth_thickness_rack1,
                 'circular_tooth_thickness_rack2':circular_tooth_thickness_rack2,
                 'gear_addendum_rack1':gear_addendum_rack1,
                 'gear_addendum_rack2':gear_addendum_rack2,
                 'gear_dedendum_rack1':gear_dedendum_rack1,
                 'gear_dedendum_rack2':gear_dedendum_rack2,
                 'root_radius_T1':root_radius_T1,
                 'root_radius_T2':root_radius_T2,
                 'root_radius_R1':root_radius_R1,
                 'root_radius_R2':root_radius_R2,
                 'transverse_pressure_angle_rack_R1':transverse_pressure_angle_rack_R1,
                 'transverse_pressure_angle_rack_R2':transverse_pressure_angle_rack_R2}
        
        self.InitX0()
        self.save=self.Xu[:]
        self.solutions=[]
        
        self.DefXU()
        self.xt=dict(list(self.xk.items())+list(self.xu.items())+list(self.xo.items()))
        self.GearAssembly=GearAssembly(**self.xt)
        
        self.coeff=6*[1]
        
        self.rim1=rim1
        self.alpha_rim1=alpha_rim1
        self.boring_diameter1=boring_diameter1
        self.thickness_rim1=thickness_rim1
        self.rim2=rim2
        self.alpha_rim2=alpha_rim2
        self.boring_diameter2=boring_diameter2
        self.thickness_rim2=thickness_rim2
        
    def Update(self,Xu):
        
        self.DefXU(Xu)
        self.xt=dict(list(self.xk.items())+list(self.xu.items())+list(self.xo.items()))
        self.DefVar(**self.xt)
        self.GearAssembly.Update(**self.xt)
        self.save=Xu[:]
        
    def GearAssemblyParam(self,Z1,Z2,center_distance,transverse_pressure_angle,
           helix_angle,coefficient_profile_shift1,
           coefficient_profile_shift2,gear_width,maximum_torque,
           material1,material2,nb_cycle1,
           rim1,rim2,alpha_rim1,alpha_rim2,boring_diameter1,boring_diameter2,thickness_rim1,thickness_rim2,
           transverse_pressure_angle_rack_T1,
           transverse_pressure_angle_rack_T2,
           circular_tooth_thickness_rack1,circular_tooth_thickness_rack2,
           gear_addendum_rack1,
           gear_addendum_rack2,gear_dedendum_rack1,gear_dedendum_rack2,
           root_radius_T1,root_radius_T2,root_radius_R1,
           root_radius_R2,transverse_pressure_angle_rack_R1,transverse_pressure_angle_rack_R2):
        
        self.Z1=Z1
        self.Z2=Z2
        self.maximum_torque=maximum_torque
        self.center_distance=[]
        self.transverse_pressure_angle=[]
        self.helix_angle=[]
        self.coefficient_profile_shift1=[]
        self.coefficient_profile_shift2=[]
        self.gear_width=[]
        self.material1=material1
        self.material2=material2
        self.nb_cycle1=nb_cycle1
        
    def DefVar(self,Z1=None,Z2=None,center_distance=None,transverse_pressure_angle=None,
           helix_angle=None,coefficient_profile_shift1=None,
           coefficient_profile_shift2=None,gear_width=None,maximum_torque=None,
           material1=None,material2=None,nb_cycle1=None,
           rim1=None,rim2=None,alpha_rim1=None,alpha_rim2=None,
           boring_diameter1=None,boring_diameter2=None,thickness_rim1=None,thickness_rim2=None,
           transverse_pressure_angle_rack_T1=None,
           transverse_pressure_angle_rack_T2=None,
           circular_tooth_thickness_rack1=None,circular_tooth_thickness_rack2=None,
           gear_addendum_rack1=None,
           gear_addendum_rack2=None,gear_dedendum_rack1=None,gear_dedendum_rack2=None,
           root_radius_T1=None,root_radius_T2=None,root_radius_R1=None,
           root_radius_R2=None,transverse_pressure_angle_rack_R1=None,transverse_pressure_angle_rack_R2=None):
        self.Z1=Z1
        self.Z2=Z2
        self.maximum_torque=maximum_torque
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.helix_angle=helix_angle
        self.coefficient_profile_shift1=coefficient_profile_shift1
        self.coefficient_profile_shift2=coefficient_profile_shift2
        self.gear_width=gear_width
    
    def DefXU(self,Xu=None):
        #Définition des inconnues utilisées dans l'optimiseur
        if type(Xu)==type(npy.array(())):
            self.center_distance=Xu[0,0]
            self.transverse_pressure_angle=Xu[1,0]
            self.helix_angle=Xu[2,0]
            self.coefficient_profile_shift1=Xu[3,0]
            self.coefficient_profile_shift2=Xu[4,0]
            self.gear_width=Xu[5,0]
            
        self.xu={'center_distance':self.center_distance,
                 'transverse_pressure_angle':self.transverse_pressure_angle,
                 'helix_angle':self.helix_angle,
                 'coefficient_profile_shift1':self.coefficient_profile_shift1,
                 'coefficient_profile_shift2':self.coefficient_profile_shift2,
                 'gear_width':self.gear_width}
        self.Xu=npy.array([[(self.center_distance),
                 (self.transverse_pressure_angle),
                 (self.helix_angle),
                 (self.coefficient_profile_shift1),
                 (self.coefficient_profile_shift2),
                 (self.gear_width)]])
        self.Xu=npy.transpose(self.Xu)
        
    def InitX0(self):
        
        self.DefXU()
        for i in self.bounds:
            self.i=((self.bounds[i][1]-self.bounds[i][0])*npy.random.random(1)+self.bounds[i][0])[0]
            self.xu[i]=self.i
        self.xt=dict(list(self.xk.items())+list(self.xu.items())+list(self.xo.items()))
        self.DefVar(**self.xt)
        self.DefXU()
        
    def CriteriaEq(self):
        
        crit=[0]
        return crit
        
    def CriteriaIneq(self):
        crit=[]
        if self.GearAssembly.Gear1.rim in 'rim_gear':
            crit.append(self.GearAssembly.Gear1.rim_diameter_int-self.GearAssembly.Gear1.boring_diameter_out)
        if self.GearAssembly.Gear2.rim in 'rim_gear':
            crit.append(self.GearAssembly.Gear2.rim_diameter_int-self.GearAssembly.Gear2.boring_diameter_out)
        return crit
        
    def Objective(self,x):

        FEQ=self.feq(x)
        FINEQ=self.fineq(x)
        obj=0
#        for k,i in enumerate(FINEQ):
#            obj+=(i/self.coeff[k])**2
        obj+=0.01*((self.GearAssembly.radial_contact_ratio-2)**2)
        obj+=1000*((self.GearAssembly.linear_backlash-0.05*10**(-3))**2)
#        obj+=(FINEQ[-1]-10)**2
#        obj+=(FINEQ[-2]-10)**2
        for i in FINEQ:
            if i < 0:
                obj+=-1000*i
            else:
                obj+=0.0001*i
        
#        for k,i in enumerate(FINEQ):
#            if i<0.1:
#                obj+=(i/self.coeff[k])**2
#            else:
#                obj+=0
        
        return obj
    
    def ObjectiveCMA(self,x):

        for i,dat in enumerate(self.Bounds):
            x[i]=(x[i]+1)/2*(dat[1]-dat[0])+dat[0]
        
        FEQ=self.feq(x)
        FINEQ=self.fineq(x)
        obj=0
#        for k,i in enumerate(FINEQ):
#            obj+=(i/self.coeff[k])**2
        obj+=0.01*((self.GearAssembly.radial_contact_ratio-2)**2)
        obj+=0.01*((self.GearAssembly.linear_backlash-0.07*10**(-3))**2)
#        obj+=(FINEQ[-1]-10)**2
#        obj+=(FINEQ[-2]-10)**2
        for i in FINEQ:
            if i < 0:
                obj+=-1000*i
            else:
                obj+=0.0001*i
        return obj
             
    def fineq(self,x):
        
#        x=npy.transpose([x])
#        self.GearAssembly.Update(x)
        if False in (npy.transpose([x])==self.save):
            x=npy.transpose([x])
            self.Update(x)
        self.save=self.Xu[:]
        
        ineq=[]
        #ineq.extend(self.GearAssembly.Gear1.CriteriaIneq())
        #ineq.extend(self.GearAssembly.Gear2.CriteriaIneq())
#        ineq.extend(self.CriteriaIneq())
        ineq.extend(self.GearAssembly.CriteriaIneq())
        #ineq.extend(self.GearAssembly.Rack1.CriteriaIneq())
        #ineq.extend(self.GearAssembly.Rack2.CriteriaIneq())
        
        return ineq
        
    def feq(self,x):
        
        eq=[0]
#        eq.extend(self.CriteriaEq())
#        eq.extend(self.GearAssembly.CriteriaEq())
        return eq
        
    def Optimize(self):
        
        boucle=20
        i=0
        arret=0
        while i<boucle and arret==0:
            print('Boucle d\'itération locale {}'.format(i))
#            dim=npy.shape(self.save)[0]
#            sol=npy.random.random(dim)
#            x0=(self.GearAssembly.bounds[:,1]-self.GearAssembly.bounds[:,0])*sol+self.GearAssembly.bounds[:,0]
            self.InitX0()
            x0=npy.array(self.Xu[:])
#            self.GearAssembly.Update(x0)
            self.Update(x0)

            cons = ({'type': 'eq','fun' : self.feq},{'type': 'ineq','fun' : self.fineq})
            cons = {'type': 'ineq','fun' : self.fineq}
            opt = {'maxiter':10000,'xtol': 1e-8, 'disp': False}
            cx = minimize(self.Objective, x0, bounds=self.Bounds, options=opt)
            xsol=cx.x
            if min(self.fineq(xsol))>-1e-7:
                self.solutions.append(xsol)
                arret=1

#            x0= 6 * [0]
#            cx = cma.fmin(self.ObjectiveCMA, [0] * 6, 0.5,options={'bounds':[-1,1]})
##                          options={'bounds':[-1,1],'maxiter':5000,'tolfun':1e-10})
#            x=cx[0]
#            for ii,dat in enumerate(self.Bounds):
#                x[ii]=(x[ii]+1)/2*(dat[1]-dat[0])+dat[0]
#            if min(self.fineq(x))>-1e-7:
#                self.solutions.append(x)
#                arret=1
            
                print('Convergence atteinte avec le statut {}, Valeur de la fonctionnelle {}'.format(cx.status,cx.fun))
            i=i+1
            
            
class GearAssemblyOptimizer:
    def __init__(self,Z1,Z2,center_distance,transverse_pressure_angle,helix_angle,coefficient_profile_shift1,coefficient_profile_shift2,gear_width,maximum_torque,material1,material2,nb_cycle1,ratio=None,rim1=None,rim2=None,alpha_rim1=None,alpha_rim2=None,boring_diameter1=None,boring_diameter2=None,thickness_rim1=None,thickness_rim2=None):
        
        self.ratio=ratio
        self.Z1=Z1
        self.Z2=Z2
        self.center_distance=center_distance
        self.transverse_pressure_angle=transverse_pressure_angle
        self.helix_angle=helix_angle
        self.coefficient_profile_shift1=coefficient_profile_shift1
        self.coefficient_profile_shift2=coefficient_profile_shift2
        self.gear_width=gear_width
        self.maximum_torque=maximum_torque
        self.material1=material1
        self.material2=material2
        self.nb_cycle1=nb_cycle1
        
        self.rim1=rim1
        self.rim2=rim2
        self.alpha_rim1=alpha_rim1
        self.alpha_rim2=alpha_rim2
        self.boring_diameter1=boring_diameter1
        self.boring_diameter2=boring_diameter2
        self.thickness_rim1=thickness_rim1
        self.thickness_rim2=thickness_rim2
        
        self.AnalyzeCombination()
        self.solutions=[]
        
    def AnalyzeCombination(self):
        
        mat1=list(npy.arange(self.Z1['min'],self.Z1['max']+1))
        mat2=list(npy.arange(self.Z2['min'],self.Z2['max']+1))
        plex=itertools.product(mat1,mat2)
        self.plex_calcul=[]
        for i in plex:
            drap=0
            ratio=i[0]/i[1]
            if self.ratio is not None:
                if ratio<=self.ratio['max'] and ratio>=self.ratio['min']:
                    drap=1
            else:
                drap=1
            if drap==1:
                Temp=[(i[0],i[0])]
                Temp1={}
                Temp1['Z1']=i[0]
                Temp.append((i[1],i[1]))
                Temp1['Z2']=i[1]
                data=[self.center_distance,self.transverse_pressure_angle,self.helix_angle,self.coefficient_profile_shift1,self.coefficient_profile_shift2,self.gear_width,self.maximum_torque]
                for j in data:
                    Temp.append((j['min'],j['max']))
                Temp1['center_distance']=(self.center_distance['min'],self.center_distance['max'])
                Temp1['transverse_pressure_angle']=(self.transverse_pressure_angle['min'],self.transverse_pressure_angle['max'])
                Temp1['helix_angle']=(self.helix_angle['min'],self.helix_angle['max'])
                Temp1['coefficient_profile_shift1']=(self.coefficient_profile_shift1['min'],self.coefficient_profile_shift1['max'])
                Temp1['coefficient_profile_shift2']=(self.coefficient_profile_shift2['min'],self.coefficient_profile_shift2['max'])
                Temp1['gear_width']=(self.gear_width['min'],self.gear_width['max'])
                Temp1['maximum_torque']=self.maximum_torque['min']
                Temp1['material1']=self.material1['nom']
                Temp1['material2']=self.material2['nom']
                Temp1['nb_cycle1']=self.nb_cycle1['nom']
                Temp1['rim1']=self.rim1['nom']
                Temp1['rim2']=self.rim2['nom']
                Temp1['alpha_rim1']=self.alpha_rim1['nom']
                Temp1['alpha_rim2']=self.alpha_rim2['nom']
                Temp1['boring_diameter1']=self.boring_diameter1['nom']
                Temp1['boring_diameter2']=self.boring_diameter2['nom']
                Temp1['thickness_rim1']=self.thickness_rim1['nom']
                Temp1['thickness_rim2']=self.thickness_rim2['nom']
                self.plex_calcul.append(Temp1)
                
        
    def Optimize(self,callback=lambda x:x):
        lpx=len(self.plex_calcul)
        for ii,i in enumerate(self.plex_calcul):
            print('{}%'.format(ii/lpx*100))
            callback(ii/lpx)
            A1=ContinuousGearAssemblyOptimizer(**i)
            A1.Optimize()
            try:
                xsol=npy.transpose([A1.solutions[-1]])
                A1.DefXU(xsol)
                xt=dict(list(A1.xk.items())+list(A1.xu.items())+list(A1.xo.items()))
                self.solutions.append(GearAssembly(**xt))
            except:
                pass    
            
class GearAssemblyOptimizerWizard:
    
    def __init__(self,data):
        
        dico={}
        for lign in data:
            dico[lign['type']]={}
            for key,val in lign.items():
                if 'type' not in key:
                    dico[lign['type']][key]=val
        
        #Analyse conformite des datas et lancement des calculs
        self.ok=self.AnalyzeDataSet(dico)
        if self.ok==True:
            dico=self.DefaultDataSet(dico)
            self.dico_gear_assembly_optimizer=dico
            
            
    def Optimize(self,callback=lambda x:x):
        
        M1=GearAssemblyOptimizer(**self.dico_gear_assembly_optimizer)
        M1.Optimize(callback)
        self.solutions=M1.solutions
        
            
    def AnalyzeDataSet(self,data):
        
        if 'ratio' in data:
            if 'nom' not in data['ratio']:
                if 'min' not in data['ratio']:
                    if 'Z1' in data:
                        if 'nom' not in data['Z1']:
                            if 'min' not in data['Z1'] or 'max' not in data['Z1']:
                                return False
                        elif 'nom' not in data['Z2']:
                            if 'min' not in data['Z2'] or 'max' not in data['Z2']:
                                return False
                elif 'max' not in data['ratio']:
                    if 'Z1' in data:
                        if 'nom' not in data['Z1']:
                            if 'min' not in data['Z1'] or 'max' not in data['Z1']:
                                return False
                        elif 'nom' not in data['Z2']:
                            if 'min' not in data['Z2'] or 'max' not in data['Z2']:
                                return False
            elif 'nom' not in data['ratio']:
                if 'min' not in data['ratio'] or 'max' not in data['ratio']:
                    return False
                
                    
        dat=['Z1','Z2','center_distance','transverse_pressure_angle','helix_angle','coefficient_profile_shift1','coefficient_profile_shift2','gear_width','maximum_torque']
        for i in dat:
            if 'nom' not in data[i]:
                if 'min' not in data[i] or 'max' not in data[i]:
                    return False
            
        return True
        
    def DefaultDataSet(self,data):
        
        if 'ratio' in data:
            if 'nom' in data['ratio'] and 'err' not in data['ratio']:
                data['ratio']['err']=0.05
            if 'nom' in data['ratio']:
                data['ratio']['min']=data['ratio']['nom']*(1-data['ratio']['err'])
                data['ratio']['max']=data['ratio']['nom']*(1+data['ratio']['err'])
                del(data['ratio']['nom'])
                del(data['ratio']['err'])
        if 'Z1' in data:
            if 'nom' not in data['Z1']:
                if 'min' not in data['Z1']:
                    data['Z1']['min']=10
                elif 'max' not in data['Z1']:
                    data['Z1']['max']=100
            else:
                data['Z1']['min']=data['Z1']['nom']
                data['Z1']['max']=data['Z1']['nom']
                del(data['Z1']['nom'])
        if 'Z2' in data:
            if 'nom' not in data['Z2']:
                if 'min' not in data['Z2']:
                    data['Z2']['min']=10
                elif 'max' not in data['Z2']:
                    data['Z2']['max']=100
            else:
                data['Z2']['min']=data['Z2']['nom']
                data['Z2']['max']=data['Z2']['nom']
                del(data['Z2']['nom'])

        def analyze1(data,key,err,mini,maxi):
            if key not in data:
                data[key]={}
            if 'nom' in data[key] and 'err' not in data[key]:
                data[key]['err']=err
            if 'nom' in data[key]:
                data[key]['min']=data[key]['nom']*(1-data[key]['err'])
                data[key]['max']=data[key]['nom']*(1+data[key]['err'])
                del(data[key]['nom'])
                del(data[key]['err'])
            if 'min' not in data[key]:
                data[key]['min']=mini
                if 'max' not in data[key]:
                    data[key]['max']=maxi
            if 'max' not in data[key]:
                data[key]['max']=maxi
            return data
            
                
        #Data par defaut en mm,deg et Nm
        data=analyze1(data,'helix_angle',2,15,25)
        data=analyze1(data,'transverse_pressure_angle',2,15,25)
        data=analyze1(data,'center_distance',0.05,40,100)
        data=analyze1(data,'coefficient_profile_shift1',0.05,-1,1)
        data=analyze1(data,'coefficient_profile_shift2',0.05,-1,1)
        data=analyze1(data,'gear_width',0.05,15,25)
        data=analyze1(data,'maximum_torque',0,100,150)
        
        liste=['material1','material2','nb_cycle1','rim1','rim2','boring_diameter1','boring_diameter2','thickness_rim1','thickness_rim2','alpha_rim1','alpha_rim2']
        for i in liste:
            if i not in data:
                data[i]={}
                
        def analyze2(data,key,nom):
            if 'nom' not in data[key]:
                data[key]['nom']=nom
            return data
        
        data=analyze2(data,'material1','hardened_alloy_steel')
        data=analyze2(data,'material2','hardened_alloy_steel')
        data=analyze2(data,'nb_cycle1',10000000)
        data=analyze2(data,'rim1','shaft_gear')
        data=analyze2(data,'rim2','shaft_gear')
        data=analyze2(data,'alpha_rim1',1)
        data=analyze2(data,'alpha_rim2',1)
        data=analyze2(data,'boring_diameter1',30*1e-3)
        data=analyze2(data,'boring_diameter2',30*1e-3)
        data=analyze2(data,'thickness_rim1',7*1e-3)
        data=analyze2(data,'thickness_rim2',7*1e-3)
        
        #changement unité
        data['center_distance']['min']= data['center_distance']['min']/1000
        data['center_distance']['max']= data['center_distance']['max']/1000
        data['transverse_pressure_angle']['min']= data['transverse_pressure_angle']['min']/180*npy.pi
        data['transverse_pressure_angle']['max']= data['transverse_pressure_angle']['max']/180*npy.pi
        data['helix_angle']['min']= data['helix_angle']['min']/180*npy.pi
        data['helix_angle']['max']= data['helix_angle']['max']/180*npy.pi
        data['gear_width']['min']= data['gear_width']['min']/1000
        data['gear_width']['max']= data['gear_width']['max']/1000
        data['boring_diameter1']['nom']= data['boring_diameter1']['nom']/1000
        data['boring_diameter2']['nom']= data['boring_diameter2']['nom']/1000
        data['thickness_rim1']['nom']= data['thickness_rim1']['nom']/1000
        data['thickness_rim2']['nom']= data['thickness_rim2']['nom']/1000
        
        return data

class GearAssemblyOptimizationResults(persistent.Persistent):
    
    def __init__(self,gear_assemblies,bounds):
        for obj in gear_assemblies:
            obj.RoundData(5)
        self.solutions=gear_assemblies
        self.input_data=bounds
        self.type='mc_gear_assembly'
#            
#    def Add(self,list_solutions,bounds,family):
#        
#        for i in list_solutions:
#            self.solutions[family]['obj'].append(i)
#            self.solutions[family]['bnds'].append(bounds)
    
    def CSVExport(self,name,opt='w',family='Famille_A'):
        if self.solutions!=[]:
            (temp1,temp2)=self.solutions[0].CSVExport()
            temp=temp1[0]
            for i in temp1[1::]:
                temp+=','+i
            if opt=='a':
                fichier=open(name,'r')
                temp1=fichier.read()
                temp1=temp1.split('\n')[0]
                temp1=temp1.split(',')
                fichier.close()
            fichier=open(name,opt)
            if not opt=='a':
                fichier.write(temp+'\n')
            for GA in self.solutions:
                (temp3,temp4)=GA.CSVExport()
                temp=''
                for i in temp1:
                    add=temp3.index(i)
                    temp+=str(temp4[add])+','
                fichier.write(temp[0:-1]+'\n')
            fichier.close()
            
    def Dict(self):
        d={}
        solutions=[]
        for ga in self.solutions:
            solutions.append(ga.Dict())
        d['solutions']=solutions
        d['input_data']=self.input_data
        return d
    
        
#class GearAssemblyDBClient(ResultsDBClient):
#    def __init__(self,address):
#        ResultsDBClient.__init__(address,'mc_gear_assembly')
