
import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm,solve,LinAlgError
from scipy.optimize import *
from scipy.interpolate import splprep, splev

import pyDOE

class Rack:
    def __init__(self,PressureAngleT,CircularPitch,ChordalThickness=None,GearAddendum=None,
                 GearDedendum=None,FiletRadiusT=None,FiletRadiusR=None,PressureAngleR=None):
        
        self.RackParam(PressureAngleT,CircularPitch,ChordalThickness,GearAddendum,GearDedendum,
                       FiletRadiusT,FiletRadiusR,PressureAngleR)
        
    def RackParam(self,PressureAngleT,CircularPitch,ChordalThickness,GearAddendum,GearDedendum,
                  FiletRadiusT,FiletRadiusR,PressureAngleR):
    
        self.PressureAngleT=PressureAngleT
        self.PressureAngleR=PressureAngleR
        self.CircularPitch=CircularPitch
        self.ChordalThickness=ChordalThickness
        self.GearAddendum=GearAddendum
        self.GearDedendum=GearDedendum
        self.FiletRadiusT=FiletRadiusT
        self.FiletRadiusR=FiletRadiusR
        self.PressureAngleR=PressureAngleR
        modulus_apparent=CircularPitch/npy.pi
        if GearAddendum==None:
            self.GearAddendum=modulus_apparent
        if GearDedendum==None:
            self.GearDedendum=1.25*modulus_apparent
        if FiletRadiusT==None:
            self.FiletRadiusT=0.38*modulus_apparent
        if ChordalThickness==None:
            self.ChordalThickness=CircularPitch/2
        if PressureAngleR==None:
            self.PressureAngleR=self.PressureAngleT
        if FiletRadiusR==None:
            self.FiletRadiusR=self.FiletRadiusT
        self.ToothSpace=CircularPitch-self.ChordalThickness
        
    def Update(self,PressureAngleT,CircularPitch,ChordalThickness=None,GearAddendum=None,
               GearDedendum=None,FiletRadiusT=None,FiletRadiusR=None,PressureAngleR=None):
        
        self.RackParam(PressureAngleT,CircularPitch,ChordalThickness,GearAddendum,GearDedendum,
                       FiletRadiusT,FiletRadiusR,PressureAngleR)
        
    def RackContours(self,number):
        p1=vm.Point2D((0,0))
        p2=p1.Translation((self.GearAddendum*npy.tan(self.PressureAngleR),self.GearAddendum))
        p4=p1.Translation((self.ChordalThickness,0))
        p3=p4.Translation((-self.GearAddendum*npy.tan(self.PressureAngleT),self.GearAddendum))
        p5=p4.Translation((self.GearDedendum*npy.tan(self.PressureAngleT),-self.GearDedendum))
        p7=p4.Translation((self.ToothSpace,0))
        p6=p7.Translation((-self.GearDedendum*npy.tan(self.PressureAngleR),-self.GearDedendum))
        L=primitives2D.RoundedLines2D([p1,p2,p3,p4,p5,p6,p7],{4:self.FiletRadiusT,5:self.FiletRadiusR},False)
        Rack_Elem=[L]
        for i in range(number):
            Rack_Elem.append(L.Translation(((i+1)*(p7.vector-p1.vector))))
        return Rack_Elem
        
    def CriteriaEq(self):
        crit=[0]
        return crit
        
    def CriteriaIneq(self):
        crit=[self.CircularPitch-self.ChordalThickness-self.GearDedendum*npy.tan(self.PressureAngleT)
        -self.GearDedendum*npy.tan(self.PressureAngleR)-(self.FiletRadiusT*npy.cos(self.PressureAngleT)
        -npy.tan(self.PressureAngleT)*self.FiletRadiusT*(1-npy.sin(self.PressureAngleT)))
        -(self.FiletRadiusR*npy.cos(self.PressureAngleR)-npy.tan(self.PressureAngleR)
        *self.FiletRadiusR*(1-npy.sin(self.PressureAngleR)))]
        crit.append(self.ChordalThickness-(self.GearAddendum*npy.tan(self.PressureAngleT)+self.GearAddendum*npy.tan(self.PressureAngleR)))
        return crit
        
    def Mass(self):
        pass
    
    def VolumeModel(self):
        pass
    
class Gear:
    def __init__(self,ToothNumber,rack_data):
        
        self.GearParam(ToothNumber,rack_data)
        self.rack=rack_data
        
    def GearParam(self,ToothNumber,rack_data):
        
        self.ToothNumber=ToothNumber
        self.pitch_diameter_factory=rack_data.CircularPitch*ToothNumber/npy.pi
        self.external_diameter=self.pitch_diameter_factory+2*rack_data.GearAddendum
        self.internal_diameter=self.pitch_diameter_factory-2*rack_data.GearDedendum
        self.ChordalThickness=rack_data.ChordalThickness
        self.ToothSpace=rack_data.ToothSpace
        self.PressureAngleT_factory=rack_data.PressureAngleT
        ChordalThickness_angle=self.ChordalThickness/(self.pitch_diameter_factory/2)
        ToothSpace_angle=self.ToothSpace/(self.pitch_diameter_factory/2)
        self.BaseDiameter=self.pitch_diameter_factory*npy.cos(self.PressureAngleT_factory)
        temp=self.pitch_diameter_factory-2*rack_data.GearDedendum+2*rack_data.FiletRadiusT-2*rack_data.FiletRadiusT*npy.sin(rack_data.PressureAngleT)
        theta_int=fsolve( (lambda theta:self.BaseDiameter/2*npy.cos(theta)+self.BaseDiameter/2*theta*npy.sin(theta)-temp/2) ,0)[0]
        theta2=fsolve( (lambda theta:self.BaseDiameter/2*npy.cos(rack_data.PressureAngleT-theta)/npy.cos(theta)-temp/2) ,0)[0]
        #print(11,self.BaseDiameter/2*npy.cos(theta_int)+self.BaseDiameter/2*theta_int*npy.sin(theta_int)-temp/2)
        self.internal_diameter_active=(temp/npy.cos(theta_int-npy.arctan(theta_int)))
        self.internal_diameter_active=self.BaseDiameter/npy.cos(theta2)
        self.alt_internal_contact=self.internal_diameter_active/2*npy.sin(rack_data.PressureAngleT-theta2)
        self.theta0=npy.tan(theta2)-rack_data.PressureAngleT
        self.theta2=theta2
        alpha_internal_diameter_active=npy.arccos(self.BaseDiameter/self.internal_diameter_active)
        self.tan_alpha_internal_diameter_active=npy.tan(alpha_internal_diameter_active)
        self.alpha_internal_diameter_active=alpha_internal_diameter_active
        self.alpha_external_diameter=npy.arccos(self.BaseDiameter/self.external_diameter)
        alpha_pitch_diameter=npy.arccos(self.BaseDiameter/self.pitch_diameter_factory)
        self.inv_alpha_pitch_diameter=npy.tan(alpha_pitch_diameter)-alpha_pitch_diameter
        self.BaseChordalThickness=self.BaseDiameter/2*(2*self.ChordalThickness/self.pitch_diameter_factory+2*self.inv_alpha_pitch_diameter)
        self.BasePitch=rack_data.CircularPitch*npy.cos(rack_data.PressureAngleT)
        
    def Update(self,ToothNumber,rack_data):
        
        self.GearParam(ToothNumber,rack_data)
        
    def CriteriaIneq(self):
        crit=[self.internal_diameter_active-self.BaseDiameter-0.3*self.BasePitch]
        crit.append(self.external_diameter-self.internal_diameter-0.3*self.BasePitch)
        #Distance top of the rack
        #DistRI=self.DistRackInvolute()
        #crit.append(DistRI)
        #Top of the gear must be positive
        crit.append(self.BaseChordalThickness*2/self.BaseDiameter-2*(npy.tan(self.alpha_external_diameter)-self.alpha_external_diameter))
        return crit
        
    def CriteriaEq(self):
        crit=[0]
        return crit
        
    def InvoluteContours(self,discret=10,number=None):
        #Analytical involute profil
        if number==None:
            number=int(self.ToothNumber)
        L=[self.InvoluteTrace(discret,0,'T')]
        L.append(self.InvoluteTrace(discret,0,'R'))
        for i in range(1,number):
            L.append(self.InvoluteTrace(discret,i,'T'))
            L.append(self.InvoluteTrace(discret,i,'R'))
        return L
        
    def InvoluteTrace(self,discret,number,ind='T'):
        
        if ind=='T':
            drap=1
        else:
            drap=-1
        theta=npy.linspace(self.tan_alpha_internal_diameter_active,npy.tan(self.alpha_external_diameter),discret)
        theta=npy.linspace(0,npy.tan(self.alpha_external_diameter),discret)
        sol=self.Involute(drap*theta)
        x=sol[0]
        y=sol[1]
        p=[vm.Point2D((x[0],y[0]))]
        for i in range(1,discret):
            p.append(vm.Point2D((x[i],y[i])))
        ref=primitives2D.RoundedLines2D(p,{},False)
        if ind=='T':
            L=ref.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.ToothNumber)
        else:
            L=ref.Rotation(vm.Point2D((0,0)),self.BaseChordalThickness*2/self.BaseDiameter)
            L=L.Rotation(vm.Point2D((0,0)),-number*2*npy.pi/self.ToothNumber)
        return L
        
    def Involute(self,alpha):
        
        sol=(self.BaseDiameter/2*npy.cos(alpha)+self.BaseDiameter/2*alpha*npy.sin(alpha),
             self.BaseDiameter/2*npy.sin(alpha)-self.BaseDiameter/2*alpha*npy.cos(alpha))
        return sol
        
    def UpdateProfilRack(self):
        
        #Initialisation rack position
        repere=self.PosInitRack()
        
        L=self.InvoluteOptimize()
        L.append(self.rack_profil0)
        return L
        
    def PosRack(self,alpha):
        
        angle=npy.tan(alpha)-self.rack.PressureAngleT
        pt0=vm.Point2D((self.pitch_diameter_factory/2,0))
        repere=[pt0,pt0.Translation((1,0)),pt0.Translation((0,1))]
        for i in range(3):
            repere[i]=repere[i].Rotation(vm.Point2D((0,0)),angle)
        ptT=vm.Point2D(self.Involute(npy.tan(alpha)))
        DistY=npy.dot(ptT.vector-repere[0].vector,repere[2].vector-repere[0].vector)
        DistX=npy.dot(ptT.vector-repere[0].vector,repere[1].vector-repere[0].vector)
        decal=DistY+self.rack.ChordalThickness-npy.tan(self.rack.PressureAngleT)*DistX
        temp=self.rack.RackContours(0)[0]
        temp=temp.Rotation(vm.Point2D((0,0)),-npy.pi/2)
        temp=temp.Translation((self.pitch_diameter_factory/2,decal))
        rack_profil=temp.Rotation(vm.Point2D((0,0)),angle)
        return repere,rack_profil
        
    def EvolRack(self,alpha,repere,rack_profil):
        
        #Angular evolution of the rack
        delta=alpha*self.pitch_diameter_factory/2
        deltaU=-(repere[2].vector-repere[0].vector)*delta
        for j in range(3):
            repere[j]=repere[j].Rotation(vm.Point2D((0,0)),alpha)
        temp=rack_profil
        temp=temp.Translation((deltaU))
        rack_profil=temp.Rotation(vm.Point2D((0,0)),alpha)
        return repere,rack_profil 
    
    def DistRackInvolute(self):
        
        repere,rack_profil=self.PosRack(self.alpha_internal_diameter_active)
        
        l1=rack_profil.primitives[5]
        p1=l1.points[0]
        p2=l1.points[1]
        p3=p2.Translation((p2.vector-p1.vector))
        l1=vm.Line2D(p1,p3)
        
        L=self.InvoluteTrace(10,1,'R')
        PT=[]
        for i in L.primitives:
            PT.append(i.points[0])

        tck1,pt1=self.BSpline(PT)
        p8=vm.Point2D.MiddlePoint(p1,p2)
        p9=vm.Point2D.MiddlePoint(p2,p3)
        tck2,pt2=self.BSpline([p1,p8,p2,p9,p3])        
        fun = (lambda t : norm(npy.array([splev(t[0],tck1)[0],splev(t[0],tck1)[1]])-npy.array([splev(t[1],tck2)[0],splev(t[1],tck2)[1]]))**2)
        bnds = ((0, 1), (0,1))
        sol=minimize(fun,[0,0], bounds=bnds)
        xsol=sol.x
        dist=fun(xsol)
        
#        L=[primitives2D.RoundedLines2D(PT,{},False)]
#        L.append(l1)
#        L.append(vm.Point2D(splev(xsol[0],tck1)))
#        L.append(vm.Point2D(splev(xsol[1],tck2)))
#        G2=vm.Contour2D(L)
#        G2.MPLPlot()
        
        return dist
        
    def TrochoidalRackEstimation(self,discret=10):
        
        #Initialisation rack position
        alpha=self.alpha_internal_diameter_active
        repere,self.rack_profilf=self.PosRack(alpha)
        alpha=self.alpha_external_diameter
        repere,self.rack_profil0=self.PosRack(alpha)
        
        alpha=2.5*npy.pi/self.ToothNumber
        repere,rack_profil=self.EvolRack(alpha,repere[:],self.rack_profil0)
        
        #Primary trochoide
        nb=discret
        liste_ligneT=[rack_profil.primitives[2]]
        liste_pointT=[liste_ligneT[0].points[0]]
        liste_ligneR=[rack_profil.primitives[5]]
        liste_pointR=[liste_ligneR[0].points[0]]
        liste_circleT=[rack_profil.primitives[-2]]
        liste_circleR=[rack_profil.primitives[-1]]
        Ltrajet=[rack_profil,rack_profil.Translation(-(self.ChordalThickness+self.ToothSpace)*(repere[2].vector-repere[0].vector))]
        for i in range(nb-1):
            Alpha=(-5*alpha)/nb
            repere,rack_profil=self.EvolRack(Alpha,repere[:],rack_profil)
            ligne_temp=rack_profil.primitives[2]
            liste_pointT.append(vm.Point2D.LinesIntersection(liste_ligneT[-1],ligne_temp))
            liste_ligneT.append(ligne_temp)
            liste_circleT.append(rack_profil.primitives[-2])
            ligne_temp=rack_profil.primitives[5]
            liste_pointR.append(vm.Point2D.LinesIntersection(liste_ligneR[-1],ligne_temp))
            liste_ligneR.append(ligne_temp)
            liste_circleR.append(rack_profil.primitives[-1])
            Ltrajet.append(rack_profil)
            Ltrajet.append(rack_profil.Translation(-(self.ChordalThickness+self.ToothSpace)*(repere[2].vector-repere[0].vector)))
        TrochoidalPrimaryT=[]
        TrochoidalPrimaryR=[]
        for i in liste_circleT:
            TrochoidalPrimaryT.append(i.center)
        for i in liste_circleR:
            TrochoidalPrimaryR.append(i.center)
        self.TrochoidalPrimaryControlT,pt1=self.BSpline(TrochoidalPrimaryT)
        self.TrochoidalPrimaryControlR,pt1=self.BSpline(TrochoidalPrimaryR)
        self.Ltrajet=Ltrajet
        
        #Secondary trochoide of the drive coast 
        cm=liste_circleT[-1].center
        cp=liste_circleT[-2].center
        u=cp.vector-cm.vector
        v=npy.array([u[1],-u[0]])/norm(u)
        liste_ptT=[vm.Point2D(cm.vector+v*self.rack.FiletRadiusT)]
        liste_ptT.append(vm.Point2D(cp.vector+v*self.rack.FiletRadiusT))
        for i in liste_circleT[::-1][2::]:
            cp=i.center
            u=cp.vector-cm.vector
            v=npy.array([u[1],-u[0]])/norm(u)
            liste_ptT.append(vm.Point2D(cm.vector+v*self.rack.FiletRadiusT))
            liste_ptT.append(vm.Point2D(cp.vector+v*self.rack.FiletRadiusT))
            cm=cp
        self.TrochoidalSecondaryControlT,pt1=self.BSpline(liste_ptT[2::])
        
        #Secondary trochoide of the rear coast 
        cm=liste_circleR[-1].center
        cp=liste_circleR[-2].center
        u=cp.vector-cm.vector
        v=npy.array([u[1],-u[0]])/norm(u)
        liste_ptR=[vm.Point2D(cm.vector+v*self.rack.FiletRadiusT)]
        liste_ptR.append(vm.Point2D(cp.vector+v*self.rack.FiletRadiusT))
        for i in liste_circleR[::-1][2::]:
            cp=i.center
            u=cp.vector-cm.vector
            v=npy.array([u[1],-u[0]])/norm(u)
            liste_ptR.append(vm.Point2D(cm.vector+v*self.rack.FiletRadiusT))
            liste_ptR.append(vm.Point2D(cp.vector+v*self.rack.FiletRadiusT))
            cm=cp
        self.TrochoidalSecondaryControlR,pt1=self.BSpline(liste_ptR[2::])
        
        #Construction arc cercle pied et tete de dent
        angle=self.BaseChordalThickness*2/self.BaseDiameter-(npy.tan(self.alpha_external_diameter)-self.alpha_external_diameter)
        liste_pointDE=[vm.Point2D((self.external_diameter/2*npy.cos(angle),self.external_diameter/2*npy.sin(angle)))]
        nb=discret
        for i in range(nb):
            liste_pointDE.append(liste_pointDE[-1].Rotation(vm.Point2D((0,0)),-angle/nb))
        self.DEControl,pt1=self.BSpline(liste_pointDE)
        
        #Estimation new internal radius
        angle=npy.pi/self.ToothNumber
        Control1=self.TrochoidalSecondaryControlT
        fun = (lambda t : norm(npy.array([splev(t,Control1)[0],splev(t,Control1)[1]])))
        sol=minimize(fun,0)
        self.internal_diameter_factory=2*fun(sol.x)
        liste_pointDI=[vm.Point2D((self.internal_diameter_factory/2*npy.cos(angle),self.internal_diameter_factory/2*npy.sin(angle)))]
        delta_angle=angle+(2*npy.pi/self.ToothNumber)
        for i in range(nb):
            liste_pointDI.append(liste_pointDI[-1].Rotation(vm.Point2D((0,0)),-delta_angle/nb))
        self.DIControl,pt1=self.BSpline(liste_pointDI)
        
        return liste_ptT
        
    def BSpline(self,pt):
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
        
    def GearContoursRack(self,discret=10):
        
        discret=50
        post=self.TrochoidalRackEstimation(discret)
        
        L=self.InvoluteTrace(10,0,'T')
        InvoluteL=L
        PT=[]
        for i in L.primitives:
            PT.append(i.points[0])
        PT.append(i.points[1])
        temp=PT
        self.InvoluteControlT,InvolutePT=self.BSpline(PT)
        L=self.InvoluteTrace(10,1,'R')
        PT=[]
        for i in L.primitives:
            PT.append(i.points[0])
        PT.append(i.points[1])
        self.InvoluteControlR,InvolutePR=self.BSpline(PT)
        
        Control1=self.TrochoidalSecondaryControlT
        Control2=self.DIControl
        fun = (lambda t : norm(npy.array([splev(t[0],Control1)[0],splev(t[0],Control1)[1]])-npy.array([splev(t[1],Control2)[0],splev(t[1],Control2)[1]])))
        bnds = ((0, 1), (0,1))
        sol=minimize(fun,[0,0], bounds=bnds, method='SLSQP', tol=1e-20)
        xsolI1=sol.x
        print(fun(xsolI1),xsolI1)
        Control1=self.TrochoidalSecondaryControlR
        Control2=self.DIControl
        fun = (lambda t : norm(npy.array([splev(t[0],Control1)[0],splev(t[0],Control1)[1]])-npy.array([splev(t[1],Control2)[0],splev(t[1],Control2)[1]])))
        bnds = ((0, 1), (0,1))
        sol=minimize(fun,[0,1], bounds=bnds, method='SLSQP', tol=1e-20)
        xsolI2=sol.x
        print(fun(xsolI2),xsolI2)
        
        Control1=self.TrochoidalSecondaryControlT
        Control2=self.InvoluteControlT
        fun = (lambda t : norm(npy.array([splev(t[0],Control1)[0],splev(t[0],Control1)[1]])-npy.array([splev(t[1],Control2)[0],splev(t[1],Control2)[1]]))**2)
        bnds = ((-0.1,xsolI1[0]), (0,1))
        sol=minimize(fun,[0,0.5],bounds=bnds, method='SLSQP', tol=1e-20)
        xsolT=sol.x
        print(fun(xsolT))
        Control1=self.TrochoidalSecondaryControlR
        Control2=self.InvoluteControlR
        fun = (lambda t : norm(npy.array([splev(t[0],Control1)[0],splev(t[0],Control1)[1]])-npy.array([splev(t[1],Control2)[0],splev(t[1],Control2)[1]]))**2)
        bnds = ((xsolI2[0],1.1), (0,1))
        sol=minimize(fun,[1,0.5], bounds=bnds, method='SLSQP', tol=1e-20)
        xsolR=sol.x
        print(fun(xsolR))
        
        Control1=self.DEControl
        Control2=self.InvoluteControlT
        fun = (lambda t : norm(npy.array([splev(t[0],Control1)[0],splev(t[0],Control1)[1]])-npy.array([splev(t[1],Control2)[0],splev(t[1],Control2)[1]]))**2)
        bnds = ((0,1), (0,1))
        sol=minimize(fun,[0,1],bounds=bnds, method='SLSQP', tol=1e-20)
        xsolE1=sol.x
        
        PT=[]
        dis=20
        pas=(xsolT[0]-xsolI1[0])/dis
        for j in range(dis+1):
            i=xsolI1[0]+j*pas
            PT.append(vm.Point2D((splev(i,self.TrochoidalSecondaryControlT))))
        L=[primitives2D.RoundedLines2D(PT,{},False)]
        PT=[]
        pas=(xsolE1[1]-xsolT[1])/dis
        for j in range(dis+1):
            i=xsolT[1]+j*pas
            PT.append(vm.Point2D((splev(i,self.InvoluteControlT))))
        L.append(primitives2D.RoundedLines2D(PT,{},False))
#        PT=[]
#        for i in npy.arange(0,1+pas,pas):
#            PT.append(vm.Point2D((splev(i,self.TrochoidalPrimaryControlT))))
#        L.append(primitives2D.RoundedLines2D(PT,{},False))
        PT=[]
        pas=(xsolI2[0]-xsolR[0])/dis
        for j in range(dis+1):
            i=xsolR[0]+j*pas
            PT.append(vm.Point2D((splev(i,self.TrochoidalSecondaryControlR))))
        L.append(primitives2D.RoundedLines2D(PT,{},False))
        PT=[]
        pas=(1-xsolR[1])/dis
        for j in range(dis+1):
            i=xsolR[1]+j*pas
            PT.append(vm.Point2D((splev(i,self.InvoluteControlR))))
        L.append(primitives2D.RoundedLines2D(PT,{},False))
#        PT=[]
#        for i in npy.arange(0,1+pas,pas):
#            PT.append(vm.Point2D((splev(i,self.TrochoidalPrimaryControlR))))
#        L.append(primitives2D.RoundedLines2D(PT,{},False))
        PT=[]
        pas=(xsolI2[1]-xsolI1[1])/dis
        for j in range(dis+1):
            i=xsolI1[1]+j*pas
            PT.append(vm.Point2D((splev(i,self.DIControl))))
        L.append(primitives2D.RoundedLines2D(PT,{},False))
        PT=[]
        pas=(xsolE1[0])/dis
        for j in range(dis+1):
            i=j*pas
            PT.append(vm.Point2D((splev(i,self.DEControl))))
        L.append(primitives2D.RoundedLines2D(PT,{},False))
        
        Lintersect=[vm.Point2D(splev(xsolT[0],self.TrochoidalSecondaryControlT))]
        Lintersect.append(vm.Point2D(splev(xsolT[1],self.InvoluteControlT)))
        Lintersect.append(vm.Point2D(splev(xsolR[0],self.TrochoidalSecondaryControlR)))
        Lintersect.append(vm.Point2D(splev(xsolR[1],self.InvoluteControlR)))
        Lintersect.append(vm.Point2D(splev(xsolI1[0],self.TrochoidalSecondaryControlT)))
        Lintersect.append(vm.Point2D(splev(xsolI1[1],self.DIControl)))
        Lintersect.append(vm.Point2D(splev(xsolI2[0],self.TrochoidalSecondaryControlR)))
        Lintersect.append(vm.Point2D(splev(xsolI2[1],self.DIControl)))
        Lintersect.append(vm.Point2D(splev(xsolE1[1],self.InvoluteControlT)))
        Lintersect.append(vm.Point2D(splev(xsolE1[0],self.DEControl)))
        Lintersect.extend(L)
        
        G2=vm.Contour2D(Lintersect)
        G2.MPLPlot()
        
        Lint=[]
        for i in range(int(self.ToothNumber)):
            angle=i*2*npy.pi/self.ToothNumber
            temp=[]
            for j in L:
                Lint.append(j.Rotation(vm.Point2D((0,0)),angle))
        
        #Definition tools path
        #Loutil=[primitives2D.RoundedLines2D(liste_int,{},False)]
#        Loutil=Ltrajet
#        Loutil.extend(liste_ptT)
#        Loutil.extend(liste_ptR)
#        Loutil.extend(primary_trochoide)
#        
#        #Definition construction line
#        Lconst=[vm.Circle2D(vm.Point2D((0,0)),self.BaseDiameter/2)]
#        Lconst.append(vm.Circle2D(vm.Point2D((0,0)),self.internal_diameter/2))
#        Lconst.append(vm.Circle2D(vm.Point2D((0,0)),self.pitch_diameter_factory/2))
#        
#        return Loutil,Lconst,Ltrajet
        return Lint
    
    def WheelContours(self):
        pass
    
    def Mass(self):
        return self.VolumeModel().Mass()
    
    def VolumeModel(self):
        p=vm.Point3D((0,0,0))
        x=vm.Vector3D((1,0,0))
        y=vm.Vector3D((0,1,0))
        z=vm.Vector3D((0,0,self.width))
        wheel3D=vm.primitives3D.ExtrudedProfile(p,x,y,self.WheelContours(),z)
        lever3D=vm.primitives3D.ExtrudedProfile(p,x,y,self.WheelContours(),z)
        return vm.VolumeModel([wheel3D,lever3D])
        
class AssemblyGear:
    def __init__(self,Z1,Z2,CenterDistance,CircularPitch,PressureAngle,PressureAngleRackT1,
                 PressureAngleRackT2,ChordalThicknessRack1,ChordalThicknessRack2,GearAddendumRack1,
                 GearAddendumRack2,GearDedendumRack1,GearDedendumRack2,FiletRadiusT1,FiletRadiusT2,FiletRadiusR1,
                 FiletRadiusR2,PressureAngleRackR1,PressureAngleRackR2):
        
        self.AssemblyGearParam(Z1,Z2,CenterDistance,CircularPitch,PressureAngle,PressureAngleRackT1,
                 PressureAngleRackT2,ChordalThicknessRack1,ChordalThicknessRack2,GearAddendumRack1,
                 GearAddendumRack2,GearDedendumRack1,GearDedendumRack2,FiletRadiusT1,FiletRadiusT2,FiletRadiusR1,
                 FiletRadiusR2,PressureAngleRackR1,PressureAngleRackR2)
        self.Rack1=Rack(self.PressureAngleRackT1,self.CircularPitchRack1,ChordalThicknessRack1,
                        GearAddendumRack1,GearDedendumRack1,FiletRadiusT1,FiletRadiusR1,PressureAngleRackR1)
        self.Rack2=Rack(self.PressureAngleRackT2,self.CircularPitchRack2,ChordalThicknessRack2,
                        GearAddendumRack2,GearDedendumRack2,FiletRadiusT2,FiletRadiusR2,PressureAngleRackR2)
        self.Gear1=Gear(self.Z1,self.Rack1)
        self.Gear2=Gear(self.Z2,self.Rack2)
        
    def AssemblyGearParam(self,Z1,Z2,CenterDistance,CircularPitch,PressureAngle,PressureAngleRackT1,
                 PressureAngleRackT2,ChordalThicknessRack1,ChordalThicknessRack2,GearAddendumRack1,
                 GearAddendumRack2,GearDedendumRack1,GearDedendumRack2,FiletRadiusT1,FiletRadiusT2,FiletRadiusR1,
                 FiletRadiusR2,PressureAngleRackR1,PressureAngleRackR2):
        
        self.Z1=Z1
        self.Z2=Z2
        self.CenterDistance=CenterDistance
        self.CircularPitch=CircularPitch
        self.PressureAngle=PressureAngle
        self.PressureAngleRackT1=PressureAngleRackT1
        self.PressureAngleRackT2=PressureAngleRackT2
        self.BasePitch=self.CircularPitch*npy.cos(self.PressureAngle)
        self.DF1=2*self.CenterDistance*self.Z1/self.Z2/(1+self.Z1/self.Z2)
        self.DF2=2*self.CenterDistance-self.DF1
        self.DB1=self.DF1*npy.cos(self.PressureAngle)
        self.DB2=self.DF2*npy.cos(self.PressureAngle)
        self.CircularPitchRack1=self.BasePitch/npy.cos(self.PressureAngleRackT1)
        self.CircularPitchRack2=self.BasePitch/npy.cos(self.PressureAngleRackT2)
        self.ChordalThickness1=ChordalThicknessRack1*npy.cos(self.PressureAngleRackT1)/npy.cos(self.PressureAngle)
        self.ChordalThickness2=ChordalThicknessRack2*npy.cos(self.PressureAngleRackT2)/npy.cos(self.PressureAngle)
        self.SpaceWidth1=self.CircularPitch-self.ChordalThickness1
        self.SpaceWidth2=self.CircularPitch-self.ChordalThickness2
        
    def Update(self,Z1,Z2,CenterDistance,CircularPitch,PressureAngle,PressureAngleRackT1,
                 PressureAngleRackT2,ChordalThicknessRack1,ChordalThicknessRack2,GearAddendumRack1,
                 GearAddendumRack2,GearDedendumRack1,GearDedendumRack2,FiletRadiusT1,FiletRadiusT2,FiletRadiusR1,
                 FiletRadiusR2,PressureAngleRackR1,PressureAngleRackR2):
        
        self.AssemblyGearParam(Z1,Z2,CenterDistance,CircularPitch,PressureAngle,PressureAngleRackT1,
                 PressureAngleRackT2,ChordalThicknessRack1,ChordalThicknessRack2,GearAddendumRack1,
                 GearAddendumRack2,GearDedendumRack1,GearDedendumRack2,FiletRadiusT1,FiletRadiusT2,FiletRadiusR1,
                 FiletRadiusR2,PressureAngleRackR1,PressureAngleRackR2)
        self.Rack1.Update(self.PressureAngleRackT1,self.CircularPitchRack1,ChordalThicknessRack1,
                        GearAddendumRack1,GearDedendumRack1,FiletRadiusT1,FiletRadiusR1,PressureAngleRackR1)
        self.Rack2.Update(self.PressureAngleRackT2,self.CircularPitchRack2,ChordalThicknessRack2,
                        GearAddendumRack2,GearDedendumRack2,FiletRadiusT2,FiletRadiusR2,PressureAngleRackR2)
        self.Gear1.Update(self.Z1,self.Rack1)
        self.Gear2.Update(self.Z2,self.Rack2)
        
    def CriteriaIneq(self):
        crit=[
              2*self.CenterDistance-self.Gear1.internal_diameter_active-self.Gear2.external_diameter,
              2*self.CenterDistance-self.Gear2.internal_diameter_active-self.Gear1.external_diameter,
              self.SpaceWidth1-self.ChordalThickness2,
              self.SpaceWidth2-self.ChordalThickness1
              ]
        return crit
        
    def CriteriaEq(self):
        
        crit=[
              self.Gear1.ToothNumber-int(self.Gear1.ToothNumber),
              self.Gear2.ToothNumber-int(self.Gear2.ToothNumber)
              ]
        return crit
        
    def InitialPosition(self,):
        
        
        
class AssemblyGearOptimize1:
    def __init__(self,x=None):
        
        if x==None:
            x=self.Bounds()
        self.AssemblyGearParam(x)
        self.save=x[:]
        self.AssemblyGear=AssemblyGear(self.Z1,self.Z2,self.CenterDistance,self.CircularPitch,
                                       self.PressureAngle,self.PressureAngleRackT1,self.PressureAngleRackT2,
                                       self.ChordalThicknessRack1,self.ChordalThicknessRack2,self.GearAddendumRack1,
                                       self.GearAddendumRack2,self.GearDedendumRack1,self.GearDedendumRack2,
                                       self.FiletRadiusT1,self.FiletRadiusT2,self.FiletRadiusR1,
                                       self.FiletRadiusR2,self.PressureAngleRackR1,self.PressureAngleRackR2)
        
        
    def Bounds(self):
        
        bounds=[(15,25)] #bound p1 x
        bounds.append((30,56))
        bounds.append((50,100))
        bounds.append((6,12))
        bounds.append((20/180*npy.pi,20/180*npy.pi))
        bounds.append((15/180*npy.pi,30/180*npy.pi))
        bounds.append((15/180*npy.pi,30/180*npy.pi))
        bounds.append((1,10))
        bounds.append((1,10))
        bounds.append((3,10))
        bounds.append((3,10))
        bounds.append((3,10))
        bounds.append((3,10))
        bounds.append((0.3,5))
        bounds.append((0.3,5))
        bounds.append((0.3,5))
        bounds.append((0.3,5))
        bounds.append((15/180*npy.pi,30/180*npy.pi))
        bounds.append((15/180*npy.pi,30/180*npy.pi))
        self.bounds=npy.array(bounds)
        dim=npy.shape(bounds)[0]
        sol=npy.random.random(dim)
        x0=(self.bounds[:,1]-self.bounds[:,0])*sol+self.bounds[:,0]
        return x0
        
    def AssemblyGearParam(self,x):
        
        self.Z1=x[0]
        self.Z2=x[1]
        self.CenterDistance=x[2]
        self.CircularPitch=x[3]
        self.PressureAngle=x[4]
        self.PressureAngleRackT1=x[5]
        self.PressureAngleRackT2=x[6]
        self.ChordalThicknessRack1=x[7]
        self.ChordalThicknessRack2=x[8]
        self.GearAddendumRack1=x[9]
        self.GearAddendumRack2=x[10]
        self.GearDedendumRack1=x[11]
        self.GearDedendumRack2=x[12]
        self.FiletRadiusT1=x[13]
        self.FiletRadiusT2=x[14]
        self.FiletRadiusR1=x[15]
        self.FiletRadiusR2=x[16]
        self.PressureAngleRackR1=x[17]
        self.PressureAngleRackR2=x[18]
        
    def Update(self,x):
        
        self.AssemblyGearParam(x)
        self.AssemblyGear.Update(self.Z1,self.Z2,self.CenterDistance,self.CircularPitch,
                                       self.PressureAngle,self.PressureAngleRackT1,self.PressureAngleRackT2,
                                       self.ChordalThicknessRack1,self.ChordalThicknessRack2,self.GearAddendumRack1,
                                       self.GearAddendumRack2,self.GearDedendumRack1,self.GearDedendumRack2,
                                       self.FiletRadiusT1,self.FiletRadiusT2,self.FiletRadiusR1,
                                       self.FiletRadiusR2,self.PressureAngleRackR1,self.PressureAngleRackR2)
        self.save=x[:]

    def CriteriaEq(self):
        
        crit=[self.AssemblyGear.Rack1.PressureAngleT-self.AssemblyGear.Rack1.PressureAngleR,
              self.AssemblyGear.Rack1.FiletRadiusT-self.AssemblyGear.Rack1.FiletRadiusR,
              self.AssemblyGear.Rack2.PressureAngleT-self.AssemblyGear.Rack2.PressureAngleR,
              self.AssemblyGear.Rack2.FiletRadiusT-self.AssemblyGear.Rack2.FiletRadiusR]
        return crit
        
    def CriteriaIneq(self):
        
        crit=[0]
        return crit
        
class Optimization:
    def __init__(self,Assembly):
        
        self.AssemblyGear=Assembly
        self.solution=[]
        
    def Objective(self,x):

        self.AssemblyGear.Update(x)
        FEQ=self.feq(x)
        FINEQ=self.fineq(x)
        obj=0
        for i in FINEQ:
            obj=obj+(i**2)

        return obj
             
    def fineq(self,x):
        
        self.AssemblyGear.Update(x)
        ineq=[]
        ineq.extend(self.AssemblyGear.AssemblyGear.Gear1.CriteriaIneq())
        ineq.extend(self.AssemblyGear.AssemblyGear.Gear2.CriteriaIneq())
        ineq.extend(self.AssemblyGear.AssemblyGear.CriteriaIneq())
        ineq.extend(self.AssemblyGear.AssemblyGear.Rack1.CriteriaIneq())
        ineq.extend(self.AssemblyGear.AssemblyGear.Rack2.CriteriaIneq())
        return ineq
    
        
    def feq(self,x):
        
        self.AssemblyGear.Update(x)
        eq=[]
        eq.extend(self.AssemblyGear.AssemblyGear.CriteriaEq())
        eq.extend(self.AssemblyGear.CriteriaEq())
        
        return eq
        
    def Optimize(self):
        
        boucle=50
        i=0
        arret=0
        while i<boucle and arret==0:
            dim=npy.shape(self.AssemblyGear.save)[0]
            sol=npy.random.random(dim)
            x0=(self.AssemblyGear.bounds[:,1]-self.AssemblyGear.bounds[:,0])*sol+self.AssemblyGear.bounds[:,0]
            self.AssemblyGear.Update(x0)
            cons = ({'type': 'eq','fun' : self.feq},{'type': 'ineq','fun' : self.fineq})
            opt = {'maxiter':50}
            try:
                #[x, _, _, imode,_] = fmin_slsqp(objective, x0,bounds=bounds, ieqcons=[const2], full_output=1,args=(p0,dR),iter=1000,iprint=0,acc=1e-1)
                #[x, _, imode,_] = root(const1, x0, method='lm')
                cx = minimize(self.Objective, x0, bounds=self.AssemblyGear.bounds,constraints=cons,options=opt)
                FEQ=self.feq(cx.x)
                FINEQ=self.fineq(cx.x)
                print(FEQ,FINEQ)
                if max(FEQ)<1e-6 and min(FEQ)>-1e-6 and min(FINEQ)>-1e-6:
                    xsol=cx.x
                    self.solution.append(xsol)
                    arret=1
            except:
                print('erreur de nom')
            i=i+1
            
    def MPLPlot(self,x,args):

        if args==0:
            #Trochoidal primary with two teeth
            self.AssemblyGear.Update(x)
            TG1=self.AssemblyGear.AssemblyGear.Gear1.InvoluteContours(10,2)
            TG2=self.AssemblyGear.AssemblyGear.Gear2.InvoluteContours(10,2)
            TrochoideT1,TrochoideR1=self.AssemblyGear.AssemblyGear.Gear1.TrochoidalRackEstimation()
            TrochoideT2,TrochoideR2=self.AssemblyGear.AssemblyGear.Gear2.TrochoidalRackEstimation()
            L=TrochoideT1
            L.extend(TrochoideR1)
            L.extend(TG1)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profil0)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profili)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profilj)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profilf)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.pti)
            L.append(self.AssemblyGear.AssemblyGear.Gear1.ptj)
            G2=vm.Contour2D(L)
            G2.MPLPlot()
            
        if args==1:
        
#            Loutil1,Lconst1,Ltrajet1=self.AssemblyGear.AssemblyGear.Gear1.GearContoursRack(10)
#            Loutil2,Lconst2,Ltrajet2=self.AssemblyGear.AssemblyGear.Gear2.GearContoursRack(10)
            L1=self.AssemblyGear.AssemblyGear.Gear1.GearContoursRack(30)
            L2=self.AssemblyGear.AssemblyGear.Gear2.GearContoursRack(30)
            TG1=self.AssemblyGear.AssemblyGear.Gear1.InvoluteContours(10,2)
            TG2=self.AssemblyGear.AssemblyGear.Gear2.InvoluteContours(10,2)
            TG=[]
            for i in L2:
                TG.append(i.Translation((self.AssemblyGear.CenterDistance,0)))
            TG.extend(L1)
            c1=vm.Contour2D(TG)
            c1.MPLPlot()
            
            G2=[self.AssemblyGear.AssemblyGear.Gear1.rack_profil0]
            G2.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profilf)
            G2.extend(TG1)
            G2.extend(self.AssemblyGear.AssemblyGear.Gear1.Ltrajet)
            G2=vm.Contour2D(G2)
            G2.MPLPlot()
            
            G2=[self.AssemblyGear.AssemblyGear.Gear2.rack_profil0]
            G2.append(self.AssemblyGear.AssemblyGear.Gear2.rack_profilf)
            G2.extend(TG2)
            G2.extend(self.AssemblyGear.AssemblyGear.Gear2.Ltrajet)
            G2=vm.Contour2D(G2)
            G2.MPLPlot()
            
            LR=self.AssemblyGear.AssemblyGear.Rack1.RackContours(10)
            R1=vm.Contour2D(LR)
            R1.MPLPlot()
            
            LR=self.AssemblyGear.AssemblyGear.Rack2.RackContours(10)
            R1=vm.Contour2D(LR)
            R1.MPLPlot()
            
            DistRI=self.AssemblyGear.AssemblyGear.Gear1.DistRackInvolute()
    
        
##Initialisation
#A1=AssemblyGearOptimize1()
#O1=Optimization(A1)
#O1.Optimize()
#print('calcul finalisé')
#O1.MPLPlot(O1.solution[-1],1)
    
sol=npy.array([ 20.        ,  49.        ,  70.55799373,   6.20167685,
         0.34906585,   0.26179939,   0.33565841,   3.17991653,
         2.5373728 ,   3.00000002,   3.        ,   3.00000008,
         3.00000001,   0.30000001,   1.08524723,   0.30000001,
         1.08524723,   0.26179939,   0.33565841])
A1=AssemblyGearOptimize1(sol)
O1=Optimization(A1)
O1.MPLPlot(sol,1)

#CenterDistance=100
#Z1=20
#Z2=46
#CircularPitch=6 #circular pitch on the pitch circle
#PressureAngle=18/180*npy.pi #real presure angle
#
#BasePitch=CircularPitch*npy.cos(PressureAngle)
#DF1=2*CenterDistance*Z1/Z2/(1+Z1/Z2)
#DF2=2*CenterDistance-DF1
#DB1=DF1*npy.cos(PressureAngle)
#DB2=DF2*npy.cos(PressureAngle)
#
##BaseDiameter=RackPitchDiameter*npy.cos(RackPressureAngle)       
#        
##Rack definition
#PressureAngleRackT1=21/180*npy.pi
#PressureAngleRackT2=20/180*npy.pi
#CircularPitchRack1=BasePitch/npy.cos(PressureAngleRackT1)
#CircularPitchRack2=BasePitch/npy.cos(PressureAngleRackT2)
#Rack1=Rack(PressureAngleRackT1,CircularPitchRack1,CircularPitchRack1*0.6,1*CircularPitchRack1/npy.pi,0.5*CircularPitchRack1/npy.pi)
#Rack2=Rack(PressureAngleRackT2,CircularPitchRack2)
#LR=Rack1.RackContours(10)
#R1=vm.Contour2D(LR)
#R1.MPLPlot()
#
#Gear1=Gear(Z1,Rack1)
#Gear2=Gear(Z2,Rack2)
#Gear1.GearContoursRack(40)

