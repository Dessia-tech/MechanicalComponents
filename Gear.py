
import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm,solve,LinAlgError
from scipy.optimize import *

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
        self.tooth_space=CircularPitch-self.ChordalThickness
        
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
        p7=p4.Translation((self.tooth_space,0))
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
        return crit
        
    def Mass(self):
        pass
    
    def VolumeModel(self):
        pass
    
class Gear:
    def __init__(self,tooth_number,rack_data):
        
        self.GearParam(tooth_number,rack_data)
        self.rack=rack_data
        
    def GearParam(self,tooth_number,rack_data):
        
        self.tooth_number=tooth_number
        self.pitch_diameter_factory=rack_data.CircularPitch*tooth_number/npy.pi
        self.external_diameter=self.pitch_diameter_factory+2*rack_data.GearAddendum
        self.internal_diameter=self.pitch_diameter_factory-2*rack_data.GearDedendum
        self.ChordalThickness=rack_data.ChordalThickness
        self.tooth_space=rack_data.tooth_space
        self.PressureAngleT_factory=rack_data.PressureAngleT
        ChordalThickness_angle=self.ChordalThickness/(self.pitch_diameter_factory/2)
        tooth_space_angle=self.tooth_space/(self.pitch_diameter_factory/2)
        self.base_diameter=self.pitch_diameter_factory*npy.cos(self.PressureAngleT_factory)
        temp=self.pitch_diameter_factory-2*rack_data.GearDedendum+2*rack_data.FiletRadiusT-2*rack_data.FiletRadiusT*npy.sin(rack_data.PressureAngleT)
        theta_int=fsolve( (lambda theta:self.base_diameter/2*npy.cos(theta)+self.base_diameter/2*theta*npy.sin(theta)-temp/2) ,0)[0]
        self.internal_diameter_active=(temp/npy.cos(theta_int-npy.arctan(theta_int)))
        alpha_internal_diameter_active=npy.arccos(self.base_diameter/self.internal_diameter_active)
        self.tan_alpha_internal_diameter_active=npy.tan(alpha_internal_diameter_active)
        self.alpha_internal_diameter_active=alpha_internal_diameter_active
        self.alpha_external_diameter=npy.arccos(self.base_diameter/self.external_diameter)
        alpha_pitch_diameter=npy.arccos(self.base_diameter/self.pitch_diameter_factory)
        self.inv_alpha_pitch_diameter=npy.tan(alpha_pitch_diameter)-alpha_pitch_diameter
        self.base_ChordalThickness=self.base_diameter/2*(2*self.ChordalThickness/self.pitch_diameter_factory+2*self.inv_alpha_pitch_diameter)
        self.BasePitch=rack_data.CircularPitch*npy.cos(rack_data.PressureAngleT)
        
    def Update(self,tooth_number,rack_data):
        
        self.GearParam(tooth_number,rack_data)
        
    def CriteriaIneq(self):
        crit=[self.internal_diameter_active-self.base_diameter-0.3*self.BasePitch]
        crit.append(self.external_diameter-self.internal_diameter-0.3*self.BasePitch)
        return crit
        
    def CriteriaEq(self):
        crit=[0]
        return crit
        
    def InvoluteOptimize(self,discret=10):
        #Analytical involute profil
        theta=npy.linspace(self.tan_alpha_internal_diameter_active,npy.tan(self.alpha_external_diameter),discret)
        sol=self.Involute(theta)
        xT=sol[0]
        yT=sol[1]
        sol=self.Involute(-theta)
        xR=sol[0]
        yR=sol[1]
        pT=[vm.Point2D((xT[0],yT[0]))]
        pR=[vm.Point2D((xR[0],yR[0]))]
        for i in range(1,discret):
            pT.append(vm.Point2D((xT[i],yT[i])))
            pR.append(vm.Point2D((xR[i],yR[i])))
        refT=primitives2D.RoundedLines2D(pT,{},False)
        refR=primitives2D.RoundedLines2D(pR,{},False)
        refR=refR.Rotation(vm.Point2D((0,0)),self.base_ChordalThickness*2/self.base_diameter)
        L=[refT,refR]
        for i in range(1,int(self.tooth_number)):
            L.append(refT.Rotation(vm.Point2D((0,0)),i*2*npy.pi/self.tooth_number))
            L.append(refR.Rotation(vm.Point2D((0,0)),i*2*npy.pi/self.tooth_number))
        return L
        
    def Involute(self,alpha):
        
        sol=(self.base_diameter/2*npy.cos(alpha)+self.base_diameter/2*alpha*npy.sin(alpha),
             self.base_diameter/2*npy.sin(alpha)-self.base_diameter/2*alpha*npy.cos(alpha))
        return sol
        
    def UpdateProfilRack(self):
        
        #Initialisation rack position
        repere=self.PosInitRack()
        
        L=self.InvoluteOptimize()
        L.append(self.rack_profil0)
        return L
        
    def IntersecContours(self,liste1,liste2,prec=1e-1):
        jm=liste1[0]
        arret=0
        liste11=[]
        liste12=[]
        liste21=[]
        liste22=[]
        for j in liste1[1::]:
            km=liste2[0]
            incr=0
            for k in liste2[1::]:
                if arret==0:
                    l1=vm.Line2D(jm,j)
                    l2=vm.Line2D(km,k)
                    pt,bl1,bl2=vm.Point2D.LinesIntersection(l1,l2,True)
                if bl1>-prec and bl1<1+prec and bl2>-prec and bl2<1+prec and arret==0:
                    arret=1
                    incr_fin=incr
                if arret==0:
                    incr+=1
                km=k
            if arret!=1 and arret!=2:
                liste11.append(jm)
            if arret==1:
                liste11.append(jm)
                liste11.append(pt)
                liste12=[pt]
                liste12.extend(liste1[npy.size(liste11)::])
                liste21=liste2[0:incr_fin+1]
                liste21.append(pt)
                liste22=[pt]
                liste22.extend(liste2[npy.size(liste21)-1::])
                arret=2
            jm=j
        return liste11,liste12,liste21,liste22
        
    def PosRack(self,alpha):
        
        angle=npy.tan(alpha)-self.rack.PressureAngleT
        pt0=vm.Point2D((self.pitch_diameter_factory/2,0))
        repere=[pt0,pt0.Translation((1,0)),pt0.Translation((0,1))]
        for i in range(3):
            repere[i]=repere[i].Rotation(vm.Point2D((0,0)),angle)
        ptT=vm.Point2D(self.Involute(alpha))
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
        
    def TrochoidalRackEstimation(self):
        
        #Initialisation rack position
        alpha=self.alpha_external_diameter
        repere,rack_profil=self.PosRack(alpha)
        alpha=2*npy.pi/self.tooth_number
        repere,rack_profil=self.EvolRack(alpha,repere[:],rack_profil)
        
        #Primary trochoide
        nb=20
        TrochoideT=[rack_profil.primitives[-2].center]
        TrochoideR=[rack_profil.primitives[-1].center]
        Dist=[]
        for i in range(2*nb-1):
            Alpha=(-4*npy.pi/self.tooth_number)/nb
            repere,rack_profil=self.EvolRack(Alpha,repere[:],rack_profil)
            TrochoideT.append(rack_profil.primitives[-2].center)
            TrochoideR.append(rack_profil.primitives[-1].center)
            Dist.append(norm(TrochoideT[-1].vector-TrochoideR[-1].vector))
            
        return TrochoideT,TrochoideR
        
    def GearContoursRack(self,discret=10):
        
#        TrochoideT,TrochoideR=self.TrochoidalRackEstimation()
#        L=TrochoideT
#        L.extend(TrochoideR)
#        G2=vm.Contour2D(L)
#        G2.MPLPlot()
        
        #Initialisation rack position
        alpha=self.alpha_internal_diameter_active
        repere,self.rack_profilf=self.PosRack(alpha)
        alpha=self.alpha_external_diameter
        repere,self.rack_profil0=self.PosRack(alpha)
        
        alpha=2*npy.pi/self.tooth_number
        repere,rack_profil=self.EvolRack(alpha,repere[:],self.rack_profil0)
        
        #Primary trochoide
        nb=discret
        liste_ligneT=[rack_profil.primitives[2]]
        liste_pointT=[liste_ligneT[0].points[0]]
        liste_ligneR=[rack_profil.primitives[5]]
        liste_pointR=[liste_ligneR[0].points[0]]
        liste_circleT=[rack_profil.primitives[-2]]
        liste_circleR=[rack_profil.primitives[-1]]
        Ltrajet=[rack_profil]
        for i in range(nb-1):
            Alpha=(-4*alpha)/nb
            repere,rack_profil=self.EvolRack(Alpha,repere[:],rack_profil)
            ligne_temp=rack_profil.primitives[2]
            liste_pointT.append(vm.Point2D.LinesIntersection(liste_ligneT[-1],ligne_temp))
            liste_ligneT.append(ligne_temp)
            liste_circleT.append(rack_profil.primitives[-2])
            ligne_temp=rack_profil.primitives[5]
            liste_pointR.append(vm.Point2D.LinesIntersection(liste_ligneR[-1],ligne_temp))
            liste_ligneR.append(ligne_temp)
            liste_circleR.append(rack_profil.primitives[-1])
        primary_trochoide=[]
        for i in liste_circleT:
            primary_trochoide.append(i.center)
        for i in liste_circleR:
            primary_trochoide.append(i.center)
        
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
            
        #Construction arc cercle pied et tete de dent
        angle=self.base_ChordalThickness*2/self.base_diameter-(npy.tan(self.alpha_external_diameter)-self.alpha_external_diameter)
        liste_pointDE=[vm.Point2D((self.external_diameter/2*npy.cos(angle),self.external_diameter/2*npy.sin(angle)))]
        nb=discret
        for i in range(nb):
            liste_pointDE.append(liste_pointDE[-1].Rotation(vm.Point2D((0,0)),-angle/nb))
        angle=npy.pi/self.tooth_number
        #Estimation new internal radius
        self.internal_diameter_factory=self.external_diameter
        for i in liste_ptT:
            if norm(i.vector)<self.internal_diameter_factory/2:
                self.internal_diameter_factory=2*norm(i.vector)
        liste_pointDI=[vm.Point2D((self.internal_diameter_factory/2*npy.cos(angle),self.internal_diameter_factory/2*npy.sin(angle)))]
        delta_angle=angle+(2*npy.pi/self.tooth_number)
        for i in range(nb):
            liste_pointDI.append(liste_pointDI[-1].Rotation(vm.Point2D((0,0)),-delta_angle/nb))
        
        #Reconstruction du profil complet
        prec=10/nb
        listeA,liste12,liste21,liste22=self.IntersecContours(liste_pointDE,liste_pointT[1::],prec)
        liste_int=[]
        liste_int.extend(listeA)
        
        liste31,liste32,liste41,liste42=self.IntersecContours(liste22,liste_ptT,prec)
        liste51,liste52,liste61,liste62=self.IntersecContours(liste22,liste_ptR,prec)
        if npy.size(liste31)>npy.size(liste51):
            listeB=liste51
            liste6=liste62
            liste7=liste_ptT
        else:
            listeB=liste31
            liste6=liste42
            liste7=liste_ptR
        liste_int.extend(listeB)

        listeC,liste72,liste81,liste82=self.IntersecContours(liste6,liste_pointDI,prec)
        liste_int.extend(listeC)
        
        liste91,liste92,listeD,liste102=self.IntersecContours(liste7,liste82,prec)
        liste_int.extend(listeD)
        
        angle=2*npy.pi/self.tooth_number-2*self.base_ChordalThickness/self.base_diameter+(npy.tan(self.alpha_external_diameter)-self.alpha_external_diameter)
        ind1=[vm.Point2D((self.external_diameter/2*npy.cos(-angle),self.external_diameter/2*npy.sin(-angle)))]
        ind1.append(ind1[-1].Rotation(vm.Point2D((0,0)),0.01))
        listeF,liste132,liste141,liste142=self.IntersecContours(liste_pointR[1::],ind1,prec)
        listeF,liste152,listeE,liste162=self.IntersecContours(listeF[-1::-1],liste92,prec)
        listeF=listeF[-1::-1]
        liste_int.extend(listeE)
        liste_int.extend(listeF)
        
        L=[primitives2D.RoundedLines2D(liste_int,{},False)]
        for i in range(int(self.tooth_number)):
            angle=2*npy.pi/self.tooth_number
            L.append(L[-1].Rotation(vm.Point2D((0,0)),angle))
        
        #Definition tools path
        Loutil=[primitives2D.RoundedLines2D(liste_int,{},False)]
        Loutil.extend(liste_ptT)
        Loutil.extend(liste_ptR)
        Loutil.extend(primary_trochoide)
        
        #Definition construction line
        Lconst=[vm.Circle2D(vm.Point2D((0,0)),self.base_diameter/2)]
        Lconst.append(vm.Circle2D(vm.Point2D((0,0)),self.internal_diameter/2))
        Lconst.append(vm.Circle2D(vm.Point2D((0,0)),self.pitch_diameter_factory/2))
        
        return L,Loutil,Lconst,Ltrajet
    
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
              self.Gear1.tooth_number-int(self.Gear1.tooth_number),
              self.Gear2.tooth_number-int(self.Gear2.tooth_number)
              ]
        return crit
        
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
        bounds.append((1,10))
        bounds.append((1,10))
        bounds.append((1,10))
        bounds.append((1,10))
        bounds.append((0.1,5))
        bounds.append((0.1,5))
        bounds.append((0.1,5))
        bounds.append((0.1,5))
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
            
    def MPLPlot(self,x):

        self.AssemblyGear.Update(x)
        #TG1=self.AssemblyGear.AssemblyGear.Gear1.InvoluteOptimize()
        #TG2=self.AssemblyGear.AssemblyGear.Gear2.InvoluteOptimize()
        TG1,Loutil1,Lconst1,Ltrajet1=self.AssemblyGear.AssemblyGear.Gear1.GearContoursRack(100)
        TG2,Loutil2,Lconst2,Ltrajet2=self.AssemblyGear.AssemblyGear.Gear2.GearContoursRack(100)
        TG=[]
        for i in TG2:
            TG.append(i.Translation((self.AssemblyGear.CenterDistance,0)))
        TG.extend(TG1)
        c1=vm.Contour2D(TG)
        c1.MPLPlot()
        
        G2=Loutil1
        G2.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profil0)
        G2.append(self.AssemblyGear.AssemblyGear.Gear1.rack_profilf)
        G2=vm.Contour2D(G2)
        G2.MPLPlot()
        
        G2=Loutil2
        G2.append(self.AssemblyGear.AssemblyGear.Gear2.rack_profil0)
        G2.append(self.AssemblyGear.AssemblyGear.Gear2.rack_profilf)
        G2=vm.Contour2D(G2)
        G2.MPLPlot()
        
        LR=self.AssemblyGear.AssemblyGear.Rack1.RackContours(10)
        R1=vm.Contour2D(LR)
        R1.MPLPlot()
        
        LR=self.AssemblyGear.AssemblyGear.Rack2.RackContours(10)
        R1=vm.Contour2D(LR)
        R1.MPLPlot()
    
        
#Initialisation
A1=AssemblyGearOptimize1()
O1=Optimization(A1)
O1.Optimize()
print('calcul finalisé')
O1.MPLPlot(O1.solution[-1])
    

#CenterDistance=100
#Z1=13
#Z2=46
#CircularPitch=10 #circular pitch on the pitch circle
#PressureAngle=20/180*npy.pi #real presure angle
#
#BasePitch=CircularPitch*npy.cos(PressureAngle)
#DF1=2*CenterDistancea*Z1/Z2/(1+Z1/Z2)
#DF2=2*CenterDistance-DF1
#DB1=DF1*npy.cos(PressureAngle)
#DB2=DF2*npy.cos(PressureAngle)
#
##BaseDiameter=RackPitchDiameter*npy.cos(RackPressureAngle)       
#        
##Rack definition
#PressureAngleRackT1=20/180*npy.pi
#PressureAngleRackT2=19/180*npy.pi
#CircularPitchRack1=BasePitch/npy.cos(PressureAngleRackT1)
#CircularPitchRack2=BasePitch/npy.cos(PressureAngleRackT2)
#Rack1=Rack(PressureAngleRackT1,CircularPitchRack1,CircularPitchRack1*0.55)
#Rack2=Rack(PressureAngleRackT2,CircularPitchRack2)
#LR=Rack1.RackContours(10)
#R1=vm.Contour2D(LR)
#R1.MPLPlot()
#
#Gear1=Gear(Z1,Rack1)
#Gear2=Gear(Z2,Rack2)
#
#TG1=Gear1.InvoluteOptimize()
#TG2=Gear2.InvoluteOptimize()
#TG=[]
#for i in TG2:
#    TG.append(i.Translation((CenterDistance,0)))
#TG.extend(TG1)
#cg=vm.Contour2D(TG)
#cg.MPLPlot()


#LG1,Loutil1,Lconst1,Ltrajet1=Gear1.GearContoursRack(50)
#LG2,Loutil2,Lconst2,Ltrajet2=Gear2.GearContoursRack(50)
#LG=[]
#for i in LG2:
#    LG.append(i.Translation((CenterDistance,0)))


#LG.extend(LG1)
#G1=vm.Contour2D(LG)
#G1.MPLPlot()
#G2=Loutil1
#G2.extend(Ltrajet1)
#G2=vm.Contour2D(G2)
#G2.MPLPlot()
#G2=Loutil2
#G2.extend(Ltrajet2)
#G2=vm.Contour2D(G2)
#G2.MPLPlot()
#
#G=Gear1.UpdateProfilRack()
#G1=vm.Contour2D(G)
#G1.MPLPlot()

#LG1,Loutil1,Lconst1,Ltrajet1=Gear1.GearContoursRack(50)
#G2=Loutil1
#G2.append(Gear1.rack_profil0)
#G2=vm.Contour2D(G2)
#G2.MPLPlot()
