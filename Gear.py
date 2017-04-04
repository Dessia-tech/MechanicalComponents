#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 17:53:41 2017

@author: steven
"""

import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import math
from scipy.linalg import norm,solve,LinAlgError
from scipy.optimize import *

import pyDOE

class Rack:
    def __init__(self,pressure_angle_apparent,pitch_apparent,tooth_thickness=None,gear_addendum=None,gear_slack=None,filet_radius=None):
        self.pressure_angle_apparent=pressure_angle_apparent
        self.pitch_apparent=pitch_apparent
        self.tooth_thickness=tooth_thickness
        self.gear_addendum=gear_addendum
        self.gear_slack=gear_slack
        self.filet_radius=filet_radius
        modulus_apparent=pitch_apparent/npy.pi
        if gear_addendum==None:
            self.gear_addendum=modulus_apparent
        if gear_slack==None:
            self.gear_slack=1.25*modulus_apparent
        if filet_radius==None:
            self.filet_radius=0.38*modulus_apparent
        if tooth_thickness==None:
            self.tooth_thickness=pitch_apparent/2
        self.tooth_space=pitch_apparent-self.tooth_thickness
    
    def RackContours(self,number):
        p1=vm.Point2D((0,0))
        p2=p1.Translation((self.gear_addendum*npy.tan(self.pressure_angle_apparent),self.gear_addendum))
        p4=p1.Translation((self.tooth_thickness,0))
        p3=p4.Translation((-self.gear_addendum*npy.tan(self.pressure_angle_apparent),self.gear_addendum))
        p5=p4.Translation((self.gear_slack*npy.tan(self.pressure_angle_apparent),-self.gear_slack))
        p7=p4.Translation((self.tooth_space,0))
        p6=p7.Translation((-self.gear_slack*npy.tan(self.pressure_angle_apparent),-self.gear_slack))
        L=primitives2D.RoundedLines2D([p1,p2,p3,p4,p5,p6,p7],{4:self.filet_radius,5:self.filet_radius},False)
        Rack_Elem=[L]
        for i in range(number):
            Rack_Elem.append(L.Translation(((i+1)*(p7.vector-p1.vector))))
        return Rack_Elem
    
    def Mass(self):
        pass
    
    def VolumeModel(self):
        pass
    
class Gear:
    def __init__(self,tooth_number,rack_data):
        self.tooth_number=tooth_number
        self.pitch_diameter_factory=rack_data.pitch_apparent*tooth_number/npy.pi
        self.rack=rack_data
        self.external_diameter=self.pitch_diameter_factory+2*rack_data.gear_addendum
        self.internal_diameter=self.pitch_diameter_factory-2*rack_data.gear_slack
        self.internal_diameter_active=self.pitch_diameter_factory-2*rack_data.gear_slack+2*rack_data.filet_radius*npy.tan(npy.pi/4-rack_data.pressure_angle_apparent/2)*npy.sin(npy.pi/2-rack_data.pressure_angle_apparent)
        self.internal_diameter_active=self.pitch_diameter_factory-2*rack_data.gear_slack+2*rack_data.filet_radius-2*rack_data.filet_radius*npy.sin(rack_data.pressure_angle_apparent)
        self.tooth_thickness=rack_data.tooth_thickness
        self.tooth_space=rack_data.tooth_space
        self.pressure_angle_apparent_factory=rack_data.pressure_angle_apparent
        tooth_thickness_angle=self.tooth_thickness/(self.pitch_diameter_factory/2)
        tooth_space_angle=self.tooth_space/(self.pitch_diameter_factory/2)
        self.base_diameter=self.pitch_diameter_factory*npy.cos(self.pressure_angle_apparent_factory)
        alpha_internal_diameter_active=npy.arccos(self.base_diameter/self.internal_diameter_active)
        self.tan_alpha_internal_diameter_active=npy.tan(alpha_internal_diameter_active)
        self.alpha_internal_diameter_active=npy.arccos(self.base_diameter/self.internal_diameter_active)
        self.alpha_external_diameter=npy.arccos(self.base_diameter/self.external_diameter)
        alpha_pitch_diameter=npy.arccos(self.base_diameter/self.pitch_diameter_factory)
        self.inv_alpha_pitch_diameter=npy.tan(alpha_pitch_diameter)-alpha_pitch_diameter
        self.base_tooth_thickness=self.base_diameter/2*(2*self.tooth_thickness/self.pitch_diameter_factory+2*self.inv_alpha_pitch_diameter)

    def InvoluteOptimize(self,discret=10):
        #Analytical involute profil
        theta=npy.linspace(self.tan_alpha_internal_diameter_active,npy.tan(self.alpha_external_diameter),discret)
        xT=self.base_diameter/2*npy.cos(theta)+self.base_diameter/2*theta*npy.sin(theta)
        yT=self.base_diameter/2*npy.sin(theta)-self.base_diameter/2*theta*npy.cos(theta)
        xR=self.base_diameter/2*npy.cos(-theta)+self.base_diameter/2*(-theta)*npy.sin(-theta)
        yR=self.base_diameter/2*npy.sin(-theta)-self.base_diameter/2*(-theta)*npy.cos(-theta)
        pT=[vm.Point2D((xT[0],yT[0]))]
        pR=[vm.Point2D((xR[0],yR[0]))]
        for i in range(1,discret):
            pT.append(vm.Point2D((xT[i],yT[i])))
            pR.append(vm.Point2D((xR[i],yR[i])))
        refT=primitives2D.RoundedLines2D(pT,{},False)
        refR=primitives2D.RoundedLines2D(pR,{},False)
        refR=refR.Rotation(vm.Point2D((0,0)),self.base_tooth_thickness*2/self.base_diameter)
        L=[refT,refR]
        for i in range(1,self.tooth_number):
            L.append(refT.Rotation(vm.Point2D((0,0)),i*2*npy.pi/self.tooth_number))
            L.append(refR.Rotation(vm.Point2D((0,0)),i*2*npy.pi/self.tooth_number))
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
        
    def GearContours(self,discret=10):
        
        #Initial position of the rack
        angle_init=npy.tan(self.alpha_external_diameter)-self.rack.pressure_angle_apparent
        angle_fin=npy.tan(self.alpha_internal_diameter_active)-self.rack.pressure_angle_apparent
        pt0=vm.Point2D((self.pitch_diameter_factory/2,0))
        repere=[pt0,pt0.Translation((1,0)),pt0.Translation((0,1))]
        for i in range(3):
            repere[i]=repere[i].Rotation(vm.Point2D((0,0)),angle_init)
        ptT=vm.Point2D((self.external_diameter/2*npy.cos(npy.tan(self.alpha_external_diameter)-self.alpha_external_diameter),self.external_diameter/2*npy.sin(npy.tan(self.alpha_external_diameter)-self.alpha_external_diameter)))
        DistY=npy.dot(ptT.vector-repere[0].vector,repere[2].vector-repere[0].vector)
        DistX=npy.dot(ptT.vector-repere[0].vector,repere[1].vector-repere[0].vector)
        self.decal=DistY+self.rack.tooth_thickness-npy.tan(self.rack.pressure_angle_apparent)*DistX
        temp=self.rack.RackContours(0)[0]
        temp=temp.Rotation(vm.Point2D((0,0)),-npy.pi/2)
        temp=temp.Translation((self.pitch_diameter_factory/2,self.decal))
        self.rack_profil=temp.Rotation(vm.Point2D((0,0)),angle_init)
        
        #Angular evolution of the rack
        avance=2*npy.pi/self.tooth_number
        angle_init=0
        angle_fin=-3*avance
        delta=avance*self.pitch_diameter_factory/2
        deltaU=-(repere[2].vector-repere[0].vector)*delta
        for j in range(3):
            repere[j]=repere[j].Rotation(vm.Point2D((0,0)),avance)
        temp=self.rack_profil
        temp=temp.Translation((deltaU))
        self.rack_profil=temp.Rotation(vm.Point2D((0,0)),avance)
            
        #Primary trochoide
        nb=discret
        liste=[self.rack_profil]
        liste_ligneT=[self.rack_profil.primitives[2]]
        liste_pointT=[liste_ligneT[0].points[0]]
        liste_ligneR=[self.rack_profil.primitives[5]]
        liste_pointR=[liste_ligneR[0].points[0]]
        liste_circleT=[self.rack_profil.primitives[-2]]
        liste_circleR=[self.rack_profil.primitives[-1]]
        for i in range(nb-1):
            angle=(angle_fin-angle_init-avance)/nb
            delta=angle*self.pitch_diameter_factory/2
            deltaU=-(repere[2].vector-repere[0].vector)*delta
            for j in range(3):
                repere[j]=repere[j].Rotation(vm.Point2D((0,0)),angle)
            export1=[]
            for j in liste:
                export1.append(j.Translation((deltaU)))
            export2=[]
            for j in export1:
                export2.append(j.Rotation(vm.Point2D((0,0)),angle))
            #for j in export2:
                #L.append(vm.Circle2D(j,self.rack.filet_radius))
                #L.append(j)
            liste=export2
            ligne_temp=liste[0].primitives[2]
            liste_pointT.append(vm.Point2D.LinesIntersection(liste_ligneT[-1],ligne_temp))
            liste_ligneT.append(ligne_temp)
            liste_circleT.append(liste[0].primitives[-2])
            ligne_temp=liste[0].primitives[5]
            liste_pointR.append(vm.Point2D.LinesIntersection(liste_ligneR[-1],ligne_temp))
            liste_ligneR.append(ligne_temp)
            liste_circleR.append(liste[0].primitives[-1])
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
        liste_ptT=[vm.Point2D(cm.vector+v*self.rack.filet_radius)]
        liste_ptT.append(vm.Point2D(cp.vector+v*self.rack.filet_radius))
        for i in liste_circleT[::-1][2::]:
            cp=i.center
            u=cp.vector-cm.vector
            v=npy.array([u[1],-u[0]])/norm(u)
            liste_ptT.append(vm.Point2D(cm.vector+v*self.rack.filet_radius))
            liste_ptT.append(vm.Point2D(cp.vector+v*self.rack.filet_radius))
            cm=cp
        
        #Secondary trochoide of the rear coast 
        cm=liste_circleR[-1].center
        cp=liste_circleR[-2].center
        u=cp.vector-cm.vector
        v=npy.array([u[1],-u[0]])/norm(u)
        liste_ptR=[vm.Point2D(cm.vector+v*self.rack.filet_radius)]
        liste_ptR.append(vm.Point2D(cp.vector+v*self.rack.filet_radius))
        for i in liste_circleR[::-1][2::]:
            cp=i.center
            u=cp.vector-cm.vector
            v=npy.array([u[1],-u[0]])/norm(u)
            liste_ptR.append(vm.Point2D(cm.vector+v*self.rack.filet_radius))
            liste_ptR.append(vm.Point2D(cp.vector+v*self.rack.filet_radius))
            cm=cp
            
        #Construction arc cercle pied et tete de dent
        angle=self.base_tooth_thickness*2/self.base_diameter-(npy.tan(self.alpha_external_diameter)-self.alpha_external_diameter)
        liste_pointDE=[vm.Point2D((self.external_diameter/2*npy.cos(angle),self.external_diameter/2*npy.sin(angle)))]
        nb=discret
        for i in range(nb):
            liste_pointDE.append(liste_pointDE[-1].Rotation(vm.Point2D((0,0)),-angle/nb))
        angle=npy.pi/self.tooth_number
        liste_pointDI=[vm.Point2D((self.internal_diameter/2*npy.cos(angle),self.internal_diameter/2*npy.sin(angle)))]
        delta_angle=angle+(2*npy.pi/self.tooth_number)
        for i in range(nb):
            liste_pointDI.append(liste_pointDI[-1].Rotation(vm.Point2D((0,0)),-delta_angle/nb))
        
        #Reconstruction du profil complet
        prec=10/nb
        listeA,liste12,liste21,liste22=self.IntersecContours(liste_pointDE,liste_pointT[1::],prec)
        liste_int=listeA
        
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
        
        angle=2*npy.pi/self.tooth_number-2*self.base_tooth_thickness/self.base_diameter+(npy.tan(self.alpha_external_diameter)-self.alpha_external_diameter)
        ind1=[vm.Point2D((self.external_diameter/2*npy.cos(-angle),self.external_diameter/2*npy.sin(-angle)))]
        ind1.append(ind1[-1].Rotation(vm.Point2D((0,0)),0.01))
        listeF,liste132,liste141,liste142=self.IntersecContours(liste_pointR[1::],ind1,prec)
        listeF,liste152,listeE,liste162=self.IntersecContours(listeF[-1::-1],liste92,prec)
        listeF=listeF[-1::-1]
        liste_int.extend(listeE)
        liste_int.extend(listeF)
        
        L=[primitives2D.RoundedLines2D(liste_int,{},False)]
        for i in range(self.tooth_number):
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
        
        return L,Loutil,Lconst
    
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

#Rack definition
Rack1=Rack(21/180*npy.pi,15,8)
LR=Rack1.RackContours(10)
R1=vm.Contour2D(LR)
R1.MPLPlot()

Gear1=Gear(25,Rack1)

LG,Loutil,Lconst=Gear1.GearContours(50)
LG.extend(Lconst)
G1=vm.Contour2D(LG)
G1.MPLPlot()
G2=vm.Contour2D(Loutil)
G2.MPLPlot()



