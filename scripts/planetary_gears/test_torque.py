#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 15:55:17 2020

@author: launay
"""



import genmechanics
import genmechanics.linkages as linkages
import genmechanics.loads as loads

import numpy as npy
import scipy.linalg as linalg

Lb=0.05# bearing width
Lgs=0.05# Gearset width

C1=4
C2=-5
w=1
r1=0.048
r2=0.047
y2=0
z2=-0
y3=0
z3=0
Ld=30

Ca=0
Cr=0
Cf=0
Cwb=0# Speed coeff for bearings
Cvgs=0# Speed coeff for gear sets

alpha_gs1=20/360*2*3.1415
beta_gs1=25/360*2*3.1415
alpha_gs2=20/360*2*3.1415
beta_gs2=25/360*2*3.1415

Z1=20
Z3=50
Z2=8
Z4=7
Z6=60
Z5=17
m1=0.0200336
m2=0.0200336
ground=genmechanics.Part('ground')
shaft1=genmechanics.Part('shaft1')
shaft2=genmechanics.Part('shaft2')
shaft3=genmechanics.Part('shaft3')
shaft4=genmechanics.Part('shaft4')
shaft5=genmechanics.Part('shaft5')
shaft6=genmechanics.Part('shaft6')

p1a=npy.array([0.,0.,0.])
p1b=npy.array([Lb+ Lgs,0,0])
p2a=npy.array([0,y2,z2])
p2b=npy.array([Lb+ Lgs,0,0])
p3a=npy.array([0,0,0])
p3b=npy.array([Lb+ Lgs,0,z3])
p6a=npy.array([Ld,y3,z3])
p6b=npy.array([Ld+Lb+Lgs,y3,z3])
v1=npy.array([5,2,0])
v2=npy.array([5,2,0.1])
v3=npy.array([4,7,0])
dir_axis=npy.array([1,0,0])
dgsv=npy.cross(v3-v2,v3-v1)
egs1=genmechanics.geometry.Direction2Euler(dgsv,dir_axis)
print(egs1)
pgs1=npy.array([0,Z1*m1*0.5,0])
pgs2=npy.array([0,Z3*m1*0.5,0])
print(pgs1)
print(pgs2)
pgs3=npy.array([0,Z2*m1*0.5+Z1*m1*0.5,0])
pgs5=npy.array([0,Z2*m1+Z1*m1*0.5,0])
pgs4=npy.array([0,Z2*m1+Z1*m1*0.5+Z4*m1*0.5,0])
pgs6=npy.array([Ld,Z2*m1+Z1*m1*0.5+Z4*m1*0.5,0])
pgs7=npy.array([Ld,Z6*m2*0.5,0])
print(pgs5)
print(pgs7)
dgs1=npy.cross(pgs3-p1a,pgs3-p1b)
egs1=genmechanics.geometry.Direction2Euler([0,0,1])
dgs2=npy.cross(p3a-p2a,p3b-p2b)
print(egs1)
egs2=genmechanics.geometry.Direction2Euler(pgs2,dir_axis)
egs2=[1.57,-1,0]
egs1=genmechanics.geometry.Direction2Euler([0,0,1])
bearing1a=linkages.BallLinkage(ground,shaft1,p1a,[0,0,0],Ca,Cr,Cwb,'bearing1a')
bearing1b=linkages.LinearAnnularLinkage(ground,shaft1,p1b,[0,0,0],Cr,Cwb,'bearing1b')
bearing2a=linkages.BallLinkage(ground,shaft2,p2a,[0,0,0],Ca,Cr,Cwb,'bearing2a')
bearing2b=linkages.LinearAnnularLinkage(ground,shaft2,p2b,[0,0,0],Cr,Cwb,'bearing2b')
bearing3a=linkages.BallLinkage(ground,shaft3,p3a,[0,0,0],Ca,Cr,Cwb,'bearing3a')
bearing3b=linkages.LinearAnnularLinkage(ground,shaft3,p3b,[0,0,0],Cr,Cwb,'bearing3b')
bearing4a=linkages.BallLinkage(ground,shaft6,p6a,[0,0,0],Ca,Cr,Cwb,'bearing4a')
bearing4b=linkages.LinearAnnularLinkage(ground,shaft6,p6b,[0,0,0],Cr,Cwb,'bearing4b')

pivot=linkages.FrictionlessRevoluteLinkage(shaft3,shaft4,pgs3,[0,0,0],'pivot1')
pivot2=linkages.FrictionlessRevoluteLinkage(shaft5,shaft3,pgs4,[0,0,0],'pivot2')
print(egs1)
gearset12=linkages.GearSetLinkage(shaft1,shaft4,pgs1,egs1,alpha_gs1,beta_gs1,Cf,Cvgs,'Gear set 1')
gearset23=linkages.GearSetLinkage(shaft4,shaft5,pgs5,egs1,alpha_gs2,beta_gs2,Cf,Cvgs,'Gear set 2')
gearset45=linkages.GearSetLinkage(shaft5,shaft2,pgs2,egs1,alpha_gs2,beta_gs2,Cf,Cvgs,'Gear set 3')
gearset46=linkages.GearSetLinkage(shaft5,shaft6,pgs7,egs1,alpha_gs2,beta_gs2,Cf,Cvgs,'Gear set 4')

load1=loads.KnownLoad(shaft6,[0,0,0],[0,0,0],[0,0,0],[C1,0,0],'input torque')
load4=loads.KnownLoad(shaft3,[0,0,0],[0,0,0],[0,0,0],[C2,0,0],'input torque2')
load2=loads.SimpleUnknownLoad(shaft1,[0,0,0],[0,0,0],[],[0],'output torque1')
load3=loads.SimpleUnknownLoad(shaft2,[0,0,0],[0,0,0],[],[0],'output torque2')




imposed_speeds=[(bearing1a,0,w),(bearing4a,0,2)]

mech=genmechanics.Mechanism([bearing1a,bearing1b,bearing2a,bearing2b,bearing3a,bearing3b,bearing4a,bearing4b,gearset12,gearset23,gearset45,gearset46,pivot,pivot2],ground,imposed_speeds,[load1,load4],[load2,load3])
list_part=[bearing1a,bearing1b,bearing2a,bearing2b,bearing3a,bearing3b,bearing4a,bearing4b,gearset12,gearset23,gearset45,gearset46,pivot,pivot2]



for l,lv in mech.static_results.items():
    
    for d,v in lv.items():
        print(l.name,d,v)

# print('Cth: ',r1/(1-r1)*r2/(1-r2)*C)# False now

for l,lv in mech.kinematic_results.items():
    for d,v in lv.items():
        print(l.name,d,v)

print('wth: ',(1-r1)/r1*(1-r2)/r2*w) # False now

mech.GlobalSankey()