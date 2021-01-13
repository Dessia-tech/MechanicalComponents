#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 13:57:07 2020

@author: launay
"""

import scipy.optimize as op
import math as m
import numpy as np
import copy    
import matplotlib.pyplot as plt

number_planet=3
def function(x,Z1,Z2,Z3,Z4,Z5,Z6):
    M1=x[6]
    X2=0
    Z=[Z1,Z2,Z3,Z4,Z5,Z6]
    Y2=M1*(Z[0]+Z[1])/2
    D_mini=11
    D_maxi=200
    X2prime=X2*m.cos(2*m.pi/number_planet)-Y2*m.sin(2*m.pi/number_planet)
    Y2prime=X2*m.sin(2*m.pi/number_planet)+Y2*m.cos(2*m.pi/number_planet)
    X3=x[0]
    Y3=x[1]

    f1 =(x[0])**2+(x[1])**2-(x[6]*((Z[4]-Z[2])/2))**2
    f2 =(X3-X2)**2 +(Y3-Y2)**2-(M1*(Z[2]+Z[1])/2)**2
    

    X3prime=X3*m.cos(2*m.pi/number_planet)-Y3*m.sin(2*m.pi/number_planet)
    Y3prime=X3*m.sin(2*m.pi/number_planet)+Y3*m.cos(2*m.pi/number_planet)

    
    X4=X2
    Y4=Y2

    M2=Y4*2/(Z[5]-Z[3])
    
    X4prime=X4*m.cos(2*m.pi/number_planet)-Y4*m.sin(2*m.pi/number_planet)
    Y4prime=X4*m.sin(2*m.pi/number_planet)+Y4*m.cos(2*m.pi/number_planet)
    min_new=M1
    max_m1=D_maxi/Z[4]
    min_x=[((Z[2])*min_new)**2,((Z[2]+Z[1])*min_new/2)**2,((Z[1])*min_new)**2,((Z[2]+Z[0])*min_new/2)**2,0,((Z[3])*min_new)**2]
    x0=(((Z[2])*min_new)**2,((Z[2]+Z[1])*min_new/2)**2,((Z[1])*min_new)**2,((Z[2]+Z[0])*min_new/2)**2,0,((Z[3])*min_new)**2)
    max_x=[np.inf,np.inf,np.inf,np.inf,((Z[4]-Z[1])*min_new/2)**2,np.inf]
    
    min_max_x=[(min_x[0],max_x[0]),(min_x[1],max_x[1]),(min_x[2],max_x[2]),(min_x[3],max_x[3]),(min_x[4],max_x[4]),(min_x[5],max_x[5])]
    
    def function_4(x_4,X2,Y2,X2prime,Y2prime,X3,Y3,X3prime,Y3prime, X4,Y4,X4prime,Y4prime):
        f6 = (Y3prime-Y3)**2+(X3prime-X3)**2-x_4[0]
        f7 = (Y3prime-Y2)**2+(X3prime-X2)**2-x_4[1]
        f8= (Y2prime-Y2)**2+(X2prime-X2)**2-x_4[2]
    
        f9=X3**2+Y3**2 -x_4[3]
        
        f10=-(X2**2+Y2**2) +x_4[4]
   
        f11= (Y4prime-Y4)**2+(X4prime-X4)**2-x_4[5]
        return f6**2+f7**2+f8**2+f9**2+f10**2+f11**2
    
    res_4= op.minimize(function_4,x0,bounds=min_max_x,args=(X2,Y2,X2prime,Y2prime,X3,Y3,X3prime,Y3prime, X4,Y4,X4prime,Y4prime))
    print(res_4)
    print(min_max_x)
    # f12=D_maxi-M1*Z[4]
    f13=(D_maxi-Z[3]*M2)**2-4*(X4**2+Y4**2)
    # print(f8)
    # print(x[3])
    #F=f13
    # for f in F:
    #     if f<0:
    #         F=function(x+100,Z1,Z2,Z3,Z4,Z5,Z6)
    #         return F
    
    print(f1)
    print(f2)

    print(f13)
    # print(f12)
    return (f1/1e-3)**2+(f2/1e-3)**2+(f13)**2
D_mini=10
D_maxi=200
Z=[60,10,11,20,89,130]
max_m1=D_maxi/Z[4]
min_m1=D_mini/Z[4]

min_x=(-np.inf,-np.inf,0,-np.inf,-np.inf,D_mini,min_m1,)

max_x=(np.inf,np.inf,np.inf,np.inf,np.inf,D_maxi,max_m1)
print(max_x)
print(min_x)
x0=(1,1,1,1,1,D_maxi,max_m1)
min_max_x=[(min_x[0],max_x[0]),(min_x[1],max_x[1]),(min_x[2],max_x[2]),(min_x[3],max_x[3]),(min_x[4],max_x[4]),(min_x[5],max_x[5]),(min_x[6],max_x[6])]
Z=[60,10,11,20,89,130]
myfactr = 1e-300000
res= op.minimize(function,x0,bounds=min_max_x,args=(Z[0],Z[1],Z[2],Z[3],Z[4],Z[5]),method='L-BFGS-B',options={'ftol' : myfactr*np.finfo(float).eps ,'gtol':1e-300000 })

print(res.x)
print(res)
x=res.x
X2=0


Y2=x[6]*(Z[0]+Z[1])/2
D_mini=10

X2prime=X2*m.cos(2*m.pi/number_planet)-Y2*m.sin(2*m.pi/number_planet)
Y2prime=X2*m.sin(2*m.pi/number_planet)+Y2*m.cos(2*m.pi/number_planet)

   
def function_2(x_2,Z1,Z2,Z3,Z4,Z5,Z6,X2,Y2,M1):
    Z=[Z1,Z2,Z3,Z4,Z5,Z6]
    f1 =(x_2[0])**2+(x_2[1])**2-(M1*((Z[4]-Z[2])/2))**2
    f2 =(x_2[0]-X2)**2 +(x_2[1]-Y2)**2-(M1*(Z[2]+Z[1])/2)**2
    return [f1,f2]
x0=[0.1,0.2]

res= op.least_squares(function_2,x0,args=(Z[0],Z[1],Z[2],Z[3],Z[4],Z[5],X2,Y2,x[6]))
x_2=res.x
M1=x[6]
f1 =(x_2[0])**2+(x_2[1])**2-(M1*((Z[4]-Z[2])/2))**2
f2 =(x_2[0]-X2)**2 +(x_2[1]-Y2)**2-(M1*(Z[2]+Z[1])/2)**2
print(f1,f2)
X3=res.x[0]
Y3=res.x[1]
X3prime=X3*m.cos(2*m.pi/number_planet)-Y3*m.sin(2*m.pi/number_planet)
Y3prime=X3*m.sin(2*m.pi/number_planet)+Y3*m.cos(2*m.pi/number_planet)
def function_3(x_3,Z1,Z2,Z3,Z4,Z5,Z6,X2,Y2):
    Z=[Z1,Z2,Z3,Z4,Z5,Z6]
    f3 = x_3[1]-X2
    f4 =x_3[2]-Y2
    f5 = x_3[2]-x_3[0]*(Z[5]-Z[3])/2
    return [f3,f4,f5]
x0=[0.05,0.1,0.2]

res_2= op.least_squares(function_3,x0,args=(Z[0],Z[1],Z[2],Z[3],Z[4],Z[5],X2,Y2))

x_3=res_2.x
f3 = x_3[1]-X2
f4 =x_3[2]-Y2
f5 = x_3[2]-x_3[0]*(Z[5]-Z[3])/2
print(f3,f4,f5)

M2=res_2.x[0]
X4=res_2.x[1]
Y4=res_2.x[2]
X4prime=X4*m.cos(2*m.pi/number_planet)-Y4*m.sin(2*m.pi/number_planet)
Y4prime=X4*m.sin(2*m.pi/number_planet)+Y4*m.cos(2*m.pi/number_planet)
print(M2)
def function_4(x_4,X2,Y2,X2prime,Y2prime,X3,Y3,X3prime,Y3prime, X4,Y4,X4prime,Y4prime):
   f6 = (Y3prime-Y3)**2+(X3prime-X3)**2-x_4[0]
   f7 = (Y3prime-Y2)**2+(X3prime-X2)**2-x_4[1]
   f8= (Y2prime-Y2)**2+(X2prime-X2)**2-x_4[2]
   
   f9=X3**2+Y3**2 -x_4[3]
   
   f10=-(X2**2+Y2**2) +x_4[4]
  
   f11= (Y4prime-Y4)**2+(X4prime-X4)**2-x_4[5]
   return [f6,f7,f8,f9,f10,f11]
min_new=x[6]

min_x=(((Z[2])*min_new)**2,((Z[2]+Z[1])*min_new/2)**2,((Z[1])*min_new)**2,((Z[2]+Z[0])*min_new/2)**2,0,((Z[3])*min_new)**2)
x0=(((Z[2])*min_new)**2,((Z[2]+Z[1])*min_new/2)**2,((Z[1])*min_new)**2,((Z[2]+Z[0])*min_new/2)**2,0,((Z[3])*min_new)**2)
max_x=(np.inf,np.inf,np.inf,np.inf,((Z[4]-Z[1])*min_new/2)**2,np.inf)
min_max_x=[min_x,max_x]
res_4= op.least_squares(function_4,x0,bounds=min_max_x,max_nfev=10000,args=(X2,Y2,X2prime,Y2prime,X3,Y3,X3prime,Y3prime, X4,Y4,X4prime,Y4prime))
x_4=res_4.x
print(min_x)
print(max_x)
print(x_4)
f6 = (Y3prime-Y3)**2+(X3prime-X3)**2-x_4[0]
f7 = (Y3prime-Y2)**2+(X3prime-X2)**2-x_4[1]
f8= (Y2prime-Y2)**2+(X2prime-X2)**2-x_4[2]

f9=X3**2+Y3**2 -x_4[3]
print(X3**2+Y3**2)
print(x_4[3])
f10=-(X2**2+Y2**2) +x_4[4]
  
f11= (Y4prime-Y4)**2+(X4prime-X4)**2-x_4[5]
print(f6,f7,f8,f9,f10,f11)
f12=D_maxi-x[6]*Z[4]
f13=(D_maxi-Z[3]*M2)**2-4*(X4**2+Y4**2)

F=[f6,f7,f8,f9,f10,f11,f12,f13]
print(F)
plt.figure()
ax=plt.subplot(aspect='equal')

ax.add_patch(plt.Circle([X3,Y3], M1*Z[2]/2,color='r',fill=False))
ax.add_patch(plt.Circle([X3prime,Y3prime],M1*Z[2]/2,color='b',fill=False))
ax.add_patch(plt.Circle([X2,Y2],M1*Z[1]/2,color='r',fill=False))
ax.add_patch(plt.Circle([X2prime,Y2prime],M1*Z[1]/2,color='r',fill=False))
ax.add_patch(plt.Circle([0,0],M1*Z[0]/2,color='y',fill=False))
ax.add_patch(plt.Circle([0,0],M1*Z[4]/2,color='y',fill=False))
ax.set_xlim((-100, 100))
ax.set_ylim((-100, 100))
plt.show()

def function_max(x,Z1,Z2,Z3,Z4,Z5,Z6):
    X2=0
    Z=[Z1,Z2,Z3,Z4,Z5,Z6]
    Y2=x[0]*(Z[0]+Z[1])/2
    D_mini=90
    
    X2prime=X2*m.cos(380/3)-Y2*m.sin(380/3)
    Y2prime=X2*m.sin(380/3)+Y2*m.cos(380/3)
    
   
    def function_2_max(x_2,Z1,Z2,Z3,Z4,Z5,Z6,X2,Y2,M1):
        Z=[Z1,Z2,Z3,Z4,Z5,Z6]
        f1 =(x_2[0])**2+(x_2[1])**2-(M1*((Z[4]-Z[2])/2))**2
        f2 =(x_2[0]-X2)**2 +(x_2[1]-Y2)**2-(M1*(Z[2]+Z[1])/2)**2
        return [f1,f2]
    x0=[0.1,0.2]

    res= op.least_squares(function_2_max,x0,args=(Z1,Z2,Z3,Z4,Z5,Z6,X2,Y2,x[0]))
    X3=res.x[0]
    Y3=res.x[1]
    X3prime=X3*m.cos(380/3)-Y3*m.sin(380/3)
    Y3prime=X3*m.sin(380/3)+Y3*m.cos(380/3)
    def function_3_max(x_3,Z1,Z2,Z3,Z4,Z5,Z6,X2,Y2):
        f3 = x_3[1]-X2
        f4 =x_3[2]-Y2
        f5 = x_3[2]-x_3[0]*(Z[5]-Z[3])/2
        return [f3,f4,f5]
    x0=[0.05,0.1,0.2]

    res_2= op.least_squares(function_3_max,x0,args=(Z1,Z2,Z3,Z4,Z5,Z6,X2,Y2))
    M2=res_2.x[0]
    X4=res_2.x[1]
    Y4=res_2.x[2]
    X4prime=X4*m.cos(380/3)-Y4*m.sin(380/3)
    Y4prime=X4*m.sin(380/3)+Y4*m.cos(380/3)
    
    f6 = (Y3prime-Y3)**2+(X3prime-X3)**2-x[1]
    f7 = (Y3prime-Y2)**2+(X3prime-X2)**2-x[2]
    f8= (Y2prime-Y2)**2+(X2prime-X2)**2-x[3]

    f9=X3**2+Y3**2 -x[4]
    
    f10=-(X2**2+Y2**2) +x[5]
   
    f11= (Y4prime-Y4)**2+(X4prime-X4)**2-x[6]

    
    F=[f6,f7,f8,f9,f10,f11,]
    
    # for f in F:
    #     if f<0:
    #         F=function(x+100,Z1,Z2,Z3,Z4,Z5,Z6)
    #         return F
    return F
D_mini=10
D_maxi=10000
Z=[60,10,19,20,89,130]
max_m1=D_maxi/Z[4]
min_m1=D_mini/Z[4]

min_x=(min_m1,((Z[2])*min_m1)**2,((Z[2]+Z[1])*min_m1/2)**2,((Z[1])*min_m1)**2,((Z[2]+Z[0])*min_m1/2)**2,((Z[4]-Z[1])*min_m1/2)**2,((Z[3])*min_m1)**2)

max_x=(max_m1,((Z[2])*max_m1)**2,((Z[2]+Z[1])*max_m1/2)**2,((Z[1])*max_m1)**2,((Z[2]+Z[0])*max_m1/2)**2,((Z[4]-Z[1])*max_m1/2)**2,((Z[3])*max_m1)**2)
print(max_x)
print(min_x)
x0=min_x
min_max_x=[min_x,max_x]
Z=[60,10,19,20,89,130]
res= op.least_squares(function_max,x0,bounds=min_max_x,xtol=0.00001,args=(Z[0],Z[1],Z[2],Z[3],Z[4],Z[5]))
print(res.x)
print(res)
x=res.x
X2=0


Y2=x[0]*(Z[0]+Z[1])/2
D_mini=90

X2prime=X2*m.cos(380/3)-Y2*m.sin(380/3)
Y2prime=X2*m.sin(380/3)+Y2*m.cos(380/3)

   
def function_2_max(x_2,Z1,Z2,Z3,Z4,Z5,Z6,X2,Y2,M1):
    Z=[Z1,Z2,Z3,Z4,Z5,Z6]
    f1 =(x_2[0])**2+(x_2[1])**2-(M1*((Z[4]-Z[2])/2))**2
    f2 =(x_2[0]-X2)**2 +(x_2[1]-Y2)**2-(M1*(Z[2]+Z[1])/2)**2
    return [f1,f2]
x0=[0.1,0.2]

res= op.least_squares(function_2_max,x0,args=(Z[0],Z[1],Z[2],Z[3],Z[4],Z[5],X2,Y2,x[0]))
x_2=res.x
M1=x[0]
f1 =(x_2[0])**2+(x_2[1])**2-(M1*((Z[4]-Z[2])/2))**2
f2 =(x_2[0]-X2)**2 +(x_2[1]-Y2)**2-(M1*(Z[2]+Z[1])/2)**2
print(f1,f2)
X3=res.x[0]
Y3=res.x[1]
X3prime=X3*m.cos(380/3)-Y3*m.sin(380/3)
Y3prime=X3*m.sin(380/3)+Y3*m.cos(380/3)
def function_3_max(x_3,Z1,Z2,Z3,Z4,Z5,Z6,X2,Y2):
    Z=[Z1,Z2,Z3,Z4,Z5,Z6]
    f3 = x_3[1]-X2
    f4 =x_3[2]-Y2
    f5 = x_3[2]-x_3[0]*(Z[5]-Z[3])/2
    return [f3,f4,f5]
x0=[0.05,0.1,0.2]

res_2= op.least_squares(function_3_max,x0,args=(Z[0],Z[1],Z[2],Z[3],Z[4],Z[5],X2,Y2))

x_3=res_2.x
f3 = x_3[1]-X2
f4 =x_3[2]-Y2
f5 = x_3[2]-x_3[0]*(Z[5]-Z[3])/2
print(f3,f4,f5)

M2=res_2.x[0]
X4=res_2.x[1]
Y4=res_2.x[2]
X4prime=X4*m.cos(380/3)-Y4*m.sin(380/3)
Y4prime=X4*m.sin(380/3)+Y4*m.cos(380/3)
print(M2)
f6 = (Y3prime-Y3)**2+(X3prime-X3)**2-x[1]
f7 = (Y3prime-Y2)**2+(X3prime-X2)**2-x[2]
f8= (Y2prime-Y2)**2+(X2prime-X2)**2-x[3]
print(Y3prime)
print(Y3)
f9=X3**2+Y3**2 -x[4]

f10=-(X2**2+Y2**2) +x[5]
   
f11= (Y4prime-Y4)**2+(X4prime-X4)**2-x[6]


F=[f6,f7,f8,f9,f10,f11]
print(F)



