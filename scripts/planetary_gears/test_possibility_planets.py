#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 11:39:40 2020

@author: launay
"""
import math as m
import scipy.optimize as op
import matplotlib.pyplot as plt

M1=0.12
X2=0
Z=[40,10,10,10,89,130]
Y2=M1*(Z[0]+Z[1])/2
D_mini=90
D_maxi=10000
X2prime=X2*m.cos(380/3)-Y2*m.sin(380/3)
Y2prime=X2*m.sin(380/3)+Y2*m.cos(380/3)

   
def function_2(x_2,Z1,Z2,Z3,Z4,Z5,Z6,X2,Y2,M1):
    Z=[Z1,Z2,Z3,Z4,Z5,Z6]
    f1 = (x_2[0])**2+(x_2[1])**2-(M1*((Z[4]-Z[2])/2))**2
    f2 = (x_2[0]-x_2[2])**2 +(x_2[1]-x_2[3])**2-(M1*(Z[2]+Z[3])/2)**2
    f3 = (x_2[2]-X2)**2 +(x_2[3]-Y2)**2-(M1*(Z[3]+Z[1])/2)**2
    return [f1,f2,f3]
x0=[-0.05,-2.5,0.1,-0.5]

res= op.least_squares(function_2,x0,xtol=0.000000000000001,args=(Z[0],Z[1],Z[2],Z[3],Z[4],Z[5],X2,Y2,0.12))
print(res)


x_2=res.x
f1 = (x_2[0])**2+(x_2[1])**2-(M1*((Z[4]-Z[2])/2))**2
f2 = (x_2[0]-x_2[2])**2 +(x_2[1]-x_2[3])**2-(M1*(Z[2]+Z[3])/2)**2
f3 = (x_2[2]-X2)**2 +(x_2[3]-Y2)**2-(M1*(Z[3]+Z[1])/2)**2
print(f1,f2,f3)

plt.figure()
ax=plt.subplot(aspect='equal')
circle_1=plt.Circle([x_2[0],x_2[1]],M1*Z[2])
ax.add_patch(plt.Circle([x_2[0],x_2[1]],M1*Z[2]/2,color='r',fill=False))
ax.add_patch(plt.Circle([x_2[2],x_2[3]],M1*Z[3]/2,color='b',fill=False))
ax.add_patch(plt.Circle([X2,Y2],M1*Z[1]/2,color='r',fill=False))
ax.add_patch(plt.Circle([0,0],M1*Z[0]/2,color='y',fill=False))
ax.add_patch(plt.Circle([0,0],M1*Z[4]/2,color='y',fill=False))
ax.set_xlim((-10, 10))
ax.set_ylim((-10, 10))
plt.show()