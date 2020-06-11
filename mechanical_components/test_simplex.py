#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 18:01:09 2020

@author: launay
"""


from cassowary import SimplexSolver,Variable
import scipy.optimize as op

reason_1=0.2
reason_2=0.4
speed_3_input=[-400,500]
speed_1_input=[-250,350]
speed_2_input=[-100,800]
# x=[speed_diff_1,speed_diff_2,speed_1,speed_2]
c=[-1,-1,0,0]
A=[[reason_1,(1+reason_1),reason_1,-(1+reason_1)],[reason_1,(1+reason_1),reason_1,1+reason_1],[1,0,1,0],[1,0,-1,0],[0,1,0,1],[0,1,0,-1]]
b=[-speed_3_input[0],speed_3_input[1],speed_1_input[1],-speed_1_input[0],speed_2_input[1],-speed_2_input[0]]
speed_diff_1_bound=(0,None)
speed_diff_2_bound=(0,None)
speed_1_bound=(None,None)
speed_2_bound=(None,None)
res=op.linprog(c,A_ub=A,b_ub=b,bounds=[speed_diff_1_bound,speed_diff_2_bound,speed_1_bound,speed_2_bound])
print(res.success)
# solver= SimplexSolver()

# speed_diff_1=Variable('speed_diff_1')
# speed_diff_2=Variable('speed_diff_2')
# speed_1=Variable('speed_1')
# speed_2=Variable('speed_2')

# speed_3_input=[350,360]
# speed_1_input=[260,270]
# speed_2_input=[140,150]
# solver.add_constraint(-reason_1*speed_diff_1-(1+reason_1)*speed_diff_2+(1+reason_1)*speed_2 - reason_1*speed_1 >= speed_3_input[0])
# solver.add_constraint(+reason_1*speed_diff_1+(1+reason_1)*speed_diff_2+(1+reason_1)*speed_2 + reason_1*speed_1 <= speed_3_input[1])
# solver.add_constraint(speed_diff_1+speed_1<speed_1_input[1])
# solver.add_constraint(-1*speed_diff_1+speed_1>=speed_1_input[0])
# solver.add_constraint(speed_diff_2+speed_2<=speed_2_input[1])
# solver.add_constraint(-1*speed_diff_2+speed_2>=speed_2_input[0])
# solver.add_constraint(speed_diff_2>0.0)
# print(speed_diff_2.value)

