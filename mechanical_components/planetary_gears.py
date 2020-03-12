#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:25:05 2020

@author: launay
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as npy
from scipy.linalg import solve


class Relation:
    
    def __init__ (self,name):
        self.name=name
        
    def _get_system_matrix(self):
        if not self._utd_system_equations:
            self._system_matrix, self._system_rhs = self.SystemEquations()
            self._utd_system_equations = True
        return self._system_matrix    
    system_matrix = property(_get_system_matrix)
    
    def _get_system_rhs(self):
        if not self._utd_system_equations:
            self._system_matrix, self._system_rhs = self.SystemEquations()
            self._utd_system_equations = True
        return self._system_rhs
    
    system_rhs = property(_get_system_rhs)
    
    def _get_n_equations(self):
        n_equations = self.system_matrix.shape[0]
        return n_equations
    
    n_equations = property(_get_n_equations)


class Planetary():
    
    def __init__ (self,name,Z,planetary_type):
        self.name=name
        self.Z=Z
        self.planetary_type=planetary_type
        self.p=0
        self.speed=0
        if planetary_type=='sun':
            self.p=1
        else:
            self.p=-1
        
      
        
class Planet():
    
    def __init__ (self,name,planet_type,Z):
        self.name=name
        self.planet_type=planet_type
        self.Z=Z
        self.speed=0
class PlanetCarrier():
    
    def __init__ (self,name):
        self.name=name
        self.speed=0   

class Gearing(Relation):
    
    def __init__(self,name,nodes):
        
        Relation.__init__(self,name)
        self.flag_planet=0
        self.Z_planets=[]
        for node in nodes:
            if type(node)== Planet:
                self.flag_planet+=1
                self.Z_planets.append(node.Z)
            else:
                self.Z_planetary= node.p*node.Z
                
    def SystemEquations(self):
        
        if self.flag_planet == 2:
            matrix=npy.array([ self.Z_planets[0], self.Z_planets[1]])
            rhs= npy.array([0])
            
        else:
            matrix=npy.array([self.Z_planetary ,self.Z_planets[0],-self.Z_planetary])
            rhs= npy.array([0])
            
        return matrix,rhs
    
class Pivot(Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        
    def SystemEquations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs
        
class Fixed (Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        
    def SystemEquations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0, 0])
        return matrix,rhs
        
class Double (Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        
    def SystemEquations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs
        
class Clutch (Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        
    def SystemEquations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs
        
class PlanetaryGears():
    
    def __init__(self,name,planetary_1,planetary_2,planets,planet_carrier):
        
        #Initialize Planetarygears elements
        self.name=name
        self.planetary_1=planetary_1
        self.planetary_2=planetary_2
        self.planets=planets
        self.planet_carrier=planet_carrier
        self.planetaries= [planetary_1,planetary_2]
        self.list_elements=[self.planetary_1]+self.planets+[self.planetary_2,self.planet_carrier]
        
        #Initialize relations
        self.clutch_1=Clutch('Cl1',self.planetary_1)
        self.clutch_2=Clutch('Cl2',self.planetary_2)
        self.clutch_3=Clutch('Cl3',self.planet_carrier)
        self.gearings=[]
        self.pivots=[]
        
        #Initialize gearing , pivots and Double Planets
        flag_double=0
        flag_gearing=0
        
        self.list_relations=[]
        self.list_nodes=[self.clutch_1,self.planetary_1] #list node is use for graph
        for i,planet in enumerate(self.planets):           
            
            if planet.planet_type == 'Double':
                flag_double+=1
                
            if flag_double == 2:
                self.double=Double('Double'+str(i+1),[self.planets[i-1],self.planets[i]])
                self.list_nodes.append(self.double)
                self.list_relations.append(self.double)
                flag_double=0
                
            else:    
                self.gearings.append(Gearing('Ge_'+str(i+1),[self.list_elements[i],self.list_elements[i+1]]))
                self.list_nodes.append(self.gearings[flag_gearing])
                self.list_relations.append(self.gearings[flag_gearing])
                flag_gearing+=1
                
            self.pivots.append(Pivot('Pv'+str(i+1),[planet,self.planet_carrier]))    
            self.list_nodes.append(planet)
            
        self.gearings.append(Gearing('Ge_'+str(i+2),[self.list_elements[i+1], self.list_elements[i+2]]))
        self.list_nodes.extend([self.gearings[flag_gearing],self.planetary_2,self.clutch_2])
        self.list_relations.append(self.gearings[flag_gearing])
        
    def plot(self):
        
        graph_planetary_gear= nx.Graph()
        
        for k in range((len(self.list_nodes)-1)): 
            graph_planetary_gear.add_edge(self.list_nodes[k].name,self.list_nodes[k+1].name)
            
        graph_planetary_gear.add_edge(self.clutch_3.name,self.planet_carrier.name)
        for k in range(len(self.planets)):
            graph_planetary_gear.add_edge(self.planet_carrier.name,self.pivots[k].name)
            graph_planetary_gear.add_edge(self.pivots[k].name,self.planets[k].name)
            
        nx.draw_kamada_kawai(graph_planetary_gear, with_labels=True)
        
        return graph_planetary_gear
    
    def SystemEquations(self):
        #initialize system matrix
        self.n_equations = len(self.list_relations)
        self.n_variables= len(self.list_elements)
        system_matrix = npy.zeros((self.n_equations+2, self.n_variables))
        rhs = npy.zeros(self.n_equations+2)
        num_satelite=0
        num_planetary=0
        for i,relation in enumerate(self.list_relations):
            matrix_relation, rhs_relation = relation.SystemEquations()
            if type(relation)== Gearing :
                if relation.flag_planet==2:
                    system_matrix[i][num_satelite+2]=matrix_relation[0] 
                    system_matrix[i][num_satelite+3]=   matrix_relation[1]
                    
                    rhs[i]=rhs_relation[0]
                    
                    num_satelite+=1
                else:
                    system_matrix[i][num_planetary]=   matrix_relation[0]         
                    system_matrix[i][num_satelite+2]=   matrix_relation[1]
                    system_matrix[i][self.n_variables-1]=   matrix_relation[2]
                    
                    rhs[i]=rhs_relation[0]
                    num_planetary+=1
                
                    
            else:
                system_matrix[i][num_satelite+2]=matrix_relation[0] 
                system_matrix[i][num_satelite+3]=   matrix_relation[1]
                
                rhs[i]=rhs_relation[0]
                
                num_satelite+=1
                
       
        return system_matrix,rhs
    
    def Solve(self,input_speed,input_composant,fixed_composant):
        system_matrix, vector_b = self.SystemEquations()
        list_elements_speed=[self.planetary_1,self.planetary_2]+self.planets+[self.planet_carrier]
        for i,element in enumerate(list_elements_speed):
            if element==input_composant:
                system_matrix[self.n_equations][i]=1
                vector_b[self.n_equations]=input_speed
                
            if element==fixed_composant:
                system_matrix[self.n_equations+1][i]=1
                vector_b[self.n_equations+1]=0
        
        print(system_matrix)          
        solution = solve(system_matrix, vector_b)
        
        print(solution)
        

       

        
    
class AssemblyPlanetaryGears():
    
    def __init__(self,name,planetary_gears,list_link):
        self.planetary_gears=planetary_gears
        self.list_link=list_link
        
    def plot(self):
        graph_planetary_assembly= nx.union(self.planetary_gears[0].plot(),self.planetary_gears[1].plot(),rename=((self.planetary_gears[0].name+'_'), (self.planetary_gears[1].name+'_')))
        
        for i in range(len(self.planetary_gears)-2):
            graph_planetary_assembly=nx.union(graph_planetary_assembly,self.planetary_gears[i+2].plot(),rename=('',self.planetary_gears[i+2].name+'_'))
            
        for i in range(len(self.list_link)):
            graph_planetary_assembly.add_nodes_from(['EC'+str(i+1)])
            graph_planetary_assembly.add_edges_from([(self.list_link[i][0],'EC'+str(i+1)),('EC'+str(i+1),self.list_link[i][1])])
            
        plt.clf() 
        plt.cla()
        option={
            'node_size' : 100,
            'font_size' : 11,
            'font_weight':'bold'
        }
        plt.figure(figsize=(10, 10))
        nx.draw_kamada_kawai(graph_planetary_assembly, with_labels=True,**option)
        plt.show()
        
        return graph_planetary_assembly
        
        