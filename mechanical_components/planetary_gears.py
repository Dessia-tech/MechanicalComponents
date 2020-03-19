#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:25:05 2020

@author: launay
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as npy
import dectree as dt
from scipy.linalg import solve


class Relation:
    
    def __init__ (self,name):
        self.name=name
        
    def _get_system_matrix(self):
        if not self._utd_system_equations:
            self._system_matrix, self._system_rhs = self.system_equations()
            self._utd_system_equations = True
        return self._system_matrix    
    system_matrix = property(_get_system_matrix)
    
    def _get_system_rhs(self):
        if not self._utd_system_equations:
            self._system_matrix, self._system_rhs = self.system_equations()
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

class GearingPlanetary(Relation):
    
    def __init__(self,name,nodes):
        
        Relation.__init__(self,name)
        self.Z_planetary=0
        for node in nodes:
            if isinstance(node,Planet):
                self.Z_planet=node.Z
            else:
                self.Z_planetary= node.p*node.Z
                
    def system_equations(self):
        
        
            
        matrix=npy.array([self.Z_planetary ,self.Z_planet,-self.Z_planetary])
        rhs= npy.array([0]) 
        return matrix,rhs
    
class GearingPlanet(Relation):
        def __init__(self,name,nodes):
        
            Relation.__init__(self,name)
            self.Z_planets=[]
            for node in nodes:
                self.Z_planets.append(node.Z)

        def system_equations(self):
                 matrix=npy.array([ self.Z_planets[0], self.Z_planets[1]])
                 rhs= npy.array([0])
                 return matrix,rhs
             

class Pivot(Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        
    def system_equations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs
        
class Fixed (Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        
    def system_equations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs
        
class Double (Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        
    def system_equations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs
        
class Clutch (Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        
    def system_equations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs

class ImposeSpeed(Relation):

    def __init__(self,name , node, input_speed):
              Relation.__init__(self,name)  
              self.input_speed=input_speed
              
    def system_equations(self):
        matrix= npy.array([1])
        rhs= npy.array([self.input_speed])
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
        self.list_elements_speed=[self.planetary_1,self.planetary_2]+self.planets+[self.planet_carrier]
        #Initialize relations
        
        self.gearings=[]
        self.pivots=[]
        self.doubles=[]
        #Initialize gearing , pivots and Double Planets
        flag_double=0
        
        
        self.list_relations=[]
        self.list_nodes=[self.planetary_1] #list node is use for graph
        
        gearing=GearingPlanetary('Ge_1',[self.list_elements[0],self.list_elements[1]])
        self.gearings.append(gearing)
        self.list_nodes.append(gearing)
        self.list_relations.append(gearing)
        
        i=0
        
        if planets[0].planet_type == 'Double':
                flag_double+=1
                
        for i,planet in enumerate(self.planets[1:]):           
            
            if planet.planet_type == 'Double':
                flag_double+=1
                
            if flag_double == 2:
                double=Double('Double'+str(i+1),[self.planets[i],self.planets[i+1]])
                self.doubles.append(double)
                self.list_nodes.append(double)
                self.list_relations.append(double)
                flag_double=0
                
            else: 
                gearing=GearingPlanet('Ge_'+str(i+1),[self.planets[i],self.planets[i+1]])
                self.gearings.append(gearing)
                self.list_nodes.append(gearing)
                self.list_relations.append(gearing)
                
                
            self.pivots.append(Pivot('Pv'+str(i+1),[planet,self.planet_carrier]))    
            self.list_nodes.append(planet)
            
        gearing=GearingPlanetary('Ge_'+str(i+2),[self.list_elements[-3], self.list_elements[-2]])
        self.gearings.append(gearing)
        self.list_nodes.extend([gearing,self.planetary_2])
        self.list_relations.append(gearing)
     
    def matrix_position(self,element):
        for k,element_speed in enumerate(self.list_elements_speed):
            if element_speed==element:
                return k
            
            
    def plot(self):
        
        graph_planetary_gear= nx.Graph()
        
        for node, next_node in zip(self.list_nodes[:-1], self.list_nodes[1:]):
            graph_planetary_gear.add_edge(node.name, next_node.name)
            
            
        
        for k,planet in enumerate(self.planets):
            graph_planetary_gear.add_edge(self.planet_carrier.name,self.pivots[k].name)
            graph_planetary_gear.add_edge(self.pivots[k].name,planet.name)
            
        nx.draw_kamada_kawai(graph_planetary_gear, with_labels=True)
        
        return graph_planetary_gear
    
    def system_equations(self):
        #initialize system matrix
        n_equations = len(self.list_relations)
        n_variables= len(self.list_elements)
        system_matrix = npy.zeros((n_equations, n_variables))
        rhs = npy.zeros(n_equations)
        num_satelite=0
        num_planetary=0
        for i,relation in enumerate(self.list_relations):
            matrix_relation, rhs_relation = relation.system_equations()
            if isinstance(relation,GearingPlanetary) :
                system_matrix[i][num_planetary]=   matrix_relation[0]         
                system_matrix[i][num_satelite+2]=   matrix_relation[1]
                system_matrix[i][n_variables-1]=   matrix_relation[2]
                    
                rhs[i]=rhs_relation[0]
                num_planetary+=1
                
                    
            else:
                system_matrix[i][num_satelite+2]=matrix_relation[0] 
                system_matrix[i][num_satelite+3]=matrix_relation[1]
                
                rhs[i]=rhs_relation[0]
                
                num_satelite+=1
                
       
        return system_matrix,rhs
    
    def solve(self,input_speed,input_composant,fixed_composant):
        system_matrix, vector_b = self.system_equations()
        n_equations = len(self.list_relations)
        n_variables= len(self.list_elements)
        
        system_matrix_solve_0=npy.zeros((2, n_variables))
        vector_b_solve_0=npy.zeros(2)
        
        system_matrix=npy.concatenate((system_matrix,system_matrix_solve_0),axis=0)
        vector_b=npy.concatenate((vector_b,vector_b_solve_0),axis=0)
        
        impose_speed=ImposeSpeed('ImposeSpeed',input_composant,input_speed)
        impose_fixed_composant= ImposeSpeed('fixed_composant',fixed_composant,0)
        for i,element in enumerate(self.list_elements_speed):
            if element==input_composant:
                system_matrix[n_equations][i]=impose_speed.system_equations()[0]
                vector_b[n_equations]=impose_speed.system_equations()[1]
                
            if element==fixed_composant:
                system_matrix[n_equations+1][i]=impose_fixed_composant.system_equations()[0]
                vector_b[n_equations+1]=impose_fixed_composant.system_equations()[1]
        
                
        solution = solve(system_matrix, vector_b)
        for i in range(len(self.list_elements_speed)):
            self.list_elements_speed[i].speed=solution[i]
        return solution
        

       

        
    
class AssemblyPlanetaryGears():
    
    def __init__(self,name,planetary_gears,list_fixed_elements):
        self.planetary_gears=planetary_gears
        self.list_fixed_elements=list_fixed_elements
        self.list_fixed=[]
        self.list_elements_speed=[]
        
        for planetary_gear in self.planetary_gears:
            self.list_elements_speed.extend([planetary_gear.planetary_1,planetary_gear.planetary_2]+planetary_gear.planets+[planetary_gear.planet_carrier])

        for i,link in enumerate(self.list_fixed_elements):
            self.list_fixed.append(Fixed('Fix'+str(i+1),[link[0][0],link[1][0]]))
    
    def matrix_position(self,element):
        num_element_planetary_tot=0
        for i,planetary_gear in enumerate(self.planetary_gears):
            
            if planetary_gear== element[1]:
                
              for k in range(len(planetary_gear.list_elements)):
                  
                  if self.list_elements_speed[num_element_planetary_tot+k]== element[0]:
                      
                      return num_element_planetary_tot+k
                  
            num_element_planetary_tot+=len(planetary_gear.list_elements)  
        
            
    def plot(self):
        
        graph_planetary_assembly= nx.union(self.planetary_gears[0].plot(),self.planetary_gears[1].plot(),rename=((self.planetary_gears[0].name+'_'), (self.planetary_gears[1].name+'_')))
        
        for planetary_gear in self.planetary_gears[2:]:
            graph_planetary_assembly=nx.union(graph_planetary_assembly,planetary_gear.plot(),rename=('',planetary_gear.name+'_'))
            
        for i,fixed_element in enumerate(self.list_fixed_elements):
            graph_planetary_assembly.add_edges_from([(fixed_element[0][1].name + '_' + fixed_element[0][0].name,self.list_fixed[i].name),
                                                     (self.list_fixed[i].name,fixed_element[1][1].name + '_' + fixed_element[1][0].name)])
            
        plt.clf() 
        plt.cla()
        option={
            'node_size' : 100,
            'font_size' : 11,
            'font_weight':'bold'
        }
        
        nx.draw_kamada_kawai(graph_planetary_assembly, with_labels=True,**option)
        
        
        return graph_planetary_assembly
        
    def system_equations(self):
        
        system_matrix_assembly,rhs_assembly=self.planetary_gears[0].system_equations()
        n_equations_assembly=system_matrix_assembly.shape[0]
        n_variables_assembly=system_matrix_assembly.shape[1]
      

        for planetary_gear in self.planetary_gears[1:]:
            matrix_planetary_gears=planetary_gear.system_equations()[0]
            
            matrix_0_assembly=npy.zeros((n_equations_assembly,len(planetary_gear.list_elements)))
            matrix_0_planetary_gears=npy.zeros((len(planetary_gear.list_relations),n_variables_assembly))
            
            matrix_concatenate_assembly=npy.concatenate((system_matrix_assembly,matrix_0_assembly), axis=1)
            matrix_concatenate_planetary_gears=npy.concatenate((matrix_0_planetary_gears,matrix_planetary_gears), axis=1)
            
            system_matrix_assembly= npy.concatenate((matrix_concatenate_assembly,matrix_concatenate_planetary_gears), axis=0)
            rhs_assembly=npy.concatenate((rhs_assembly,planetary_gear.system_equations()[1]), axis=0)
            
            
            
            
            n_equations_assembly=system_matrix_assembly.shape[0]            
            n_variables_assembly=system_matrix_assembly.shape[1]
            
        return system_matrix_assembly,rhs_assembly

    def solve(self,input_speed,input_composant,list_composant_speed_0):
        
        system_matrix_assembly, vector_b = self.system_equations()
        impose_speed=ImposeSpeed('impose_speed',input_composant,input_speed)
        impose_composants_speed_0= []
        
        for i,composant_speed_0 in enumerate(list_composant_speed_0):
             impose_composants_speed_0.append(ImposeSpeed('composant_speed_0',composant_speed_0,0))
             matrix_composant_speed_0=npy.zeros((1,system_matrix_assembly.shape[1]))
            
             matrix_composant_speed_0[0][self.matrix_position(composant_speed_0)]=impose_composants_speed_0[i].system_equations()[0]          
             system_matrix_assembly=npy.concatenate((system_matrix_assembly,matrix_composant_speed_0),axis=0)
             vector_b=npy.concatenate((vector_b,impose_composants_speed_0[i].system_equations()[1]))
             
        for i,fixed_element in enumerate(self.list_fixed_elements):
            matrix_fixed_elements=npy.zeros((1,system_matrix_assembly.shape[1]))
            matrix_fixed_elements[0][self.matrix_position(fixed_element[0])]=self.list_fixed[i].system_equations()[0][0]
            matrix_fixed_elements[0][self.matrix_position(fixed_element[1])]=self.list_fixed[i].system_equations()[0][1]
            system_matrix_assembly=npy.concatenate((system_matrix_assembly,matrix_fixed_elements),axis=0)
            vector_b=npy.concatenate((vector_b,self.list_fixed[i].system_equations()[1]))
        
        matrix_fixed_composant=npy.zeros((1,system_matrix_assembly.shape[1]))
        matrix_fixed_composant[0][self.matrix_position(input_composant)]=impose_speed.system_equations()[0]          
        system_matrix_assembly=npy.concatenate((system_matrix_assembly,matrix_fixed_composant),axis=0)
        vector_b=npy.concatenate((vector_b,impose_speed.system_equations()[1]))

                               
            
        
        print(system_matrix_assembly) 
        print(vector_b)         
        solution = solve(system_matrix_assembly, vector_b)
        for i in range(len(self.list_elements_speed)):
            self.list_elements_speed[i].speed=solution[i]
        return solution, system_matrix_assembly
            
def cas_vitesse_1planetary(input_composant,fixed_composant,speed_input,output_composant,speed_output,numbers_planet_planetary_gears):
        list_tree=[]
        for i in range(numbers_planet_planetary_gears+2):
            list_tree.append(74)
        
        tree=dt.RegularDecisionTree(list_tree)
        Planets=[]
        Planet_Carrier=PlanetCarrier('PlanetCarrier')
        
        while not tree.finished:
            if len(tree.current_node)==numbers_planet_planetary_gears+2:
                Planetary_1=Planetary('Planetary_1',tree.current_node[0]+7,'Sun')
                Planetary_2=Planetary('Planetary_2',tree.current_node[1]+7,'Ring')
                list_element={'Planet_Carrier': Planet_Carrier, 'Planetary_1' : Planetary_1,'Planetary_2':Planetary_2}
                Planets=[]
                
                for i in range(numbers_planet_planetary_gears):
                    Planets.append(Planet('Planet'+str(i),'Double',tree.current_node[i+2]+7))
                    
                PlanetaryGears1=PlanetaryGears('PlanetaryGears1',Planetary_1,Planetary_2,Planets,Planet_Carrier)
                speed_output_PlanetaryGears1=PlanetaryGears1.solve(speed_input,list_element[input_composant],list_element[fixed_composant])[PlanetaryGears1.matrix_position(list_element[output_composant])]
            
                if (speed_output_PlanetaryGears1<speed_output*1.00001) and (speed_output_PlanetaryGears1>speed_output*0.99999):
                        
                        print(tree.current_node)
            tree.NextNode(True)
        
   
    
    
    
    
    
    
    