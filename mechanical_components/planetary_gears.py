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
import time
import math as m

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
        if planetary_type=='Sun':
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
        self.nodes= nodes
        for node in self.nodes:
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
            self.nodes= nodes
            for node in self.nodes:
                self.Z_planets.append(node.Z)

        def system_equations(self):
                 matrix=npy.array([ self.Z_planets[0], self.Z_planets[1]])
                 rhs= npy.array([0])
                 return matrix,rhs
             

class Pivot(Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        self.nodes= nodes
    def system_equations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs
        
class Fixed (Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        self.nodes= nodes
    def system_equations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs
        
class Double (Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        self.nodes= nodes
    def system_equations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs
        
class Clutch (Relation):
    
    def __init__(self,name,nodes):
        Relation.__init__(self,name)
        self.nodes=nodes
    def system_equations(self):
        matrix= npy.array([1, -1])
        rhs= npy.array([0])
        return matrix,rhs

class ImposeSpeed(Relation):

    def __init__(self,name , node, input_speed):
              Relation.__init__(self,name)  
              self.input_speed=input_speed
              self.node=node
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
        self.list_nodes.extend([gearing,self.planets[0]])
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
            
    
    def test_assembly_condition(self, number_planet):
        if self.doubles:
            planet_1=self.doubles[0].nodes[0]
            planet_2=self.doubles[0].nodes[1]
            
        else:
            planet_1=self.planets[0]
            
        basic_ratio_planet_1=1
        basic_ratio_planet_2=1
        
        for i, node_1 in enumerate(self.list_nodes):
            if node_1==planet_1:
                
                inv_list_nodes=self.list_nodes[:i+1]
                inv_list_nodes=inv_list_nodes[::-1]
               
                for j,node_2 in enumerate(inv_list_nodes):
                    
                    if isinstance(node_2,GearingPlanetary) or isinstance(node_2,GearingPlanet):
                        
                        basic_ratio_planet_1=basic_ratio_planet_1*(-inv_list_nodes[j-1].Z/inv_list_nodes[j+1].Z)
                        if isinstance(inv_list_nodes[j+1],Planetary):
                            if inv_list_nodes[j+1].planetary_type=='Ring':
                            
                                basic_ratio_planet_1=basic_ratio_planet_1 *-1
                          
                for j,node_2 in enumerate(self.list_nodes[i:]):
                    
                    if isinstance(node_2,GearingPlanetary) or isinstance(node_2,GearingPlanet):
                        
                        basic_ratio_planet_2=basic_ratio_planet_2*(-self.list_nodes[i+j-1].Z/self.list_nodes[i+j+1].Z)
                        
                        if isinstance(self.list_nodes[i+j+1],Planetary):
                            if self.list_nodes[i+j+1].planetary_type=='Ring':
                            
                                basic_ratio_planet_2=basic_ratio_planet_2 *-1 
        if self.doubles:
            equation=(1/number_planet)*(1/basic_ratio_planet_1-1/basic_ratio_planet_2)*(planet_1.Z*planet_2.Z) 
        else:
            equation=(1/number_planet)*(1/basic_ratio_planet_1-1/basic_ratio_planet_2)*(planet_1.Z)

        if int(equation)==equation:                   
            return True 
                        
                        
                
                
        
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
    
    def solve(self,input_speeds_and_composants):
        system_matrix, vector_b = self.system_equations()
        n_equations = len(self.list_relations)
        n_variables= len(self.list_elements)
        
        system_matrix_solve_0=npy.zeros((2, n_variables))
        vector_b_solve_0=npy.zeros(2)
        
        system_matrix=npy.concatenate((system_matrix,system_matrix_solve_0),axis=0)
        vector_b=npy.concatenate((vector_b,vector_b_solve_0),axis=0)
        impose_speeds=[]
        for composant in input_speeds_and_composants:
            impose_speeds.append(ImposeSpeed('ImposeSpeed1',composant,input_speeds_and_composants[composant]))
        for impose_speed in impose_speeds:                                                                                                                                          
            for i,element in enumerate(self.list_elements_speed):
                if element==impose_speed.node:
                    system_matrix[n_equations][i]=impose_speed.system_equations()[0]
                    vector_b[n_equations]=impose_speed.system_equations()[1]
                    n_equations+=1
     
        
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

    def solve(self,input_speeds_and_composants):
        
        system_matrix_assembly, vector_b = self.system_equations()
        impose_speeds=[]

        
        for i,composant in enumerate(input_speeds_and_composants):
             impose_speeds.append(ImposeSpeed('impose_speed',composant,input_speeds_and_composants[composant][0]))
             matrix_composant_speed_0=npy.zeros((1,system_matrix_assembly.shape[1]))
             matrix_composant_speed_0[0][self.matrix_position([composant,input_speeds_and_composants[composant][1]])]=impose_speeds[-1].system_equations()[0]          
             system_matrix_assembly=npy.concatenate((system_matrix_assembly,matrix_composant_speed_0),axis=0)
             vector_b=npy.concatenate((vector_b,impose_speeds[-1].system_equations()[1]))
  
        for i,fixed_element in enumerate(self.list_fixed_elements):
            matrix_fixed_elements=npy.zeros((1,system_matrix_assembly.shape[1]))
            matrix_fixed_elements[0][self.matrix_position(fixed_element[0])]=self.list_fixed[i].system_equations()[0][0]
            matrix_fixed_elements[0][self.matrix_position(fixed_element[1])]=self.list_fixed[i].system_equations()[0][1]
            system_matrix_assembly=npy.concatenate((system_matrix_assembly,matrix_fixed_elements),axis=0)
            vector_b=npy.concatenate((vector_b,self.list_fixed[i].system_equations()[1]))
        

                               
            
        
        print(system_matrix_assembly) 
        print(vector_b)         
        solution = solve(system_matrix_assembly, vector_b)
        for i in range(len(self.list_elements_speed)):
            self.list_elements_speed[i].speed=solution[i]
        return solution, system_matrix_assembly
  




def test_speed_equal_speed_output_and_assembly_condition(node,solutions,planetary_1,planetary_2,planets,planet_carrier,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet):
     planetary_1.Z=node[1]+Z_range[0]
     planetary_2.Z=node[2]+Z_range[0]
     
     for i,planet in enumerate(planets):
          planet.Z=node[i+3]+Z_range[0]

     list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}        
     planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
     real_input_speed_and_composants={}
     for composant in input_speeds_and_composants:
         real_input_speed_and_composants[list_element[composant]]=input_speeds_and_composants[composant]
     if  planetary_gears_1.test_assembly_condition(number_planet):
         speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]                               
         if (speed_output_planetary_gears_1<speed_output*(1+precision/100)) and (speed_output_planetary_gears_1>speed_output*(1-precision/100)):
             solutions.append(planetary_gears_1)
             print(node) 
             return True
         else:
             return False
     else:
         return False


def test_ratio_max_ratio_min(node,planetary_1,planetary_2,planet_carrier,planets,input_speeds_and_composants,output_composant,speed_output,list_element,Z_range):
     
     planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
     
     for composant in input_speeds_and_composants:
        list_element[composant].speed=input_speeds_and_composants[composant] 
                        
     list_element[output_composant].speed=speed_output
                          
     ratio_goal=(abs(planetary_2.speed-planet_carrier.speed)/abs(planetary_1.speed-planet_carrier.speed))

     previous_composant_Z_min=node[1]+ Z_range[0]
     previous_composant_Z_max=previous_composant_Z_min
                          
     ratio_min=1
     ratio_max=1
                          
     for planet in planetary_gears_1.planets:
         gearing=True
                              
         for planets_double in planetary_gears_1.doubles:
                                  
             if planet==planets_double.nodes[1]:
                 gearing=False
                 previous_composant_Z_min=Z_range[0]
                 previous_composant_Z_max=Z_range[1]
                                      
         if gearing:    
                                  
            for i in range(Z_range[0],Z_range[1]+1):
                if m.gcd(previous_composant_Z_min,i)==1:
                    planet_Z_min=i
                    
                    ratio_max=ratio_max*previous_composant_Z_max/planet_Z_min
         
                    previous_composant_Z_min=planet_Z_min
                    break

            for i in range(Z_range[1],Z_range[0],-1):
                if m.gcd(previous_composant_Z_max,i)==1:
                    
                    planet_Z_max=i
                    ratio_min=ratio_max*previous_composant_Z_min/planet_Z_max
                    previous_composant_Z_max=planet_Z_max
                    break
                                      
     for i in range(Z_range[0],Z_range[1]+1):
        if m.gcd(node[2]+ Z_range[0],i)==1:
            ratio_min=ratio_min*i/(node[2]+ Z_range[0])

            break

     for i in range(Z_range[1],Z_range[0],-1):
        if m.gcd(node[2]+ Z_range[0],i)==1:
            ratio_max=ratio_max*i/(node[2]+ Z_range[0])

            break

     if (ratio_goal<ratio_min) or (ratio_goal>ratio_max):
         
        return False
     else:
         return True

def test_GCD_planet(nodes,planets,planetary_gears_1):

    for i,node in enumerate(nodes[3:]):
        planets[i].Z=node 
    for i,node in enumerate(nodes[4:]):
         gearing=True
                              
         for planets_double in planetary_gears_1.doubles:
                                  
             if planets[i+1]==planets_double.nodes[1]:
                 gearing=False
         if gearing:
             if m.gcd(planets[i].Z,planets[i+1].Z)!=1:

                 return False
             
    return True
             
        
def cas_vitesse_1_planetary_gears(input_speeds_and_composants,output_composant,speed_output,numbers_succesives_planet_planetary_gears,precision,number_planet,Z_range):
        debut=time.time()
        list_tree=[4]
        
        for i in range(numbers_succesives_planet_planetary_gears+2):
            list_tree.append(Z_range[1]-Z_range[0]+1)
        
        tree=dt.RegularDecisionTree(list_tree)
        
        
        if speed_output<0:
            precision=-precision
        
        solutions=[]
        inv_2=False
        inv_3=False
        inv_4=False
        while not tree.finished:

            valid=True
            ## Test GCD planetary_1 and planetary_2##
 
            if len(tree.current_node)>=4:
                
                if m.gcd(tree.current_node[1]+Z_range[0],tree.current_node[3]+Z_range[0])!=1:
                    valid= False
                    
                # for i in range(len(tree.current_node)-4):
                #     if m.gcd(tree.current_node[3+i]+7,tree.current_node[4+i]+7)!=1:
                #         valid= False
                        
                if len(tree.current_node)==numbers_succesives_planet_planetary_gears+3:
                    if m.gcd(tree.current_node[2]+Z_range[0],tree.current_node[-1]+Z_range[0])!=1:
                        valid= False
                       
            ## Planetary Gears Type 1 ##    
            if tree.current_node[0] == 0 and valid:
                    planet_carrier=PlanetCarrier('PlanetCarrier')
                    planetary_1=Planetary('Planetary_1',7,'Sun')
                    planetary_2=Planetary('Planetary_2',7,'Ring')
                    list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
                    planets=[]
                    for i in range(numbers_succesives_planet_planetary_gears):
                        planets.append(Planet('Planet'+str(i),'Simple',7))    
                    
                    ## Test Inverse Speed ##
                    if len(tree.current_node)==1:
                      real_input_speed_and_composants={}
                      planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                      
                      for composant in input_speeds_and_composants:
                          real_input_speed_and_composants[list_element[composant]]=input_speeds_and_composants[composant]
                          
                      speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]
                      
                      if (speed_output<0)!=(speed_output_planetary_gears_1<0):
                          valid=False
                    elif len(tree.current_node)==3: 
                        
                        ## Test Z ring >> last Z ring of solution ##
                        if solutions:
                            if tree.current_node[2]<(solutions[-1].planetary_2.Z+Z_range[0]):

                                valid= False
                        
                        ##Test Z ring >> Z sun ##
                        if tree.current_node[2]>tree.current_node[1] and valid:
                            
                            ## Test speed = speed output  and assembly condition## 
                            node=tree.current_node
                            for i in range(numbers_succesives_planet_planetary_gears):
                                node= node+ [Z_range[0]]
                                
                            valid=test_speed_equal_speed_output_and_assembly_condition(node,solutions,planetary_1,planetary_2,planets,planet_carrier,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet)
                            valid=False    
                                
                        else:
                            valid=False
                            
                            
            ## Planetary Gears Type 2 ## 
            elif tree.current_node[0]==1 and valid and numbers_succesives_planet_planetary_gears>=2:
                
                planet_carrier=PlanetCarrier('PlanetCarrier')
                planetary_1=Planetary('Planetary_1',7,'Sun')
                planetary_2=Planetary('Planetary_2',7,'Ring')
                list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
                planets=[]
                    
                for i in range(numbers_succesives_planet_planetary_gears):
                    planets.append(Planet('Planet'+str(i),'Double',7))                  
                
                if inv_2:
                        planets[-2].planet_type='Simple'
                        
                ## Test GCD planet##        
                if len(tree.current_node)>4:
                    
                    
                    planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                    valid=test_GCD_planet(tree.current_node,planets,planetary_gears_1)
                
                
                ## Test Inverse Speed ##   
                if len(tree.current_node)==1:
                      real_input_speed_and_composants={}
                      planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                      
                      for composant in input_speeds_and_composants:
                          real_input_speed_and_composants[list_element[composant]]=input_speeds_and_composants[composant]
                          
                      speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]
                      
                      if (speed_output<0)!=(speed_output_planetary_gears_1<0):
                          if numbers_succesives_planet_planetary_gears<4:
                             
                              valid=False
                          else:
                              inv_2=True
                              
                ## Test ratio max/ratio min for Z planetary 1 and Z planetary 2##               
                elif len(tree.current_node)==3:
                    
                    planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                    
                    valid=test_ratio_max_ratio_min(tree.current_node,planetary_1,planetary_2,planet_carrier,planets,input_speeds_and_composants,output_composant,speed_output,list_element,Z_range)
                         
                
                
                elif len(tree.current_node)==numbers_succesives_planet_planetary_gears+3:
                    
                    ## Test Z Ring >> Z planet ##
                    if (tree.current_node[2]>=tree.current_node[-1]):
                        
                        ## Test speed = speed output and assembly condition## 
                        valid=test_speed_equal_speed_output_and_assembly_condition(tree.current_node,solutions,planetary_1,planetary_2,planets,planet_carrier,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet)
                    else:
                        valid=False
                        
            ## Planetary Gears Type 3 ##                   
            elif tree.current_node[0]==2 and valid and numbers_succesives_planet_planetary_gears>=2:
                
                      planet_carrier=PlanetCarrier('PlanetCarrier')
                      planetary_1=Planetary('Planetary_1',7,'Ring')
                      planetary_2=Planetary('Planetary_2',7,'Ring')
                      list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
                      planets=[]
                    
                      for i in range(numbers_succesives_planet_planetary_gears):
                        planets.append(Planet('Planet'+str(i),'Double',7))  
                        
                      if inv_3:
                          planets[-2].planet_type='Simple' 
                          
                      ##Test GCD planet##    
                      if len(tree.current_node)>4:
                          planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                          valid=test_GCD_planet(tree.current_node,planets,planetary_gears_1)
                          
                      ## Test Inverse Speed ## 
                      if len(tree.current_node)==1:
                          real_input_speed_and_composants={}
                          
                          planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                          
                          for composant in input_speeds_and_composants:
                              real_input_speed_and_composants[list_element[composant]]=input_speeds_and_composants[composant]
                              
                          speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]
                          
                          if (speed_output<0)!=(speed_output_planetary_gears_1<0):
                             if numbers_succesives_planet_planetary_gears<4: 
                                 valid=False
                             else:
                                 inv_3=True
                      ## Test ratio max/ratio min for Z planetary 1 and Z planetary 2##           
                      elif len(tree.current_node)==3:
                          
                          valid=test_ratio_max_ratio_min(tree.current_node,planetary_1,planetary_2,planet_carrier,planets,input_speeds_and_composants,output_composant,speed_output,list_element,Z_range)       
                                 
                      ## Test Z Ring1 >> Z planet1 ##
                      elif len(tree.current_node)==4:
                          if (tree.current_node[1]<=tree.current_node[3]):
                              valid=False
                              
                               
                      elif len(tree.current_node)==numbers_succesives_planet_planetary_gears+3:
                          
                          ## Test Z Ring2 >> Z planet2 ##
                          if (tree.current_node[2]>=tree.current_node[-1]):
                              
                              ## Test speed = speed output and assembly condition ## 
                              valid=test_speed_equal_speed_output_and_assembly_condition(tree.current_node,solutions,planetary_1,planetary_2,planets,planet_carrier,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet)
                          else:
                             valid=False
            
            ## Planetary Gears Type 4 ##      
            elif tree.current_node[0]==3 and valid and numbers_succesives_planet_planetary_gears>=2:
                
                      planet_carrier=PlanetCarrier('PlanetCarrier')
                      planetary_1=Planetary('Planetary_1',7,'Sun')
                      planetary_2=Planetary('Planetary_2',7,'Sun')
                      list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
                      planets=[]
                    
                      for i in range(numbers_succesives_planet_planetary_gears):
                        planets.append(Planet('Planet'+str(i),'Double',7))  
                      
                      if inv_4:
                          planets[-2].planet_type='Simple'
                      ## Test GCD planet##        
                      if len(tree.current_node)>4:
                          planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                          valid=test_GCD_planet(tree.current_node,planets,planetary_gears_1)
                        
                        
                      ## Test Inverse Speed ##  
                      if len(tree.current_node)==1:

                          real_input_speed_and_composants={}
                          
                          for composant in input_speeds_and_composants:
                              real_input_speed_and_composants[list_element[composant]]=input_speeds_and_composants[composant]
                              
                          planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                          speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]
                      
                          if (speed_output<0)!=(speed_output_planetary_gears_1<0):
                             if numbers_succesives_planet_planetary_gears<4: 
                                 valid=False
                             else:
                                inv_4=True
                                
                      ## Test ratio max/ratio min for Z planetary 1 and Z planetary 2##
                      elif len(tree.current_node)==3:
                          
                          
                          
                          valid=test_ratio_max_ratio_min(tree.current_node,planetary_1,planetary_2,planet_carrier,planets,input_speeds_and_composants,output_composant,speed_output,list_element,Z_range)
                              
                      elif len(tree.current_node)==numbers_succesives_planet_planetary_gears+3:
                          
                          ## Test speed = speed output ##
                           valid=test_speed_equal_speed_output_and_assembly_condition(tree.current_node,solutions,planetary_1,planetary_2,planets,planet_carrier,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet)
                         
            tree.NextNode(valid)
        fin=time.time()
        print(debut-fin)
        return solutions
    
    

    
    
    
    