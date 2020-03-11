#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:25:05 2020

@author: launay
"""
import networkx as nx
import matplotlib.pyplot as plt
class Planetary():
    def __init__ (self,name,Z,planetary_type):
        self.name=name
        self.Z=Z
        self.planetary_type=planetary_type
        self.p=0
        if planetary_type=='sun':
            self.p=2
        else:
            self.p=1
        
      
        
class Planet():
    def __init__ (self,name,planet_type,Z):
        self.name=name
        self.planet_type=planet_type
        self.Z=Z
    
class PlanetCarrier():
    def __init__ (self,name):
        self.name=name
        

class PlanetaryGears():
    def __init__(self,name,planetary_1,planetary_2,planets,planet_carrier,):
        self.name=name
        self.planetary_1=planetary_1
        self.planetary_2=planetary_2
        self.planets=planets
        self.planet_carrier=planet_carrier
        self.planetaries= [planetary_1,planetary_2]
        
    def plot(self):
        graph_planetary_gear= nx.Graph()
        list_graph_1= ['BT1','Pv1',self.planetary_1.name]
        list_graph_2= ['BT3','Pv3',self.planet_carrier.name]
        graph_planetary_gear.add_nodes_from([self.planetary_1.name],object_type = self.planetary_1, Z=self.planetary_1.Z, speed=0, p=self.planetary_1.p)
        graph_planetary_gear.add_nodes_from([self.planetary_2.name],Z=self.planetary_2.Z, speed=0,p=self.planetary_2.p)
        graph_planetary_gear.add_nodes_from([self.planet_carrier.name], speed=0)
        flag_double=0
        for i,planet in enumerate(self.planets):
            
            if planet.planet_type == 'Double':
                flag_double+=1
                
            if flag_double == 2:
                list_graph_1.extend(('Double'+str(i+1),planet.name))
                flag_double=0
            else:    
                list_graph_1.extend(('E'+str(i+1),planet.name))
                graph_planetary_gear.add_nodes_from(['E'+str(i+1)],gearing='YES')

            graph_planetary_gear.add_nodes_from([planet.name], Z=planet.Z, planet_type=planet.planet_type, speed=0,speed_planetary_carrier=0)
            graph_planetary_gear.add_edges_from([(planet.name,'Pv'+str(i+4)),('Pv'+str(i+4),self.planet_carrier.name)])  
            
            
        list_graph_1.extend(('E'+str(i+2),self.planetary_2.name,'Pv2','BT2'))
        graph_planetary_gear.add_nodes_from(['E'+str(i+2)],gearing='YES')
        for k in range((len(list_graph_1)-1)): 
            graph_planetary_gear.add_edge(list_graph_1[k],list_graph_1[k+1])
            
        for k in range((len(list_graph_2)-1)):
            graph_planetary_gear.add_edge(list_graph_2[k],list_graph_2[k+1]) 
        
        
        nx.draw_kamada_kawai(graph_planetary_gear, with_labels=True)
        return graph_planetary_gear
    
    # def ratio(self):
    #     ratio_graph=self.plot()
    #     ratio_graph.remove_node(self.planet_carrier.name)
    #     list_gear=list(nx.shortest_path(ratio_graph,self.planetary_1.name,self.planetary_2.name))
    #     R=1
    #     Z=nx.get_node_attributes(self.plot(),'Z')
    #     list_gearing=list(nx.get_node_attributes(self.plot(),'gearing'))
    #     gearing=nx.get_node_attributes(self.plot(),'gearing')
    #     # for i in range(len(list_gear)):
    #     #         if :
    #     #             R=R* -(Z(list_gear[i-1])/Z(list_gear[i+1]))
    #     return R
                
    # def speed(self,input_composant,output_composant,input_speed,block_composant):
    #     graph_planetary_gear= self.plot()
    #     return 0
        

        
    
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
        
        