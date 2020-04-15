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
import copy
import volmdlr as vm
import volmdlr.primitives3D as p3d
import volmdlr.primitives2D as p2d
import mechanical_components.meshes as meshes
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
        self.module=0
        self.d=0
        if planetary_type=='Sun':
            self.p=1
        else:
            self.p=-1
    def volume_plot(self,xy_position,z_position,module,length):
         self.module=module
         self.d=module*self.Z
         radius=self.Z*module
         x = vm.Vector3D((1,0,0))
         y = vm.Vector3D((0,1,0))
         z = vm.Vector3D((0,0,1))
         rack=meshes.Rack(0.34,module)
         meshes_1=meshes.Mesh(self.Z,radius,0.01,rack)
         Gears3D={0:meshes_1.Contour(3)}
         export=[]
         center_2=(xy_position[0],xy_position[1])
         center = vm.Point2D(center_2)
         model_trans =Gears3D[0][0].Translation(center)
         model_trans_rot = model_trans.Rotation(center, 0.1)
         Gears3D_Rotate=[model_trans_rot]
         export=[]
         for (i,center,k) in zip([Gears3D[0]],[center_2],[-1]):
                        model_export=[]
                        
                        for m in i:
                            
                            center = vm.Point2D(center)
                            model_trans = m.Translation(center)
                            model_trans_rot = model_trans.Rotation(center, k)
                            model_export.append(model_trans_rot)
                        export.append(model_export)
         Gears3D_Rotate=export
         vect_x = z_position*z 
         extrusion_vector1 = length*z
         C1=vm.Contour2D(Gears3D_Rotate[0])
         i=vm.Vector3D(vect_x)
         t1=p3d.ExtrudedProfile(vm.Vector3D(vect_x), x, y, C1, [], vm.Vector3D(extrusion_vector1))               
         return t1
         



class Planet():

    def __init__ (self,name,planet_type,Z):
        self.name=name
        self.planet_type=planet_type
        self.Z=Z
        self.speed=0
        self.module=0
        
    def volume_plot(self,xy_position,z_position,module,length):
         self.module=module
         self.d=module*self.Z
         radius=self.Z*module
         x = vm.Vector3D((1,0,0))
         y = vm.Vector3D((0,1,0))
         z = vm.Vector3D((0,0,1))
         rack=meshes.Rack(0.34,module)
         meshes_1=meshes.Mesh(self.Z,radius,0.01,rack)
         Gears3D={0:meshes_1.Contour(3)}
         export=[]
         center_2=(xy_position[0],xy_position[1])
         center = vm.Point2D(center_2)
         model_trans =Gears3D[0][0].Translation(center)
         model_trans_rot = model_trans.Rotation(center, 0.1)
         Gears3D_Rotate=[model_trans_rot]
         export=[]
         for (i,center,k) in zip([Gears3D[0]],[center_2],[-1]):
                        model_export=[]
                        
                        for m in i:
                            
                            center = vm.Point2D(center)
                            model_trans = m.Translation(center)
                            model_trans_rot = model_trans.Rotation(center, k)
                            model_export.append(model_trans_rot)
                        export.append(model_export)
         Gears3D_Rotate=export
         vect_x = z_position*z 
         extrusion_vector1 = length*z
         C1=vm.Contour2D(Gears3D_Rotate[0])
         i=vm.Vector3D(vect_x)
         t1=p3d.ExtrudedProfile(vm.Vector3D(vect_x), x, y, C1, [], vm.Vector3D(extrusion_vector1))               
         return t1
class PlanetCarrier():

    def __init__ (self,name):
        self.name=name
        self.speed=0

class GearingPlanetary(Relation):

    def __init__(self,name,nodes):

        Relation.__init__(self,name)
        self.Z_planetary=0
        self.nodes= nodes

    def system_equations(self):
        for node in self.nodes:
            if isinstance(node,Planet):
                Z_planet=node.Z
            else:
                Z_planetary= node.p*node.Z


        matrix=npy.array([Z_planetary ,Z_planet,-Z_planetary])
        rhs= npy.array([0])
        return matrix,rhs

class GearingPlanet(Relation):
        def __init__(self,name,nodes):

            Relation.__init__(self,name)
            self.nodes= nodes
 

        def system_equations(self):
                 Z_planets=[]
                 for node in self.nodes:
                     Z_planets.append(node.Z)
                 matrix=npy.array([ Z_planets[0], Z_planets[1]])
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
    def volume_plot(self,xy_position,z_position,radius,length):
         
         pos = vm.Point3D((xy_position[0],xy_position[1],z_position))
         axis = vm.Vector3D((0,0,1))
         cylinder = p3d.Cylinder(pos, axis, radius, length)
         return cylinder
         

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

    def __init__(self,name,planetary_1,planetary_2,planets,planet_carrier,numbers_planets):

        #Initialize Planetarygears elements
        self.name=name
        self.planetary_1=planetary_1
        self.planetary_2=planetary_2
        self.planets=planets
        self.planet_carrier=planet_carrier
        self.list_elements=[self.planetary_1]+self.planets+[self.planetary_2,self.planet_carrier]
        self.list_elements_speed=[self.planetary_1,self.planetary_2]+self.planets+[self.planet_carrier]
        #Initialize relations


        self.pivots=[]
        self.doubles=[]
        #Initialize gearing , pivots and Double Planets
        flag_double=0

        gearing=GearingPlanetary('Ge_1',[self.list_elements[0],self.list_elements[1]])
        self.list_nodes=[self.planetary_1,gearing,self.planets[0]] #list node is use for graph


        self.gearings=[gearing]
        self.list_relations=[gearing]

        
        i=0
        if planets[0].planet_type == 'Double':
                flag_double+=1
        list_nodes_append=self.list_nodes.append
        list_relations_append=self.list_relations.append 
        gearins_append=self.gearings.append
        doubles_append=self.doubles.append
        pivots_append=self.pivots.append
        for i,planet in enumerate(self.planets[1:]):

            if planet.planet_type == 'Double':
                flag_double+=1

            if flag_double == 2:
                double=Double('Double'+str(i+1),[self.planets[i],self.planets[i+1]])
                doubles_append(double)
                list_nodes_append(double)
                list_relations_append(double)
                flag_double=0

            else:
                gearing=GearingPlanet('Ge_'+str(i+1),[self.planets[i],self.planets[i+1]])
                gearins_append(gearing)
                list_nodes_append(gearing)
                list_relations_append(gearing)


            pivots_append(Pivot('Pv'+str(i+1),[planet,self.planet_carrier]))
            list_nodes_append(planet)

        gearing=GearingPlanetary('Ge_'+str(i+2),[self.list_elements[-3], self.list_elements[-2]])
        gearins_append(gearing)
        self.list_nodes.extend([gearing,self.planetary_2])
        list_relations_append(gearing)
    
    def __str__(self):
        Z_planets=[]
        for planet in self.planets:
            Z_planets.append(planet.Z)
        return 'Name:' + self.name + '\n\n' + \
               'Z_planetary_1:' + str(self.planetary_1.Z) + '\n' + \
               'Planetary_1_type:'+ self.planetary_1.planetary_type + '\n\n' + \
               'Z_planetary_2:' + str(self.planetary_2.Z) + '\n' + \
               'Planetary_2_type:' + str(self.planetary_2.planetary_type) + '\n\n' + \
               'Planets_Number:' + str(len(self.planets)) + '\n' + \
               'Planets_Double_Number:' + str(len(self.doubles)) + '\n' + \
               'Z_Planets:' + str(Z_planets) + '\n\n\n'
        
    def matrix_position(self,element):
                return self.list_elements_speed.index(element)


    def test_assembly_condition(self, number_planet):
        doubles_2=self.doubles
        if doubles_2:
            planet_1=doubles_2[0].nodes[0]
            planet_2=doubles_2[0].nodes[1]

        else:
            planet_1=self.planets[0]

        basic_ratio_planet_1=1
        basic_ratio_planet_2=1
        list_nodes_2=self.list_nodes
        position_planet= list_nodes_2.index(planet_1)


        inv_list_nodes=list_nodes_2[:position_planet+1]
        inv_list_nodes=inv_list_nodes[::-1]

        for j,node_2 in enumerate(inv_list_nodes):

            if isinstance(node_2,(GearingPlanetary,GearingPlanet)):

                basic_ratio_planet_1=basic_ratio_planet_1*(-inv_list_nodes[j-1].Z/inv_list_nodes[j+1].Z)
                if isinstance(inv_list_nodes[j+1],Planetary):
                    if inv_list_nodes[j+1].planetary_type=='Ring':

                        basic_ratio_planet_1=basic_ratio_planet_1 *-1

        for j,node_2 in enumerate(list_nodes_2[position_planet:]):

            if isinstance(node_2,(GearingPlanetary,GearingPlanet)):

                basic_ratio_planet_2=basic_ratio_planet_2*(-list_nodes_2[position_planet+j-1].Z/list_nodes_2[position_planet+j+1].Z)

                if isinstance(list_nodes_2[position_planet+j+1],Planetary):
                    if list_nodes_2[position_planet+j+1].planetary_type=='Ring':

                        basic_ratio_planet_2=basic_ratio_planet_2 *-1
        if doubles_2:
            equation=(1/number_planet)*(1/basic_ratio_planet_1-1/basic_ratio_planet_2)*(planet_1.Z*planet_2.Z)

        else:
            equation=(1/number_planet)*(1/basic_ratio_planet_1-1/basic_ratio_planet_2)*(planet_1.Z)

        return (int(equation)==equation)




    def plot(self):

        graph_planetary_gear= nx.Graph()

        for node, next_node in zip(self.list_nodes[:-1], self.list_nodes[1:]):
            graph_planetary_gear.add_edge(node.name, next_node.name)



        for k,planet in enumerate(self.planets):
            graph_planetary_gear.add_edge(self.planet_carrier.name,self.pivots[k].name)
            graph_planetary_gear.add_edge(self.pivots[k].name,planet.name)

        nx.draw_kamada_kawai(graph_planetary_gear, with_labels=True)

        return graph_planetary_gear
    
    
    def volume_plot(self,xy_position,z_position,modules,lenght_gear,length_double,radius_double):
        previous_xy_position=0
        primitives=[]
        n_plane_gearing=0
        next_z_position=z_position
        for node in self.list_nodes:
             if isinstance(node,Planetary):
                 next_z_position=next_z_position
                 primitives.append(node.volume_plot(copy.copy(xy_position),copy.copy(next_z_position), modules[n_plane_gearing],lenght_gear))
                
                 
                 next_xy_position=[xy_position[0],xy_position[1]+node.d/2+modules[n_plane_gearing]/2]
                 
                 
             elif isinstance(node,Planet):
                 if next_xy_position!=previous_xy_position:
                     next_xy_position[1]+=(node.Z*modules[n_plane_gearing])/2 +modules[n_plane_gearing]/2
                 primitives.append(node.volume_plot(copy.copy(next_xy_position),copy.copy(next_z_position), modules[n_plane_gearing],lenght_gear))
                 previous_xy_position=copy.copy(next_xy_position)
                 next_xy_position[1]+=node.d/2 +modules[n_plane_gearing]/2
                 next_z_position=next_z_position
                 
             elif isinstance(node,Double)   :

                 primitives.append(node.volume_plot(previous_xy_position,next_z_position+lenght_gear+(length_double)/2,radius_double ,length_double))
                 next_z_position+=lenght_gear+length_double
                 next_xy_position=previous_xy_position
                 n_plane_gearing+=1
        return primitives
             
         

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





def test_speed_equal_speed_output_and_assembly_condition(node,node_data,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet,list_element):
     planetary_gears_1.planetary_1.Z=node_data[node[:2]]
     planetary_gears_1.planetary_2.Z=node_data[node[:3]]

     for i,planet in enumerate(planetary_gears_1.planets):
          planet.Z=node_data[node[:i+4]]




     if  planetary_gears_1.test_assembly_condition(number_planet):
         key_input_speeds_and_composant=list(input_speeds_and_composants.keys())
         real_input_speed_and_composants={list_element[key_input_speeds_and_composant[0]]:input_speeds_and_composants[key_input_speeds_and_composant[0]],
                                          list_element[key_input_speeds_and_composant[1]]:input_speeds_and_composants[key_input_speeds_and_composant[1]]}
         speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]

         if (speed_output_planetary_gears_1<speed_output*(1+precision/100)) and (speed_output_planetary_gears_1>speed_output*(1-precision/100)):

             return True
     return False

def test_ratio_max_ratio_min(node,node_data,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,Z_range):



    list_element={'Planet_Carrier': planetary_gears_1.planet_carrier, 'Planetary_1' : planetary_gears_1.planetary_1,'Planetary_2':planetary_gears_1.planetary_2}
    for composant in input_speeds_and_composants:
       list_element[composant].speed=input_speeds_and_composants[composant]

    list_element[output_composant].speed=speed_output

    ratio_goal=(abs(planetary_gears_1.planetary_2.speed-planetary_gears_1.planet_carrier.speed)/abs(planetary_gears_1.planetary_1.speed-planetary_gears_1.planet_carrier.speed))

    previous_composant_Z_min=node_data[node[:2]]
    previous_composant_Z_max=previous_composant_Z_min

    ratio_min=1
    ratio_max=1

    for j,planet in enumerate(planetary_gears_1.planets):
        gearing=True

        for planets_double in planetary_gears_1.doubles:

            if planet==planets_double.nodes[1]:

                gearing=False
                previous_composant_Z_min=Z_range[1]
                previous_composant_Z_max=Z_range[0]


        if gearing:

           for i in range(Z_range[0],Z_range[1]+1):
               if m.gcd(previous_composant_Z_min,i)==1:

                   planet_Z_min=i
                   ratio_max=ratio_max*previous_composant_Z_min/planet_Z_min
                   previous_composant_Z_min=planet_Z_min
                   break

           for i in range(Z_range[1],Z_range[0],-1):
               if m.gcd(previous_composant_Z_max,i)==1:

                   planet_Z_max=i
                   ratio_min=ratio_min*previous_composant_Z_max/planet_Z_max
                   previous_composant_Z_max=planet_Z_max
                   break


    ratio_min=ratio_min*previous_composant_Z_max/(node_data[node[:3]])
    ratio_max=ratio_max*previous_composant_Z_min/(node_data[node[:3]])
    
    if ratio_min==ratio_max:
        return int(abs(1/(ratio_goal/node_data[node[:2]])))

    if (ratio_goal<ratio_min) or (ratio_goal>ratio_max):
       
       return False
       
    return True


def test_GCD_planet(nodes,node_data,planets,planetary_gears_1):

    for i in enumerate(nodes[3:]):
        planets[i].Z=node_data[nodes[:4+i]]
    for i,node in enumerate(nodes[4:]):
         gearing=True

         for planets_double in planetary_gears_1.doubles:

             if planets[i+1]==planets_double.nodes[1]:

                 gearing=False
         if gearing:
             if m.gcd(planets[i].Z,planets[i+1].Z)!=1:

                 return False

    return True


def decision_tree_planetary_gears(input_speeds_and_composants, output_composant,
                                  speed_output, numbers_succesives_planet_planetary_gears,
                                  precision, number_planet, Z_range,module_range,module_increment):
        debut=time.time()
        # list_tree=[4]

        # for i in range(numbers_succesives_planet_planetary_gears+2):
        #     list_tree.append(Z_range[1]-Z_range[0]+1)
        list_Z=[i for i in range(Z_range[0],Z_range[1]+1)]    
        tree=dt.DecisionTree()
        node=list(tree.current_node)
        tree.SetCurrentNodeNumberPossibilities(4)
        node=tree.NextNode(True)
        if speed_output<0:
            precision=-precision

        solutions=[]
        inv=False

        
        while not tree.finished:
            # if len(node)<=numbers_succesives_planet_planetary_gears+3:
            #     tree.SetCurrentNodeNumberPossibilities(Z_range[1]-Z_range[0]+1)
            
            if len(node)==0:
                tree.SetCurrentNodeNumberPossibilities(4)
                node=tree.NextNode(True)
                continue
            if len(node)==1:
                inv=False
            if len(node)==numbers_succesives_planet_planetary_gears+3:
                tree.SetCurrentNodeNumberPossibilities(0)
            
            valid=True
            planet_carrier=PlanetCarrier('PlanetCarrier')
            node_data=tree.data
            node=tuple(node)
            ## Planetary Gears Type 1 ##
            if node[0] == 0 :
                     
                planetary_1=Planetary('Planetary_1',7,'Sun')
                planetary_2=Planetary('Planetary_2',7,'Ring')
                list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
                planets=[]
                
                for i in range(numbers_succesives_planet_planetary_gears):
                    planets.append(Planet('Planet'+str(i),'Simple',7))
                    
                planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                
                if len(node)==2:
                    node_2=node+tuple([Z_range[1]])
                    node_data[node_2]=Z_range[1]
                    Z_planetary_2=test_ratio_max_ratio_min(node_2,node_data,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,Z_range)
                    
                    Z_possibilities=[Z_planetary_2-1,Z_planetary_2,Z_planetary_2+1]
                    if Z_planetary_2 and Z_planetary_2<Z_range[1]:
                        posibilities=[]
                        for k in Z_possibilities:
                            
                            ## Test Z ring >> last Z ring of solution ##
                            if solutions:
                                if k<(solutions[-1].planetary_2.Z):
                                    valid=False
                                    ## Test speed = speed output  and assembly condition##
                            if valid:
                                
                                node_2=node+tuple([k])
                                node_data[node_2]=k
                                for j in range(numbers_succesives_planet_planetary_gears):
                                    node_2= node_2+ tuple([Z_range[0]])
                                    node_data[node_2]=Z_range[0]
                                    
                                    
                                if test_speed_equal_speed_output_and_assembly_condition(node_2,node_data,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet,list_element):
                                    
                                    posibilities.append(k) 
                        tree.SetCurrentNodeDataPossibilities(posibilities)
                    else:
                        tree.SetCurrentNodeNumberPossibilities(0)
                    
                elif len(node)==numbers_succesives_planet_planetary_gears+3 and valid:
                        planetary_gears_1.planetary_1.Z=node_data[node[:2]]
                        planetary_gears_1.planetary_2.Z=node_data[node[:3]]
                        for i in range(numbers_succesives_planet_planetary_gears):
                            planetary_gears_1.planets[i].Z=node_data[node[:4+i]]

                        if (node_data[node[:2]])/2+sum(node_data[node[:4+i]] for i in range(len(node[3:])))==(node_data[node[:3]])/2:
                            
                            solutions.append(planetary_gears_1)
                            print(planetary_gears_1)
                        tree.SetCurrentNodeNumberPossibilities(0)
             
                                


                    
                    
            ## Planetary Gears Type 2 ##   
            elif node[0]==1  and numbers_succesives_planet_planetary_gears>=2:
                
                planetary_1=Planetary('Planetary_1',7,'Sun')
                planetary_2=Planetary('Planetary_2',7,'Ring')
                list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
                planets=[]
                
                for i in range(numbers_succesives_planet_planetary_gears):
                    planets.append(Planet('Planet'+str(i),'Double',7))
    
                if inv:
                        planets[-2].planet_type='Simple'
    
                planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                if len(node)==numbers_succesives_planet_planetary_gears+3:

                        ## Test speed = speed output and assembly condition##
                        valid=test_speed_equal_speed_output_and_assembly_condition(node,node_data,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet,list_element)
                        
                        ## Test module##
                        if valid and len(planetary_gears_1.doubles)<2:
                            module_calcul_Z_1= (node_data[node[:2]]+node_data[node[:4]])/2

                            module_calcul_Z_2= (node_data[node[:5]])/2 +sum(node_data[node[:6+i]] for i in enumerate(node[5:]))
                            
                            for i in npy.arange (module_range[0],module_range[1],module_increment):
                                valid=False
                                for j in npy.arange (module_range[0],module_range[1],module_increment):
                                    if (i*module_calcul_Z_1+j*module_calcul_Z_2)==(j*node_data[node[:3]])/2:

                                        valid=True
                                        break
                                if valid:
                                    planetary_gears_1.planetary_1.module= i
                                    planetary_gears_1.planets[0].module=i
                                    planetary_gears_1.planetary_2.module=j
                                    for planet in planetary_gears_1.planets[:1]:
                                        planet.module=j
                                    solutions.append(planetary_gears_1)
                                    print(planetary_gears_1)    
                                    break
                        elif valid:
                            solutions.append(planetary_gears_1)
                            print(planetary_gears_1)
                        tree.SetCurrentNodeNumberPossibilities(0)
            
            ## Planetary Gears Type 3 ##
            elif node[0]==2  and numbers_succesives_planet_planetary_gears>=2:   
                
                planetary_1=Planetary('Planetary_1',7,'Ring')
                planetary_2=Planetary('Planetary_2',7,'Ring')
                list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}

                planets=[]

                for i in range(numbers_succesives_planet_planetary_gears):
                  planets.append(Planet('Planet'+str(i),'Double',7))

                if inv:
                    planets[-2].planet_type='Simple'
                planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                
                if len(node)==numbers_succesives_planet_planetary_gears+3:


                              ## Test speed = speed output and assembly condition ##
                              
                              valid=test_speed_equal_speed_output_and_assembly_condition(node,node_data,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet,list_element)
                              ## Test module##
                              if valid and len(planetary_gears_1.doubles)<2:
                                  module_calcul_Z_1= (node_data[node[:2]]-node_data[node[:4]])/2
                                  module_calcul_Z_2= (node_data[node[:3]])/2-(node_data[node[:5]])/2 -sum(node_data[node[:6+i]] for i in enumerate(node[5:]))
                            
                                  for i in npy.arange (module_range[0],module_range[1],module_increment):
                                      valid=False
                                      for j in npy.arange (module_range[0],module_range[1],module_increment):
                                          if i*module_calcul_Z_1==j*module_calcul_Z_2:

                                              valid=True
                                              break
                                      if valid:
                                              planetary_gears_1.planetary_1.module= i
                                              planetary_gears_1.planets[0].module=i
                                              planetary_gears_1.planetary_2.module=j
                                              for planet in planetary_gears_1.planets[:1]:
                                                  planet.module=j
                                              solutions.append(planetary_gears_1)
                                              print(planetary_gears_1)
                              elif valid:
                                  solutions.append(planetary_gears_1)
                                  print(planetary_gears_1)
                              tree.SetCurrentNodeNumberPossibilities(0)
                
            ## Planetary Gears Type 4 ##
            elif node[0]==3  and numbers_succesives_planet_planetary_gears>=2: 
                
                planetary_1=Planetary('Planetary_1',7,'Sun')
                planetary_2=Planetary('Planetary_2',7,'Sun')
                list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
                planets=[]

                for i in range(numbers_succesives_planet_planetary_gears):
                  planets.append(Planet('Planet'+str(i),'Double',7))

                if inv:
                    planets[-2].planet_type='Simple'
                planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                
                if len(node)==numbers_succesives_planet_planetary_gears+3:

                          ## Test speed = speed output ##
                           valid=test_speed_equal_speed_output_and_assembly_condition(node,node_data,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet,list_element)
                           ## Test module##
                           if valid and len(planetary_gears_1.doubles)<2:
                                module_calcul_Z_1= (node_data[node[:2]]+node_data[node[:4]])/2 
                                module_calcul_Z_2= (node_data[node[:3]])/2+(node_data[node[:5]])/2 +sum(node_data[node[:6+i]] for i in enumerate(node[5:]))
                          
                                for i in npy.arange (module_range[0],module_range[1],module_increment):
                                    valid=False
                                    for j in npy.arange (module_range[0],module_range[1],module_increment):
                                        if i*module_calcul_Z_1==j*module_calcul_Z_2:
    
                                            valid=True
                                            break
                                    if valid:
                                            planetary_gears_1.planetary_1.module= i
                                            planetary_gears_1.planets[0].module=i
                                            planetary_gears_1.planetary_2.module=j
                                            for planet in planetary_gears_1.planets[:1]:
                                                planet.module=j
                                            solutions.append(planetary_gears_1)
                                            print(planetary_gears_1)
                                            if len(solutions)==10:
                                                print(len(solutions))
                                                return solutions
                           elif valid:
                                solutions.append(planetary_gears_1)
                                print(planetary_gears_1)

                           tree.SetCurrentNodeNumberPossibilities(0)
                                
                
            ## Test_Inverse_Speed##        
            if len(node)==1:
                if node[0]!=3:
                    tree.SetCurrentNodeNumberPossibilities(0)
                elif numbers_succesives_planet_planetary_gears==1 and node[0]!=0:
                    tree.SetCurrentNodeNumberPossibilities(0)
                else:    
                    key_input_speeds_and_composant=list(input_speeds_and_composants.keys())
                    real_input_speed_and_composants={list_element[key_input_speeds_and_composant[0]]:input_speeds_and_composants[key_input_speeds_and_composant[0]],
                                        list_element[key_input_speeds_and_composant[1]]:input_speeds_and_composants[key_input_speeds_and_composant[1]]}
    
                    speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]

                    if (speed_output<0)!=(speed_output_planetary_gears_1<0):
                        if numbers_succesives_planet_planetary_gears<4:
                            tree.SetCurrentNodeNumberPossibilities(0)
                            valid=False
                        else:
                            tree.SetCurrentNodeDataPossibilities(list_Z)
                            inv=True
                    else:
                        tree.SetCurrentNodeDataPossibilities(list_Z) 
                    
            ## Test GCD ##          
            elif len(node)==3:
                possibilities=[]
                if node[0]==2:
                    Z_max=node_data[node[:2]]
                else:
                    Z_max=Z_range[1]
                for i in range(Z_range[0],Z_max):
                    
                    if m.gcd(node_data[node[:2]],i)==1:
                        possibilities.append(i)
                
                tree.SetCurrentNodeDataPossibilities(possibilities)
                
            elif len(node)>3 and len(node)<numbers_succesives_planet_planetary_gears+2:
                possibilities=[]
                for i in list_Z:
                    node_2=node+tuple([i])
                    node_date[node_2]=i
                    if test_GCD_planet(node_2,node_data,planets,planetary_gears_1) :
                        possibilities.append(i)
                        
                tree.SetCurrentNodeDataPossibilities(possibilities)    
                
            if len(node)==numbers_succesives_planet_planetary_gears+2:
                possibilities=[]
                if node[0]!=3:
                    Z_max=node_data[node[:3]]
                else:
                    Z_max=Z_range[1]
                for i in range(Z_range[0],Z_max):
                    
                    if m.gcd(node_data[node[:3]],i)==1:
                        possibilities.append(i)
                
                tree.SetCurrentNodeDataPossibilities(possibilities)
                
            ##Test Ratio min max##    
            if len(node)==2 and node[0]!=0:
                possibilities=[]
                for i in list_Z:
                    node_2=node+tuple([i])
                    node_data[node_2]=i
                    if test_ratio_max_ratio_min(node_2,node_data,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,Z_range):
                        possibilities.append(i)
                tree.SetCurrentNodeDataPossibilities(possibilities) 
            
            node=tree.NextNode(valid)  
            # ## Test GCD planetary_1 and planetary_2##   
                
            # elif len(node)==3:
            #     possibilities=[]
            #     for i in range(Z_range[0],Z_range[1]+1):
                    
            #         if m.gcd(node[1]+Z_range[0],i)==1:
            #             possibilities.append(i)
                
            #     tree.SetCurrentNodeDataPossibilities(possibilities)        
                        

            # if len(node)>=4:

            #     # if m.gcd(node[1]+Z_range[0],node[3]+Z_range[0])!=1:
            #     #     valid= False
            #     #     tree.SetCurrentNodeNumberPossibilities(0)
                    


            #     if len(node)==numbers_succesives_planet_planetary_gears+3:
            #         tree.SetCurrentNodeNumberPossibilities(0)
            #         if m.gcd(node[2]+Z_range[0],node[-1]+Z_range[0])!=1:
            #             valid= False

            # if len(node)>4:



            #         valid=test_GCD_planet(node,planets,planetary_gears_1)
                        
            # ## Planetary Gears Type 1 ##
            # if valid and node[0] == 0 :
                    
            #         planet_carrier=PlanetCarrier('PlanetCarrier')
            #         planetary_1=Planetary('Planetary_1',7,'Sun')
            #         planetary_2=Planetary('Planetary_2',7,'Ring')
            #         list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
            #         planets=[]
            #         for i in range(numbers_succesives_planet_planetary_gears):
            #             planets.append(Planet('Planet'+str(i),'Simple',7))
            #         planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
                    
            #         ##Test GCD planet##
                    
            #         if len(node)==4 and numbers_succesives_planet_planetary_gears>1:
            #             tree.SetCurrentNodeNumberPossibilities(Z_range[1]-Z_range[0]+1)
            #         if len(node)>4:

            #               valid=test_GCD_planet(node,planets,planetary_gears_1)
                          
                          
            #         ## Test Inverse Speed ##      
            #         if len(node)==1:
            #           key_input_speeds_and_composant=list(input_speeds_and_composants.keys())
            #           real_input_speed_and_composants={list_element[key_input_speeds_and_composant[0]]:input_speeds_and_composants[key_input_speeds_and_composant[0]],
            #                               list_element[key_input_speeds_and_composant[1]]:input_speeds_and_composants[key_input_speeds_and_composant[1]]}

            #           speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]

            #           if (speed_output<0)!=(speed_output_planetary_gears_1<0):
            #               tree.SetCurrentNodeNumberPossibilities(0)
            #               valid=False
            #           else:
                          
            #               tree.SetCurrentNodeDataPossibilities(list_Z)
            #         elif len(node)==3:

            #             ## Test Z ring >> last Z ring of solution ##
            #             if solutions:
            #                 if node[2]<(solutions[-1].planetary_2.Z+Z_range[0]):

            #                     valid= False

            #             ##Test Z ring >> Z sun ##
            #             if valid: #and node[2]>node[1]:

            #                 ## Test speed = speed output  and assembly condition##
            #                 node_2=node
            #                 for i in range(numbers_succesives_planet_planetary_gears):
            #                     node_2= node_2+ [Z_range[0]]

            #                 valid= test_speed_equal_speed_output_and_assembly_condition(node_2,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet,list_element)
            #                 tree.SetCurrentNodeNumberPossibilities(Z_range[1]-Z_range[0]+1)

            #             else:            
            #                 valid=False
            #                 tree.SetCurrentNodeNumberPossibilities(0)

            #         elif len(node)==numbers_succesives_planet_planetary_gears+3 and valid:
            #             planetary_gears_1.planetary_1.Z=node[1]+Z_range[0]
            #             planetary_gears_1.planetary_2.Z=node[2]+Z_range[0]
            #             for i in range(numbers_succesives_planet_planetary_gears):
            #                 planetary_gears_1.planets[i].Z=node[3+i]+Z_range[0]
            #             if (node[1]+Z_range[0])/2+(sum(node_2+Z_range[0] for node_2 in node[3:]))==(node[2]+Z_range[0])/2:
            #                 solutions.append(planetary_gears_1)
            #                 print(planetary_gears_1)
                           



            # ## Planetary Gears Type 2 ##
            # elif valid and node[0]==1  and numbers_succesives_planet_planetary_gears>=2:
                
            #     if len(node)>3 and len(node)<numbers_succesives_planet_planetary_gears+2:
            #               num_planet=len(node)-4
            #               tree.SetCurrentNodeNumberPossibilities(Z_planet_max[num_planet]-Z_planet_min[num_planet]+1)
                          
            #     if len(node)==numbers_succesives_planet_planetary_gears+2:
            #         if (node[2]+Z_range[0]-Z_planet_min[-1]>0) :
            #               tree.SetCurrentNodeNumberPossibilities(node[2]+Z_range[0]-Z_planet_min[-1])
                          
            #         else:
            #             tree.SetCurrentNodeNumberPossibilities(0)
                
                    
            #     planet_carrier=PlanetCarrier('PlanetCarrier')
            #     planetary_1=Planetary('Planetary_1',7,'Sun')
            #     planetary_2=Planetary('Planetary_2',7,'Ring')
            #     list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
            #     planets=[]
                
            #     for i in range(numbers_succesives_planet_planetary_gears):
            #         planets.append(Planet('Planet'+str(i),'Double',7))

            #     if inv:
            #             planets[-2].planet_type='Simple'

            #     planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
            #     ## Test GCD planet##
            #     if len(node)>4:



            #         valid=test_GCD_planet(node,planets,planetary_gears_1)


            #     ## Test Inverse Speed ##
            #     if len(node)==1:
            #           inv=False                   
            #           key_input_speeds_and_composant=list(input_speeds_and_composants.keys())
            #           real_input_speed_and_composants={list_element[key_input_speeds_and_composant[0]]:input_speeds_and_composants[key_input_speeds_and_composant[0]],
            #                               list_element[key_input_speeds_and_composant[1]]:input_speeds_and_composants[key_input_speeds_and_composant[1]]}

            #           speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]

            #           if (speed_output<0)!=(speed_output_planetary_gears_1<0):
            #               if numbers_succesives_planet_planetary_gears<4:
            #                   tree.SetCurrentNodeNumberPossibilities(0)
            #                   valid=False
            #               else:
            #                   tree.SetCurrentNodeDataPossibilities(list_Z)
            #                   inv=True
            #           else:
            #               tree.SetCurrentNodeDataPossibilities(list_Z)

            #     ## Test ratio max/ratio min for Z planetary 1 and Z planetary 2##
            #     elif len(node)==3:



            #         valid,Z_planet_min,Z_planet_max=test_ratio_max_ratio_min(node,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,Z_range)
                    
            #         if valid:
            #                 tree.SetCurrentNodeNumberPossibilities(Z_planet_max[0]-Z_planet_min[0]+1)
 
                        


            #     elif len(node)==numbers_succesives_planet_planetary_gears+3:

            #         ## Test Z Ring >> Z planet ##
            #         # if (node[2]>=node[-1]):
            #             ## Test speed = speed output and assembly condition##
            #             valid=test_speed_equal_speed_output_and_assembly_condition(node,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet,list_element)
                        
            #             ## Test module##
            #             if valid and len(planetary_gears_1.doubles)<2:
            #                 module_calcul_Z_1= (node[1]+node[3])/2+Z_range[0]
            #                 module_calcul_Z_2= (node[4]+Z_range[0])/2 +sum(node_planet+Z_range[0] for node_planet in node[5:])
                            
            #                 for i in npy.arange (module_range[0],module_range[1],module_increment):
            #                     valid=False
            #                     for j in npy.arange (module_range[0],module_range[1],module_increment):
            #                         if (i*module_calcul_Z_1+j*module_calcul_Z_2)==(j*node[2]+Z_range[0])/2:

            #                             valid=True
            #                             break
            #                     if valid:
            #                         planetary_gears_1.planetary_1.module= i
            #                         planetary_gears_1.planets[0].module=i
            #                         planetary_gears_1.planetary_2.module=j
            #                         for planet in planetary_gears_1.planets[:1]:
            #                             planet.module=j
            #                         solutions.append(planetary_gears_1)
            #                         print(planetary_gears_1)    
            #                         break

                            
       



            #         # else:
            #         #     valid=False

            # ## Planetary Gears Type 3 ##
            # elif valid and node[0]==2  and numbers_succesives_planet_planetary_gears>=2:
            #           if len(node)>3 and len(node)<numbers_succesives_planet_planetary_gears+2:
            #               num_planet=len(node)-4
            #               tree.SetCurrentNodeNumberPossibilities(Z_planet_max[num_planet]-Z_planet_min[num_planet]+1)
                          
            #           if len(node)==numbers_succesives_planet_planetary_gears+2 :
            #               if(node[2]+Z_range[0]-Z_planet_min[-1]>0) :
            #                   tree.SetCurrentNodeNumberPossibilities(node[2]+Z_range[0]-Z_planet_min[-1])
            #               else:
            #                   tree.SetCurrentNodeNumberPossibilities(0)
                          
                          
            #           planet_carrier=PlanetCarrier('PlanetCarrier')
            #           planetary_1=Planetary('Planetary_1',7,'Ring')
            #           planetary_2=Planetary('Planetary_2',7,'Ring')
            #           list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}

            #           planets=[]

            #           for i in range(numbers_succesives_planet_planetary_gears):
            #             planets.append(Planet('Planet'+str(i),'Double',7))

            #           if inv:
            #               planets[-2].planet_type='Simple'
            #           planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
            #           ##Test GCD planet##
            #           if len(node)>4:

            #               valid=test_GCD_planet(node,planets,planetary_gears_1)

            #           ## Test Inverse Speed ##
            #           if len(node)==1:
            #               inv=False
            #               key_input_speeds_and_composant=list(input_speeds_and_composants.keys())
            #               real_input_speed_and_composants={list_element[key_input_speeds_and_composant[0]]:input_speeds_and_composants[key_input_speeds_and_composant[0]],
            #                               list_element[key_input_speeds_and_composant[1]]:input_speeds_and_composants[key_input_speeds_and_composant[1]]}

            #               speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]

            #               if (speed_output<0)!=(speed_output_planetary_gears_1<0):
            #                  if numbers_succesives_planet_planetary_gears<4:
            #                      valid=False
            #                      tree.SetCurrentNodeNumberPossibilities(0)
            #                  else:
            #                      inv=True
            #                      tree.SetCurrentNodeNumberPossibilities(Z_range[1]-Z_range[0]+1)
            #               else:
            #                   tree.SetCurrentNodeNumberPossibilities(Z_range[1]-Z_range[0]+1)
            #           ## Test ratio max/ratio min for Z planetary 1 and Z planetary 2##
            #           elif  len(node)==3:
                          
            #               valid,Z_planet_min,Z_planet_max=test_ratio_max_ratio_min(node,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,Z_range)
            #               if valid and (node[1]+Z_range[0]-Z_planet_min[0]>0):
                              
            #                 tree.SetCurrentNodeNumberPossibilities(node[1]+Z_range[0]-Z_planet_min[0])
                            
            #               else:
            #                 tree.SetCurrentNodeNumberPossibilities(0)

                          

            #           ## Test Z Ring1 >> Z planet1 ##
            #           # elif  len(node)==4:
            #           #     if (node[1]<node[3]):
            #           #         valid=False


            #           elif len(node)==numbers_succesives_planet_planetary_gears+3:

            #               # ## Test Z Ring2 >> Z planet2 ##
            #               # if (node[2]>node[-1]):

            #                   ## Test speed = speed output and assembly condition ##
                              
            #                   valid=test_speed_equal_speed_output_and_assembly_condition(node,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet,list_element)
            #                   ## Test module##
            #                   if valid and len(planetary_gears_1.doubles)<2:
            #                       module_calcul_Z_1= (node[1]-node[3])/2
            #                       module_calcul_Z_2= (node[2])/2-(node[4])/2 -sum(node_planet+Z_range[0] for node_planet in node[5:])
                            
            #                       for i in npy.arange (module_range[0],module_range[1],module_increment):
            #                           valid=False
            #                           for j in npy.arange (module_range[0],module_range[1],module_increment):
            #                               if i*module_calcul_Z_1==j*module_calcul_Z_2:

            #                                   valid=True
            #                                   break
            #                           if valid:
            #                                   planetary_gears_1.planetary_1.module= i
            #                                   planetary_gears_1.planets[0].module=i
            #                                   planetary_gears_1.planetary_2.module=j
            #                                   for planet in planetary_gears_1.planets[:1]:
            #                                       planet.module=j
            #                                   solutions.append(planetary_gears_1)
            #                                   print(planetary_gears_1)    
                                                  
            #               # else:
            #               #                         valid=False

            # ## Planetary Gears Type 4 ##
            # elif valid and node[0]==3  and numbers_succesives_planet_planetary_gears>=2:
                      
            #           if len(node)>3 and len(node)<numbers_succesives_planet_planetary_gears+3:
            #               num_planet=len(node)-4
            #               tree.SetCurrentNodeNumberPossibilities(Z_planet_max[num_planet]-Z_planet_min[num_planet]+1)
                          
            #           planet_carrier=PlanetCarrier('PlanetCarrier')
            #           planetary_1=Planetary('Planetary_1',7,'Sun')
            #           planetary_2=Planetary('Planetary_2',7,'Sun')
            #           list_element={'Planet_Carrier': planet_carrier, 'Planetary_1' : planetary_1,'Planetary_2':planetary_2}
            #           planets=[]

            #           for i in range(numbers_succesives_planet_planetary_gears):
            #             planets.append(Planet('Planet'+str(i),'Double',7))

            #           if inv:
            #               planets[-2].planet_type='Simple'
            #           planetary_gears_1=PlanetaryGears('PlanetaryGears1',planetary_1,planetary_2,planets,planet_carrier)
            #           ## Test GCD planet##
            #           if len(node)>4:

            #               valid=test_GCD_planet(node,planets,planetary_gears_1)


            #           ## Test Inverse Speed ##
            #           if len(node)==1:
            #               inv=False
            #               key_input_speeds_and_composant=list(input_speeds_and_composants.keys())
            #               real_input_speed_and_composants={list_element[key_input_speeds_and_composant[0]]:input_speeds_and_composants[key_input_speeds_and_composant[0]],
            #                               list_element[key_input_speeds_and_composant[1]]:input_speeds_and_composants[key_input_speeds_and_composant[1]]}


            #               speed_output_planetary_gears_1=planetary_gears_1.solve(real_input_speed_and_composants)[planetary_gears_1.matrix_position(list_element[output_composant])]

            #               if (speed_output<0)!=(speed_output_planetary_gears_1<0):
            #                  if numbers_succesives_planet_planetary_gears<4:
            #                      tree.SetCurrentNodeNumberPossibilities(0)
            #                      valid=False
            #                  else:
            #                     inv=True
            #                     tree.SetCurrentNodeNumberPossibilities(Z_range[1]-Z_range[0]+1)
            #               else:      
            #                     tree.SetCurrentNodeNumberPossibilities(Z_range[1]-Z_range[0]+1)
                                        
            #           ## Test ratio max/ratio min for Z planetary 1 and Z planetary 2##
            #           elif len(node)==3:


                          
            #               valid,Z_planet_min,Z_planet_max=test_ratio_max_ratio_min(node,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,Z_range)
            #               if valid:
                                
            #                     tree.SetCurrentNodeNumberPossibilities(Z_planet_max[0]-Z_planet_min[0]+1)


            #               else:
                              
            #                   tree.SetCurrentNodeNumberPossibilities(0)
                              
            #           elif len(node)==numbers_succesives_planet_planetary_gears+3:

            #               ## Test speed = speed output ##
            #                valid=test_speed_equal_speed_output_and_assembly_condition(node,planetary_gears_1,input_speeds_and_composants,output_composant,speed_output,precision,Z_range,number_planet,list_element)
            #                ## Test module##
            #                if valid and len(planetary_gears_1.doubles)<2:
            #                     module_calcul_Z_1= (node[1]+node[3])/2 +Z_range[0]
            #                     module_calcul_Z_2= (node[2])/2+(node[4])/2+ Z_range[0] +sum(node_planet+Z_range[0] for node_planet in node[5:])
                          
            #                     for i in npy.arange (module_range[0],module_range[1],module_increment):
            #                         valid=False
            #                         for j in npy.arange (module_range[0],module_range[1],module_increment):
            #                             if i*module_calcul_Z_1==j*module_calcul_Z_2:
    
            #                                 valid=True
            #                                 break
            #                         if valid:
            #                                 planetary_gears_1.planetary_1.module= i
            #                                 planetary_gears_1.planets[0].module=i
            #                                 planetary_gears_1.planetary_2.module=j
            #                                 for planet in planetary_gears_1.planets[:1]:
            #                                     planet.module=j
            #                                 solutions.append(planetary_gears_1)
            #                                 print(planetary_gears_1)
                                            

            
        fin=time.time()
        print(debut-fin)
        return solutions
        
class OptimizerPlanetaryGears():
    def __init__(self,input_speeds,diameter_cylinder,length_cylinder, position_cylinder):
        self.input_speeds=input_speeds
        self.diameter=diameter_cylinder
        self.length=length_cylinder
    
    def fact(self,n):
   
        if n<2:
            return 1
        else:
            return n*self.fact(n-1)
    
    def multiplication_possibility_assembly_planetary_gears(self,list_1,list_2):
        list_multiplication=[]
        for element_1 in list_1 :
            for element_2 in list_2:
                if type(element_1[0])==list:
                    element_multiplication=copy.copy(element_1)
                else:
                    element_multiplication=copy.copy([element_1]) 
                if type(element_2[0])==list:
                    for element_3 in element_2:
                        element_multiplication.append(element_3)
                else:
                    element_multiplication.append(element_2)
                   
                list_multiplication.append(element_multiplication)
        return list_multiplication
    
    def multiplication_possibility_speed(self,list_1,n,element_multiplication,list_multiplication):
        
        
        for i in range(len(list_1)):
            
            element_multiplication_2=copy.copy(element_multiplication)
            if  not list_1[i] in element_multiplication_2:
                element_multiplication_2.append(list_1[i])
                if n!=len(list_1)-1:
                    self.multiplication_possibility_speed(list_1,n+1,element_multiplication_2,list_multiplication)
                    

                else:
                    if not element_multiplication_2 in list_multiplication:
                        list_multiplication.append(element_multiplication_2)
                        
                        

         
        
    
    def decission_tree_(self):
        planet_carrier=PlanetCarrier('PlanetCarrier')
        planetary_1=Planetary('Planetary_1',7,'Sun')
        planetary_2=Planetary('Planetary_2',7,'Ring')
        planet_1=Planet('Planet_1',7,'Simple')
        list_1_element_fixed=[[planetary_1,planetary_1],[planetary_1,planet_carrier],[planet_carrier,planet_carrier]]
        list_2_element_fixed=[[[planetary_1,planetary_1],[planetary_2,planetary_2]],
                              [[planetary_1,planetary_1],[planet_carrier,planetary_1]],
                              [[planetary_1,planetary_1],[planet_carrier,planet_carrier]],
                              [[planetary_1,planet_carrier],[planet_carrier,planet_carrier]]]
        number_input=len(self.input_speeds)
       
        if number_input>3:
            even=(number_input%2==0)
            n=int(number_input-5)
            if n>0:
                number_planetary_gears=int(number_input/2)

                number_possibility_assembly_planetary_gears=3**(number_planetary_gears-2)*(3+even)
                list_possibility_assembly=list_1_element_fixed
                for i in range(number_planetary_gears-3):
                    list_possibility_assembly=self.multiplication_possibility_assembly_planetary_gears(list_possibility_assembly,list_1_element_fixed)
                     
                if even:
                    list_possibility_assembly=self.multiplication_possibility_assembly_planetary_gears(list_possibility_assembly,list_2_element_fixed)
                else:
                    list_possibility_assembly=self.multiplication_possibility_assembly_planetary_gears(list_possibility_assembly,list_1_element_fixed)
                  
                    
            else:
                number_planetary_gears=2
                number_possibility_assembly_planetary_gears=(3+even)
                
                if even:
                    
                    list_possibility_assembly=list_2_element_fixed
                else:
                    list_possibility_assembly=[[[planetary_1,planetary_1]],[[planetary_1,planet_carrier]],[[planet_carrier,planet_carrier]]]
        else:
            number_planetary_gears=1
            number_possibility_assembly_planetary_gears=1
            
        number_possibility_speed=self.fact(number_input)
        list_possibility_speed=[]
        self.multiplication_possibility_speed(self.input_speeds,0,[],list_possibility_speed)
        
        list_tree=[number_possibility_assembly_planetary_gears,number_possibility_speed]    
        tree=dt.RegularDecisionTree(list_tree)
        while not tree.finished:
            node=tree.current_node

            if len(node)==1:
                planetary_gears=[]

                for i in range(number_planetary_gears):

                    planetary_gears_1=PlanetaryGears('Planetary_gears_1',planetary_1,planetary_2,[planet_1],planet_carrier)
                    planetary_gears.append(planetary_gears_1)
                    
                list_element_assembly=list_possibility_assembly[node[0]]

                list_fixed_elements=[]

                for i in range (len(list_element_assembly)-1):   
                    list_fixed_elements.append([[list_element_assembly[i][0],planetary_gears[i]],[list_element_assembly[i][1],planetary_gears[i+1]]])
            
                list_fixed_elements.append([[list_element_assembly[-1][0],planetary_gears[-2]],[list_element_assembly[-1][1],planetary_gears[-1]]])
                assembly_planetary_gears=AssemblyPlanetaryGears('Assembly_planetary_gears_1', planetary_gears,list_fixed_elements)    
                   
                                                  
                
            elif len(node)==2:
                
                list_element_speed=list_possibility_speed[node[1]]
                
            tree.NextNode(True)
                
            
            

            
                
            
            
        
    
               






