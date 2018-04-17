#!/usr/bin/env python
# -*- coding: utf-8 -*-

from bs4 import BeautifulSoup
import requests
import time
import random
import os

from bearing_skf import *

chemin='Data_SKF/'
if not os.path.isdir(chemin):
    os.mkdir(chemin)

# =============================================================
# Extraction de la base html
# =============================================================
    
soup = BeautifulSoup(html_doc)   
liste_link=[]
for p in soup.find_all('a'):
    temp=p.get("href")
    try:
        if '''http://www.toroidalrollerbearing.org/fr/products/bearings-units-housings/roller-bearings/cylindrical-roller-bearings/single-row-cylindrical-roller-bearings/single-row/''' in temp:
            if not temp in liste_link:
                liste_link.append(temp)
    except:
        pass
for i,link in enumerate(liste_link):
    name=link.split('=')[-1]
    name=name.replace(' ','_')
    name=name.replace('/','_')
    if not os.path.exists(chemin+name+'.html'):
        success=False
        while success==False:
            try:
                r=requests.get(link)
                success=True
            except requests.ConnectionError:
                time.sleep(10)
                print('wait')
        if r.status_code==200:
            
            fichier=open(chemin+name+'.html','w')
            fichier.write(r.text)
            fichier.close()
        else:
            print(r.status_code)
        print(i,name)
        time.sleep(0.1+0.5*random.random())
    else:
        print('exist')
        
# =============================================================
# Extraction de la base brut et retraitement
# =============================================================
    
liste_data={}
for fich in os.listdir(chemin):
    temp_data={}
    if fich.split('.')[-1]=='html':
        name_rlts=fich.split('.')[0]
        with open(chemin+fich,'r') as fichier:
            html_doc=fichier.read()
            
        soup = BeautifulSoup(html_doc)   
        for p in soup.find_all('table'):
            for elem in p.find_all('tr'):
                drap=0
        #        del name,data,unit
                for i,el in enumerate(elem.find_all('td')):
                    try:
                        if i==0:
                            name=el.get_text()
                        if i==1:
                            temp=el.get_text()
                            if temp in ['min.','max.']:
                                name=name+'_'+temp[0:-1]
                                name=name.replace(',','')
                        if i==2:
                            data=float(el.string)
                        if i==3:
                            unit=el.string
                    except ValueError:
                        pass
                    try:
                        if i==0:
                            name=el.get_text()
                            name=name.replace('é','e').replace(':','').replace(' ','_').replace('\'','')
                        if i==3:
                            try:
                                data=float(el.string)
                            except:
                                temp=el.find_all('a')
                                data=temp[0].get_text()
                        if i==4:
                            unit=el.string
                    except:
                        pass
                name=name.replace('Charge_dynamique_de_base','C')
                name=name.replace('Charge_statique_de_base','C0')
                name=name.replace('Vitesse_limite','Vmax')
                name=name.replace('Vitesse_de_reference','Vref')
                name=name.replace('Masse_du_roulement','Mass')
                name=name.replace('Limite_de_fatigue','Pu')
                name=name.replace('Vitesse_de_reference','Vref')
                if name=='Coefficient_de_calcul':
                    name='kr'
                name=name.replace('Calcul_de_la_charge__Coefficient_de_calcul','Y')
                name=name.replace('Calcul_de_la_charge__Valeur_limite','e')
                if name=='F' or name=='E':
                    name='EF'
                temp_data[name]=data
        liste_data[name_rlts]=temp_data

# =============================================================
# Export de la base au format CSV
# =============================================================
        
tableau_rlts={}
for key,item in liste_data.items():
    if not 'name' in tableau_rlts.keys():
        tableau_rlts['name']=[] 
    for k1,it in item.items():
        if not k1 in tableau_rlts.keys():
            tableau_rlts[k1]=[]
for key,item in liste_data.items():
    tableau_rlts['name'].append(key)
    for k1 in tableau_rlts.keys():
        if not k1=='name':
            try:
                tableau_rlts[k1].append(item[k1])
            except KeyError:
                tableau_rlts[k1].append(0)
    
import pandas
tableau_pandas=pandas.DataFrame(tableau_rlts)

tableau_pandas[['B','D','D1','D3','Da_max', 'Da_min', 'Db_min', 'EF','b', 'b1', 'd', 'd1', 'da_max', 'da_min', 'db_min','r12_min', 'r34_min', 'ra_max', 'rb_max', 's_max']]=tableau_pandas[['B','D','D1','D3','Da_max', 'Da_min', 'Db_min', 'EF','b', 'b1', 'd', 'd1', 'da_max', 'da_min', 'db_min','r12_min', 'r34_min', 'ra_max', 'rb_max', 's_max']]*1e-3
tableau_pandas[['C','C0','Pu']]=tableau_pandas[['C','C0','Pu']]*1e3
tableau_pandas[['Vmax','Vref']]=tableau_pandas[['Vmax','Vref']]*2*npy.pi/60

tableau_pandas.to_csv('tableau_rlts_SKF_brut.csv')
    