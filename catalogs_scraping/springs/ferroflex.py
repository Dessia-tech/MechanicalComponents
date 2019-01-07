#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:32:58 2018

@author: Steven Masfaraud masfaraud@dessia.tech
"""

import requests
import os
import time
import random
from bs4 import BeautifulSoup
import pandas as pd
import numpy as npy

if not os.path.isdir('ferroflex'):
    os.mkdir('ferroflex')
if not os.path.isdir('ferroflex/pages'):
    os.mkdir('ferroflex/pages')
if not os.path.isdir('ferroflex/products'):
    os.mkdir('ferroflex/products')
    
# =============================================================================
#   Getting list of products
# =============================================================================
for i in range(454):
    sc=404
    if not os.path.isfile('ferroflex/pages/ressorts_de_compression-{}.html'.format(i+1)):
        print(i)
        while sc!=200:
            r=requests.get('https://www.ferroflex.fr/fr/produits/ressorts_de_compression.html?p={}'.format(i+1))
            time.sleep(0.2+0.4*random.random())
            sc=r.status_code
        with open('ferroflex/pages/ressorts_de_compression-{}.html'.format(i+1),'w') as f:
            f.write(r.text)
        
        time.sleep(0.2+0.4*random.random())
    
# =============================================================================
#  Reading list of products
# =============================================================================

links={}
for fname in os.listdir('ferroflex/pages'):
    with open('ferroflex/pages/{}'.format(fname),'r') as f:
        s=f.read()
    soup = BeautifulSoup(s)
    for article in soup.find_all('article'):
        a=article.find('a')
        href=a.attrs['href']
        spring=href.split('/')[-1]
        spring=spring.replace('.html','')
        links[spring]=href
print(len(links.keys()))
        
# =============================================================================
# Getting each product params       
# =============================================================================
for i, (name, part_link) in enumerate(links.items()):
    if not os.path.isfile('ferroflex/products/'+name):
        print(i, ', ', name)
        full_link = 'https://www.ferroflex.fr' + part_link
        sc = 404
        while sc !=200:
            r = requests.get(full_link)
            time.sleep(0.2+0.4*random.random())
            sc = r.status_code
        with open('ferroflex/products/'+name, 'w') as f:
            f.write(r.text)
        time.sleep(0.2+0.4*random.random())
    
# =============================================================================
# Reading each product
# =============================================================================
list_parameters = ['d', 'D', 'L0', 'Mat', 'n', 'S', 'R',
                   'Fndyn', 'Lndyn', 'shdyn', 'Gew', 'PG', 'Fdn']
headers = ['Quantité', 'Prix unitaire [EUR]']
products = {'name' : [], 'prices' : []}
for key in list_parameters:
    products[key] = []
for fname in os.listdir('ferroflex/products'):
    print(fname)
    with open('ferroflex/products/{}'.format(fname), 'r') as f:
        s = f.read()
    soup = BeautifulSoup(s)
    td = soup.find(id='technical_data')
    p = soup.find(id='price')
    products['name'].append(fname)
    product_corr = {}
    for span in td.find_all('span'):
        key = span.next
        if key in list_parameters:
            strValue = span.next_sibling.next_sibling.next
            units = span.next_sibling.next_sibling.next_sibling.next_sibling.next
            try:
                value = float(strValue.replace(',', '.'))
                if units == 'mm':
                     value = value/1000
                elif units == 'N/mm':
                     value = value*1000
            except ValueError:
                value = strValue
                
            if key in ['Fndyn', 'shdyn', 'Lndyn'] and value == '\n':
                value = 0
            
            if value != '-':
                products[key].append(value)
            
    price = []
    for span in p.find_all('span')[2::2]:
        n = int(span.next)
        strValue = span.next_sibling.next_sibling.next[:4]
        value = strValue.replace(',', '.')
        price.append((n, float(value)))
            
    products['prices'].append(dict(price))
    
products_dataframe = pd.DataFrame(products)
products_dataframe.to_csv('ferroflex/catalog_SI')