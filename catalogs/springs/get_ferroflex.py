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

if not os.path.isdir('ferroflex'):
    os.mkdir('ferroflex')

# =============================================================================
#   Getting list of products
# =============================================================================
#for i in range(3):
#    print(i)
#    sc=404
#    while sc!=200:
#        r=requests.get('https://www.ferroflex.fr/fr/produits/ressorts_de_compression.html?p={}'.format(i+1))
#        time.sleep(0.2+0.4*random.random())
#        sc=r.status_code
#    with open('ferroflex/ressorts_de_compression-{}.html'.format(i+1),'w') as f:
#        f.write(r.text)
#    
#    time.sleep(0.2+0.4*random.random())
    
# =============================================================================
#  Reading list of products
# =============================================================================

links={}
for fname in os.listdir('ferroflex'):
    with open('ferroflex/{}'.format(fname),'r') as f:
        s=f.read()
#        print(s)
    soup = BeautifulSoup(s)
    for article in soup.find_all('article'):
        a=article.find('a')
#        href=a.find('href')
        href=a.attrs['href']
        spring=href.split('/')[-1]
        spring=spring.replace('.html','')
        links[spring]=href
        
        