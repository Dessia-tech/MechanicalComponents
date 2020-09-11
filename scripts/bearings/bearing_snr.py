#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 17:57:19 2018

@author: Pierrem
"""
import pandas
from xml.dom import minidom
import numpy as npy

chemin_catalogs='../../mechanical_components/catalogs/'
doc1 = minidom.parse('../../../SeaFiles/Dossiers_Clients/SNR/POC_Roulement/Doc_SNR/Catalogue_Industrie_E-Shop/Catalogue_Industrie_Standards_PRD.xml')
doc2 = minidom.parse('../../../SeaFiles/Dossiers_Clients/SNR/POC_Roulement/Doc_SNR/Catalogue_Industrie_E-Shop/Catalogue_Industrie_Paliers_PRD.xml')
doc3 = minidom.parse('../../../SeaFiles/Dossiers_Clients/SNR/POC_Roulement/Doc_SNR/Catalogue_Industrie_E-Shop/Catalogue_Expert_Tools_PRD.xml')

### Analyse des KeyLOV
dict_KeyLOV={}
for doc in [doc3, doc2, doc1]:
    for element in doc.getElementsByTagName('KeyLOV'):
        keyLOVId = element.getAttributeNode('keyLOVId').childNodes[0].data
        if keyLOVId not in dict_KeyLOV.keys():
            dict_KeyLOV[keyLOVId]={}
            for item in element.childNodes:
                if 'Name' in item.nodeName:
                    Name = element.getElementsByTagName('Name')[0].childNodes[0].data
                    dict_KeyLOV[keyLOVId]['Name'] = Name
                elif 'Value' in item.nodeName:
                    Key = item.getElementsByTagName('Key')[0].childNodes[0].data
                    Value = item.getElementsByTagName('Value')[0].childNodes[0].data
                    dict_KeyLOV[keyLOVId][Key] = Value

### Analyse des keys
dict_var={}
for doc in [doc3, doc2, doc1]:
    for element in doc.getElementsByTagName('ICSDictionaryAttribute'):
        key = element.getAttributeNode('attributeId').childNodes[0].data
        if key not in dict_var.keys():
            dict_var[key] = {}
            for item in element.childNodes:
                if 'ShortName' in item.nodeName:
                    ShortName = element.getElementsByTagName('ShortName')[0].childNodes[0].data
                    dict_var[key]['ShortName'] = ShortName
                elif 'Name' in item.nodeName:
                    Name = element.getElementsByTagName('Name')[0].childNodes[0].data
                    dict_var[key]['Name'] = Name
                elif 'Description' in item.nodeName:
                    Description = element.getElementsByTagName('Description')[0].childNodes
                    if Description != []:
                        dict_var[key]['Description'] = Description[0].data
                elif 'Comment' in item.nodeName:
                    Comment = element.getElementsByTagName('Comment')[0].childNodes
                    if Comment != []:
                        dict_var[key]['Comment'] = Comment[0].data
                elif 'Format' in item.nodeName:
                    for item_format in item.childNodes:
                        if 'ICSFormatKeyLov' in item_format.nodeName:
                            ICSFormatKeyLov = element.getElementsByTagName('ICSFormatKeyLov')[0].childNodes[0].data
                            dict_var[key]['ICSFormatKeyLov'] = ICSFormatKeyLov
                        elif 'Unit' in item_format.nodeName:
                            Unit = element.getElementsByTagName('Unit')[0].childNodes[0].data
                            dict_var[key]['Unit'] = Unit
    
compt = 0
dict_rlts = {'radial_roller_bearing':{},'radial_ball_bearing':{},
             'tapered_roller_bearing':{},'angular_ball_bearing':{},
             'spherical_ball_bearing':{}}

for doc in [doc3, doc2, doc1]:
    for element in doc.getElementsByTagName('ICSInstance'):
        ICSId = element.getElementsByTagName('ICSId')[0].childNodes[0].data
        ClassId = element.getElementsByTagName('ClassId')[0].childNodes[0].data
        UnitSystem = element.getElementsByTagName('UnitSystem')[0].childNodes[0].data
        dict_temp = {'Class':ClassId, 'Unit':UnitSystem}
        for item in element.childNodes:
            if 'Property' in item.nodeName:
                attributeId = item.getAttributeNode('attributeId').childNodes[0].data
                if attributeId in dict_var.keys():
                    name_var = dict_var[attributeId]['Name']
                    value = item.getElementsByTagName('Value')[0].childNodes[0].data
                    if 'ICSFormatKeyLov' in dict_var[attributeId].keys():
                        ICSFormatKeyLov=dict_var[attributeId]['ICSFormatKeyLov']
                        if value in dict_KeyLOV[ICSFormatKeyLov].keys():
                            val = dict_KeyLOV[ICSFormatKeyLov][value]
                        else:
                            val = value
                        dict_temp[name_var] = val
                    else:
                        dict_temp[name_var] = value
                    if name_var == 'Arborescence Catalogue Electronique Industrie':
                        ICSFormatKeyLov=dict_var[attributeId]['ICSFormatKeyLov']
                        val = dict_KeyLOV[ICSFormatKeyLov][value]
                        dict_temp['typ'] = val
                        
                        
        if 'Diamètre Extérieur (D)' in dict_temp.keys():
            dict_temp['D'] = dict_temp.pop('Diamètre Extérieur (D)')
            dict_temp['D'] = float(dict_temp['D'])*1e-3
        if 'Diamètre Intérieur (d)' in dict_temp.keys():
            dict_temp['d'] = dict_temp.pop('Diamètre Intérieur (d)')
            dict_temp['d'] = float(dict_temp['d'])*1e-3
        if 'Largeur du roulement ou de la bague intérieure (B)' in dict_temp.keys():
            dict_temp['B'] = dict_temp.pop('Largeur du roulement ou de la bague intérieure (B)')
            dict_temp['B'] = float(dict_temp['B'])*1e-3
        if 'Capacité charge dynamique (Cr)' in dict_temp.keys():
            dict_temp['Cr'] = dict_temp.pop('Capacité charge dynamique (Cr)')
            dict_temp['Cr'] = float(dict_temp['Cr'])*1e3
        if 'Capacité charge dynamique, C' in dict_temp.keys():
            dict_temp['Cr'] = dict_temp.pop('Capacité charge dynamique, C')
            dict_temp['Cr'] = float(dict_temp['Cr'])*1e3
        if 'Capacité charge statique (C0r)' in dict_temp.keys():
            dict_temp['C0r'] = dict_temp.pop('Capacité charge statique (C0r)')
            dict_temp['C0r'] = float(dict_temp['C0r'])*1e3
        if 'Capacité Charge Statique C0' in dict_temp.keys():
            dict_temp['C0r'] = dict_temp.pop('Capacité Charge Statique C0')
            dict_temp['C0r'] = float(dict_temp['C0r'])*1e3
        if 'Charge limite à la fatigue (Cu)' in dict_temp.keys():
            dict_temp['Cu'] = dict_temp.pop('Charge limite à la fatigue (Cu)')
            dict_temp['Cu'] = float(dict_temp['Cu'])*1e3
        
        if 'Temp maxi de Fonctionnement (T Max)' in dict_temp.keys():
            dict_temp['Tmax'] = dict_temp.pop('Temp maxi de Fonctionnement (T Max)')
            dict_temp['Tmax'] = float(dict_temp['Tmax'])
        if 'Temp mini de Fonctionnement (T min)' in dict_temp.keys():
            dict_temp['Tmin'] = dict_temp.pop('Temp mini de Fonctionnement (T min)')
            dict_temp['Tmin'] = float(dict_temp['Tmin'])
        if 'Temp maxi admissible (selon TTh)' in dict_temp.keys():
            dict_temp['Tmax'] = dict_temp.pop('Temp maxi admissible (selon TTh)')
            dict_temp['Tmax'] = float(dict_temp['Tmax'])
        if 'Temp mini admissible (selon TTh)' in dict_temp.keys():
            dict_temp['Tmin'] = dict_temp.pop('Temp mini admissible (selon TTh)')
            dict_temp['Tmin'] = float(dict_temp['Tmin'])
            
        if 'Longueur des rouleaux' in dict_temp.keys():
            dict_temp['Lw'] = dict_temp.pop('Longueur des rouleaux')
            dict_temp['Lw'] = float(dict_temp['Lw'])*1e-3
        if 'Longueur effective rouleau' in dict_temp.keys():
            dict_temp['Le'] = dict_temp.pop('Longueur effective rouleau')
            dict_temp['Le'] = float(dict_temp['Le'])*1e-3
        if 'Diamètre des Corps Roulants (Dw)' in dict_temp.keys():
            dict_temp['Dw'] = dict_temp.pop('Diamètre des Corps Roulants (Dw)')
            dict_temp['Dw'] = float(dict_temp['Dw'])*1e-3
        if 'Diamètre maxi des corps roulants Dw' in dict_temp.keys():
            dict_temp['Dw'] = dict_temp.pop('Diamètre maxi des corps roulants Dw')
            dict_temp['Dw'] = float(dict_temp['Dw'])*1e-3
        if 'Diamètre mini des corps roulants Bw' in dict_temp.keys():
            dict_temp['Bw'] = dict_temp.pop('Diamètre mini des corps roulants Bw')
            dict_temp['Bw'] = float(dict_temp['Bw'])*1e-3
        if 'Diamètre Primitif Dp' in dict_temp.keys():
            dict_temp['Dp'] = dict_temp.pop('Diamètre Primitif Dp')
            dict_temp['Dp'] = float(dict_temp['Dp'])*1e-3
        if 'Nombre de Corps Roulants/Rangée (Z)' in dict_temp.keys():
            dict_temp['Z'] = dict_temp.pop('Nombre de Corps Roulants/Rangée (Z)')
            dict_temp['Z'] = int(dict_temp['Z'])
        if 'Cote sous Rouleaux (F)' in dict_temp.keys():
            dict_temp['F'] = dict_temp.pop('Cote sous Rouleaux (F)')
            dict_temp['F'] = float(dict_temp['F'])*1e-3
        if 'Cote sur rouleaux (E)' in dict_temp.keys():
            dict_temp['E'] = dict_temp.pop('Cote sur rouleaux (E)')
            dict_temp['E'] = float(dict_temp['E'])*1e-3
        if 'Diamètre ext. Collet bague int. / Diamètre ext. bague HJ (d1)' in dict_temp.keys():
            dict_temp['d1'] = dict_temp.pop('Diamètre ext. Collet bague int. / Diamètre ext. bague HJ (d1)')
            dict_temp['d1'] = float(dict_temp['d1'])*1e-3
        if 'Diamètre ext. bague int. d1' in dict_temp.keys():
            dict_temp['d1'] = dict_temp.pop('Diamètre ext. bague int. d1')
            dict_temp['d1'] = float(dict_temp['d1'])*1e-3
        if 'Diamètre int. bague ext. D1' in dict_temp.keys():
            dict_temp['D1'] = dict_temp.pop('Diamètre int. bague ext. D1')
            dict_temp['D1'] = float(dict_temp['D1'])*1e-3
        if 'Diamètre int.Collet Bague ext. (D1)' in dict_temp.keys():
            dict_temp['D1'] = dict_temp.pop('Diamètre int.Collet Bague ext. (D1)')
            dict_temp['D1'] = float(dict_temp['D1'])*1e-3
        if 'Angle de contact, alpha' in dict_temp.keys():
            dict_temp['alpha'] = dict_temp.pop('Angle de contact, alpha')
            dict_temp['alpha'] = float(dict_temp['alpha'])/180*npy.pi
        if 'Demi-angle rouleaux beta' in dict_temp.keys():
            dict_temp['beta'] = dict_temp.pop('Demi-angle rouleaux beta')
            dict_temp['beta'] = float(dict_temp['beta'])/180*npy.pi
            
        if 'Diamètre mini épaulement BI (da min)' in dict_temp.keys():
            dict_temp['da_min'] = dict_temp.pop('Diamètre mini épaulement BI (da min)')
            dict_temp['da_min'] = float(dict_temp['da_min'])*1e-3
        if 'Diamètre de Passage mini BINJ et BIN dc' in dict_temp.keys():
            dict_temp['da_min'] = dict_temp.pop('Diamètre de Passage mini BINJ et BIN dc')
            dict_temp['da_min'] = float(dict_temp['da_min'])*1e-3
        if 'Diamètre mini épaulement BI (db min)' in dict_temp.keys():
            dict_temp['da_min'] = dict_temp.pop('Diamètre mini épaulement BI (db min)')
            dict_temp['da_min'] = float(dict_temp['da_min'])*1e-3
        if 'Diamètre maxi épaulement petite face da' in dict_temp.keys():
            dict_temp['da_max'] = dict_temp.pop('Diamètre maxi épaulement petite face da')
            dict_temp['da_max'] = float(dict_temp['da_max'])*1e-3
        if 'Diamètre maxi épaulement BI (da max)' in dict_temp.keys():
            dict_temp['da_max'] = dict_temp.pop('Diamètre maxi épaulement BI (da max)')
            dict_temp['da_max'] = float(dict_temp['da_max'])*1e-3
        if 'Diamètre maxi épaulement BE (Da max)' in dict_temp.keys():
            dict_temp['Da_max'] = dict_temp.pop('Diamètre maxi épaulement BE (Da max)')
            dict_temp['Da_max'] = float(dict_temp['Da_max'])*1e-3
        if 'Diamètre mini épaulement BE Da' in dict_temp.keys():
            dict_temp['Da_min'] = dict_temp.pop('Diamètre mini épaulement BE Da')
            dict_temp['Da_min'] = float(dict_temp['Da_min'])*1e-3
        if 'Diamètre mini épaulement BE Db' in dict_temp.keys():
            dict_temp['Da_min'] = dict_temp.pop('Diamètre mini épaulement BE Db')
            dict_temp['Da_min'] = float(dict_temp['Da_min'])*1e-3
        if 'Diamètre de passage maxi BEN Dc' in dict_temp.keys():
            dict_temp['Da_max'] = dict_temp.pop('Diamètre de passage maxi BEN Dc')
            dict_temp['Da_max'] = float(dict_temp['Da_max'])*1e-3
        if 'Diamètre mini épaulement petite face Db' in dict_temp.keys():
            dict_temp['da_min'] = dict_temp.pop('Diamètre mini épaulement petite face Db')
            dict_temp['da_min'] = float(dict_temp['da_min'])*1e-3
        if 'Osculation Chemin BI' in dict_temp.keys():
            dict_temp['osculation_bi'] = dict_temp.pop('Osculation Chemin BI')
            dict_temp['osculation_bi'] = float(dict_temp['osculation_bi'])*1e-3
        if 'Osculation Chemin BE' in dict_temp.keys():
            dict_temp['osculation_be'] = dict_temp.pop('Osculation Chemin BE')
            dict_temp['osculation_be'] = float(dict_temp['osculation_be'])*1e-3
        if 'Largeur de la bague extérieure (C )' in dict_temp.keys():
            dict_temp['C'] = dict_temp.pop('Largeur de la bague extérieure (C )')
            dict_temp['C'] = float(dict_temp['C'])*1e-3
        if 'Largeur totale (T)' in dict_temp.keys():
            dict_temp['T'] = dict_temp.pop('Largeur totale (T)')
            dict_temp['T'] = float(dict_temp['T'])*1e-3
        if 'Position Point Application Charge a' in dict_temp.keys():
            dict_temp['pos_a'] = dict_temp.pop('Position Point Application Charge a')
            dict_temp['pos_a'] = float(dict_temp['pos_a'])*1e-3
        if 'Dégagement mini Ca' in dict_temp.keys():
            dict_temp['Ca_min'] = dict_temp.pop('Dégagement mini Ca')
            dict_temp['Ca_min'] = float(dict_temp['Ca_min'])*1e-3
        if 'Dégagement mini Cb' in dict_temp.keys():
            dict_temp['Ca_min'] = dict_temp.pop('Dégagement mini Cb')
            dict_temp['Ca_min'] = float(dict_temp['Ca_min'])*1e-3
        
        if 'Rayon maxi de raccordement r1a' in dict_temp.keys():
            dict_temp['r1a'] = dict_temp.pop('Rayon maxi de raccordement r1a')
            dict_temp['r1a'] = float(dict_temp['r1a'])*1e-3
        if 'Rayon maxi de raccordement ra' in dict_temp.keys():
            dict_temp['ra'] = dict_temp.pop('Rayon maxi de raccordement ra')
            dict_temp['ra'] = float(dict_temp['ra'])*1e-3
            
        if 'Rayon mini de Raccordement (r1s min)' in dict_temp.keys():
            dict_temp['r1s_min'] = dict_temp.pop('Rayon mini de Raccordement (r1s min)')
            dict_temp['r1s_min'] = float(dict_temp['r1s_min'])*1e-3
        if 'Rayon mini de Raccordement (rs min)' in dict_temp.keys():
            dict_temp['rs_min'] = dict_temp.pop('Rayon mini de Raccordement (rs min)')
            dict_temp['rs_min'] = float(dict_temp['rs_min'])*1e-3
        if 'Fréquence propre Corps Roulants (60 t./min.)' in dict_temp.keys():
            dict_temp['fq_cr'] = dict_temp.pop('Fréquence propre Corps Roulants (60 t./min.)')
            dict_temp['fq_cr'] = float(dict_temp['fq_cr'])
        if 'Fréquence propre BE (60 t./min.)' in dict_temp.keys():
            dict_temp['fq_be'] = dict_temp.pop('Fréquence propre BE (60 t./min.)')
            dict_temp['fq_be'] = float(dict_temp['fq_be'])
        if 'Fréquence propre BI (60 t./min.)' in dict_temp.keys():
            dict_temp['fq_bi'] = dict_temp.pop('Fréquence propre BI (60 t./min.)')
            dict_temp['fq_bi'] = float(dict_temp['fq_bi'])
        if 'Fréquence propre Cage (60 t./min.)' in dict_temp.keys():
            dict_temp['fq_ca'] = dict_temp.pop('Fréquence propre Cage (60 t./min.)')
            dict_temp['fq_ca'] = float(dict_temp['fq_ca'])
        if 'Vitesse limite Lub huile' in dict_temp.keys():
            dict_temp['vit_max_oil'] = dict_temp.pop('Vitesse limite Lub huile')
            dict_temp['vit_max_oil'] = float(dict_temp['vit_max_oil'])
        if 'Vitesse Limite Mécanique' in dict_temp.keys():
            dict_temp['vit_max'] = dict_temp.pop('Vitesse Limite Mécanique')
            dict_temp['vit_max'] = float(dict_temp['vit_max'])
        if 'Vitesse thermique de référence' in dict_temp.keys():
            dict_temp['vit_max_ref'] = dict_temp.pop('Vitesse thermique de référence')
            dict_temp['vit_max_ref'] = float(dict_temp['vit_max_ref'])
        if 'Coef e' in dict_temp.keys():
            dict_temp['e'] = dict_temp.pop('Coef e')
            dict_temp['e'] = float(dict_temp['e'])
        if 'Coef charge axiale supérieur Y2' in dict_temp.keys():
            dict_temp['Y2'] = dict_temp.pop('Coef charge axiale supérieur Y2')
            dict_temp['Y2'] = float(dict_temp['Y2'])
        if 'Coef charge statique axiale Y0' in dict_temp.keys():
            dict_temp['Y0'] = dict_temp.pop('Coef charge statique axiale Y0')
            dict_temp['Y0'] = float(dict_temp['Y0'])
            
        if 'Vitesse limite Lub graisse' in dict_temp.keys():
            dict_temp['vit_max_grease'] = dict_temp.pop('Vitesse limite Lub graisse')
            dict_temp['vit_max_grease'] = float(dict_temp['vit_max_grease'])
        if 'Poids (sans bague HJ)' in dict_temp.keys():
            dict_temp['mass'] = dict_temp.pop('Poids (sans bague HJ)')
            dict_temp['mass'] = float(dict_temp['mass'])
        if 'Poids' in dict_temp.keys():
            dict_temp['mass'] = dict_temp.pop('Poids')
            dict_temp['mass'] = float(dict_temp['mass'])
        if '''Nb de segments d'arrêt''' in dict_temp.keys():
            dict_temp['nb_section_stop'] = dict_temp.pop('''Nb de segments d'arrêt''')
            dict_temp['nb_section_stop'] = int(dict_temp['nb_section_stop'])
        if 'Nb. de rainures sur BE' in dict_temp.keys():
            dict_temp['nb_groove_be'] = dict_temp.pop('Nb. de rainures sur BE')
            dict_temp['nb_groove_be'] = int(dict_temp['nb_groove_be'])
        if 'Nb de déflecteurs' in dict_temp.keys():
            dict_temp['nb_deflector'] = dict_temp.pop('Nb de déflecteurs')
            dict_temp['nb_deflector'] = int(dict_temp['nb_deflector'])
        if 'Nb de joints' in dict_temp.keys():
            dict_temp['nb_seal'] = dict_temp.pop('Nb de joints')
            dict_temp['nb_seal'] = int(dict_temp['nb_seal'])
        
        if 'Type de CRB' in dict_temp.keys():
            dict_temp['mounting'] = dict_temp.pop('Type de CRB')
        if 'Marque/Brand' in dict_temp.keys():
            dict_temp['brand'] = dict_temp.pop('Marque/Brand')
        if 'Libellé / Name' in dict_temp.keys():
            dict_temp['name'] = dict_temp.pop('Libellé / Name')
        if 'Class' in dict_temp.keys():
            dict_temp['class'] = dict_temp.pop('Class')
        if 'Classe de Jeu Radial' in dict_temp.keys():
            dict_temp['class_radial'] = dict_temp.pop('Classe de Jeu Radial')
        if 'Classe de précision' in dict_temp.keys():
            dict_temp['class_precision'] = dict_temp.pop('Classe de précision')
        if 'Fournisseurs NTN' in dict_temp.keys():
            dict_temp['FNR'] = dict_temp.pop('Fournisseurs NTN')
        if 'Matière Cage' in dict_temp.keys():
            dict_temp['material_cage'] = dict_temp.pop('Matière Cage')
        if 'Matière Corps Roulants' in dict_temp.keys():
            dict_temp['material_rolling'] = dict_temp.pop('Matière Corps Roulants')
        if 'Matière des Bagues' in dict_temp.keys():
            dict_temp['material_ring'] = dict_temp.pop('Matière des Bagues')    
        if '''Type d'alésage''' in dict_temp.keys():
            dict_temp['bore'] = dict_temp.pop('''Type d'alésage''')
        if 'Type de graisse' in dict_temp.keys():
            dict_temp['grease'] = dict_temp.pop('Type de graisse')
        if 'Type de Joint' in dict_temp.keys():
            dict_temp['seal'] = dict_temp.pop('Type de Joint')
        if 'Spécificité du joint' in dict_temp.keys():
            dict_temp['seal_ref'] = dict_temp.pop('Spécificité du joint')
            
        if 'Roulement avec / sans cage' in dict_temp.keys():
            dict_temp['cage'] = dict_temp.pop('Roulement avec / sans cage')
            if dict_temp['cage'] == 'Avec':
                dict_temp['cage'] = True
            else:
                dict_temp['cage'] = False
        
        if 'Unit' in dict_temp.keys():
            del dict_temp['Unit']
        if 'A prendre en compte dans le catalogue papier' in dict_temp.keys():
            del dict_temp['A prendre en compte dans le catalogue papier']
        if 'A prendre en compte dans le catalogue électronique' in dict_temp.keys():
            del dict_temp['A prendre en compte dans le catalogue électronique']
        if 'Arborescence Catalogue Electronique Industrie' in dict_temp.keys():
            del dict_temp['Arborescence Catalogue Electronique Industrie']
        if 'Référence système dimensionnel' in dict_temp.keys():
            del dict_temp['Référence système dimensionnel']
        if 'Diamètre de passage mini BINU db' in dict_temp.keys():
            del dict_temp['Diamètre de passage mini BINU db']
        if 'Coefficient A2 matière' in dict_temp.keys():
            del dict_temp['Coefficient A2 matière']
        if 'Image' in dict_temp.keys():
            del dict_temp['Image']
        if 'Code Fournisseur SNR/Supplier Code SNR' in dict_temp.keys():
            del dict_temp['Code Fournisseur SNR/Supplier Code SNR']
        if 'Profil Chemin bague intérieure' in dict_temp.keys():
            del dict_temp['Profil Chemin bague intérieure']
        if 'Profil chemin bague extérieure' in dict_temp.keys():
            del dict_temp['Profil chemin bague extérieure']
        if 'Profil Corps rouleaux' in dict_temp.keys():
            del dict_temp['Profil Corps rouleaux']
        if 'Référence Bague HJ' in dict_temp.keys():
            del dict_temp['Référence Bague HJ']
        if 'Déport bague HJ (B1)' in dict_temp.keys():
            del dict_temp['Déport bague HJ (B1)']
        if 'Largeur Totale Bague HJ (B2)' in dict_temp.keys():
            del dict_temp['Largeur Totale Bague HJ (B2)']
        if 'poids (bague HJ)' in dict_temp.keys():
            del dict_temp['poids (bague HJ)']
        if 'Rayon mini de Raccordement rNs' in dict_temp.keys():
            del dict_temp['Rayon mini de Raccordement rNs']
        if '''Réf du segment d'arrêt''' in dict_temp.keys():
            del dict_temp['''Réf du segment d'arrêt''']
        if 'Position mini Rainure a min' in dict_temp.keys():
            del dict_temp['Position mini Rainure a min']
        if 'Position maxi rainure a max' in dict_temp.keys():
            del dict_temp['Position maxi rainure a max']
        if 'Position mini segment Ca min' in dict_temp.keys():
            del dict_temp['Position mini segment Ca min']
        if 'Position maxi segment Ca max' in dict_temp.keys():
            del dict_temp['Position maxi segment Ca max']
        if 'Diamètre maxi fond de gorge rainure D3' in dict_temp.keys():
            del dict_temp['Diamètre maxi fond de gorge rainure D3']
        if 'Largeur mini rainure b' in dict_temp.keys():
            del dict_temp['Largeur mini rainure b']
        if 'Largeur maxi rainure b' in dict_temp.keys():
            del dict_temp['Largeur maxi rainure b']
        if 'Rayon maxi fond de rainure r0' in dict_temp.keys():
            del dict_temp['Rayon maxi fond de rainure r0']
        if '''Diamètre extérieur maxi anneau d'arrêt assemblé D4''' in dict_temp.keys():
            del dict_temp['''Diamètre extérieur maxi anneau d'arrêt assemblé D4''']
        if '''Epaisseur anneau d'arrêt f''' in dict_temp.keys():
            del dict_temp['''Epaisseur anneau d'arrêt f''']
        if 'Rayon maxi de raccordement côté ségment rNa' in dict_temp.keys():
            del dict_temp['Rayon maxi de raccordement côté ségment rNa']
        if '''Diamètre mini logement segment d'arrêt DbN''' in dict_temp.keys():
            del dict_temp['''Diamètre mini logement segment d'arrêt DbN''']
        if 'Capacité charge dynamique ULTAGE (Cultage)' in dict_temp.keys():
            del dict_temp['Capacité charge dynamique ULTAGE (Cultage)']
        if 'Coefficient f0' in dict_temp.keys():
            del dict_temp['Coefficient f0']
        if 'Rayon maxi de raccordement arbre & logement ra' in dict_temp.keys():
            del dict_temp['Rayon maxi de raccordement arbre & logement ra']
        if '''Diamètre mini logement segment d'arrêt Db''' in dict_temp.keys():
            del dict_temp['''Diamètre mini logement segment d'arrêt Db''']
        if 'Référence du manchon associé' in dict_temp.keys():
            del dict_temp['Référence du manchon associé']
        if 'Poids de graisse' in dict_temp.keys():
            del dict_temp['Poids de graisse']
        if 'R Base' in dict_temp.keys():
            del dict_temp['R Base']
        if 'Diamètre épaulement BI da (TRB Inch)' in dict_temp.keys():
            del dict_temp['Diamètre épaulement BI da (TRB Inch)']
        if 'Diamètre épaulement BI db (TRB Inch)' in dict_temp.keys():
            del dict_temp['Diamètre épaulement BI db (TRB Inch)']
        if 'Diamétre épaulement BE Da (TRB Inch)' in dict_temp.keys():
            del dict_temp['Diamétre épaulement BE Da (TRB Inch)']
        if 'Diamètre épaulement BE Db (TRB Inch)' in dict_temp.keys():
            del dict_temp['Diamètre épaulement BE Db (TRB Inch)']
        if 'Libellé ISO355' in dict_temp.keys():
            del dict_temp['Libellé ISO355']
        if 'Débordement maxi cage a2 (droite)' in dict_temp.keys():
            del dict_temp['Débordement maxi cage a2 (droite)']
        if 'Débordement maxi cage a1 (gauche)' in dict_temp.keys():
            del dict_temp['Débordement maxi cage a1 (gauche)']
        if 'Diamètre extérieur sphérique' in dict_temp.keys():
            del dict_temp['Diamètre extérieur sphérique']
        
        if 'typ' in dict_temp.keys():
            if 'Roulements à rouleaux cylindriques à 1 rangée' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=1
                dict_rlts['radial_roller_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à rouleaux cylindriques à 2 rangées' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=2
                dict_rlts['radial_roller_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à rouleaux cylindriques à 3 rangées' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=3
                dict_rlts['radial_roller_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à rouleaux cylindriques à 4 rangées' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=4
                dict_rlts['radial_roller_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à billes à 1 rangée' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=1
                dict_rlts['radial_ball_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à billes' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=1
                dict_rlts['radial_ball_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à Billes Radiaux' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=1
                dict_rlts['radial_ball_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à billes à contact oblique à 1 rangée ou appairables' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=1
                dict_rlts['angular_ball_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à Billes à Contact Oblique Haute Précision' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=1
                dict_rlts['angular_ball_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à Billes à Contact Oblique' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=1
                dict_rlts['angular_ball_bearing'][ICSId] = dict_temp
                compt += 1 
        if 'typ' in dict_temp.keys():
            if 'Roulements à billes à 4 points de contact' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=1
                dict_rlts['radial_ball_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à rotule sur billes' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=1
                dict_rlts['spherical_ball_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à billes à contact oblique à 2 rangées' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=2
                dict_rlts['angular_ball_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à billes à 2 rangées' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=2
                dict_rlts['radial_ball_bearing'][ICSId] = dict_temp
                compt += 1
        if 'typ' in dict_temp.keys():
            if 'Roulements à rouleaux coniques à 1 rangée' in dict_temp['typ']:
                del dict_temp['typ']
                dict_temp['i']=1
                dict_rlts['tapered_roller_bearing'][ICSId] = dict_temp
                compt += 1

list_key = ['typ_bearing','name_bearing']
for typ_bearing,dict_bearings in dict_rlts.items():
    for name_bearing,dict_bearing in dict_bearings.items():
        for key,val in dict_bearing.items():
            if key not in list_key:
                list_key.append(key)

pd_liste = []           
for typ_bearing,dict_bearings in dict_rlts.items():
    for name_bearing,dict_bearing in dict_bearings.items():
        liste = [typ_bearing,name_bearing]
        for key in list_key:
            if (key not in ['typ_bearing','name_bearing']) and (key in dict_bearing.keys()):
                liste.append(dict_bearing[key])
            elif key not in ['typ_bearing','name_bearing']:
                liste.append(npy.nan)
        pd_liste.append(liste)
        
dict_rlts_pandas = pandas.DataFrame(pd_liste, columns=list_key)
dict_rlts_pandas.to_csv('tableau_rlts_SNR.csv',index=False)
    
    