# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 14:50:23 2018

@author: Neal
"""

def unit_converter_met(vectortype,value):
    if vectortype == 'force':
        return value * 4.448
    
    if vectortype == 'temperature':
        return value * 5/9
    
    if vectortype == 'pressure':
        return value * 6894.757

def unit_converter_imp(vectortype,value):    
    if vectortype == 'velocity':
        return value * 3.281

    if vectortype == 'massflow':
        return value * 2.205
    
    if vectortype == 'temperature':
        return value * 1.8
    
    if vectortype == 'pressure':
        return value/6894.757
    
    if vectortype == 'area':
        return value * 1550
    
    if vectortype == 'length':
        return value * 39.37
    
    if vectortype == 'volume':
        return value * 61023.7
