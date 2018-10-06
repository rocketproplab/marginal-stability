# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 20:02:18 2018

@author: Neal
"""
from engine_calcs import engine_excel_writer

MW = 0
gamma = 0
Tc = 0
pe = 0
pa = 0
pc = 0
area_ratio = 0
thrust = 0
MR = 0
contraction_ratio = 0
Lstar = 0
alphad = 15
adj_coeff = 0.5
filename = 'engine_output.xlsx'

engine_excel_writer(MW,gamma,Tc,pe,pa,pc,area_ratio,thrust,MR,\
                      contraction_ratio, Lstar, alphad, adj_coeff, \
                      filename)