# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 20:02:18 2018

@author: Neal
"""
from engine_calcs import engine_excel_writer

MW = 21.15
gamma = 1.224
Tc = 6116.67
pe = 15
pa = 15
pc = 400
area_ratio = 5.5
thrust = 5000
MR = 2.23
contraction_ratio = 3.7946
Lstar = 45
alphad = 15
betad = 60
adj_coeff = 0.8
filename = 'TSD_1.xlsx'

engine_excel_writer(MW,gamma,Tc,pe,pa,pc,area_ratio,thrust,MR,\
                      contraction_ratio, Lstar, alphad,betad, adj_coeff, \
                      filename)