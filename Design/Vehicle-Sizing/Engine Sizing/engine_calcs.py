# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 20:30:33 2018

@author: Neal
"""
import numpy as np
import xlsxwriter
from unit_conversion import unit_converter_met as ucm
from unit_conversion import unit_converter_imp as uci

Rbar = 8314 #universal gas constant - J/(kmol*K)
g = 9.81 #gravity - meters/seconds^2

###############################################################################

def imp2metric(thrust,Tc,pc,pa,pe,Lstar):
    
    met_thrust = ucm('force',thrust)
    met_Tc = ucm('temperature',Tc)
    met_pc = ucm('pressure',pc)
    met_pa = ucm('pressure',pa)
    met_pe = ucm('pressure',pe)
    met_Lstar = ucm('length',Lstar)
    
    return met_thrust, met_Tc, met_pc, met_pa, met_pe, met_Lstar

###############################################################################
    
def metric2imp(cstar,mdot,mdot_lox,mdot_fuel,Tthroat,pthroat, \
                Athroat, rthroat, Achamber, rchamber, Aexit, rexit, rn, rb,
                Vol_chamber, l_nozzle, l_cyl):  

    imp_cstar = uci('velocity',cstar)
    imp_mdot = uci('massflow',mdot)
    imp_mdot_lox = uci('massflow',mdot_lox)
    imp_mdot_fuel = uci('massflow',mdot_fuel)
    imp_Tthroat = uci('temperature',Tthroat)
    imp_pthroat = uci('pressure',pthroat)
    imp_athroat = uci('area',Athroat)
    imp_rthroat = uci('length',rthroat)
    imp_achamber = uci('area',Achamber)
    imp_rchamber = uci('length',rchamber)
    imp_aexit = uci('area',Aexit)
    imp_rexit = uci('length',rexit)
    imp_rn = uci('length',rn)
    imp_rb = uci('length',rb)
    imp_vc = uci('volume',Vol_chamber)
    imp_lnozzle = uci('length',l_nozzle)
    imp_lcyl = uci('length',l_cyl)

    return imp_cstar, imp_mdot, imp_mdot_lox, imp_mdot_fuel, imp_Tthroat, \
           imp_pthroat, imp_athroat, imp_rthroat, imp_achamber, \
           imp_rchamber, imp_aexit, imp_rexit, imp_rn, imp_rb, imp_vc, \
           imp_lnozzle, imp_lcyl
           
###############################################################################
    
def constants(MW,gamma,Tc):
    
    R = Rbar/MW
    
    cstar = np.sqrt((R * Tc/gamma)*((gamma+1)/2)**((gamma+1)/(gamma-1)))
    
    return R, cstar

###############################################################################
    
def performance(gamma,pe,pc,pa,area_ratio,R,cstar):
    Ct = np.sqrt(((2*gamma**2)/(gamma-1))*(2/(gamma+1))**((gamma+1)/(gamma-1))* \
                 (1 - (pe/pc)**((gamma-1)/gamma))) + (pe-pa)/pc * area_ratio
                 
    Isp = Ct * cstar/g
    
    return Ct, Isp

###############################################################################
    
def massflow(thrust,Isp,MR):
    
    massflow_total = thrust/(Isp*g)            
    
    massflow_oxidizer = MR/(MR+1)*massflow_total
    
    massflow_fuel = 1/(MR+1)*massflow_total
    
    return massflow_total, massflow_oxidizer, massflow_fuel

###############################################################################
    
def throat_properties(Tc,pc,gamma):
    
    Tthroat = Tc/(1 + (gamma-1)/2)
    
    pthroat = pc * (2/(gamma+1))**(gamma/(gamma-1))
    
    return Tthroat, pthroat

###############################################################################
    
#def exit_properties(gamma,R,Tc,pc,pa,pe):
    
#    uexhaust = np.sqrt(((2*gamma*R*Tc)/(gamma-1))* (1-(pe/pc)**((gamma-1)/gamma)))
    
#    asound = np.sqrt(gamma*R*Tc)
    
#    Ma_exit = uexhaust/asound
    
#    return uexhaust, Ma_exit

###############################################################################
    
def throat_geometry(cstar,mdot,pc):
    
    Athroat = cstar * mdot / pc
    
    rthroat = np.sqrt(Athroat/np.pi)
    
    dthroat = 2 * rthroat

    return Athroat, rthroat, dthroat

###############################################################################
    
def chamber_geometry(Athroat,contraction_ratio):
    Achamber = Athroat * contraction_ratio
    
    rchamber = np.sqrt(Achamber/np.pi)
    
    dchamber = 2 * rchamber

    return Achamber, rchamber, dchamber

###############################################################################
    
def exit_geometry(gamma,Athroat,pe,pc):
    
    expansion_ratio = (2/(gamma+1))**(1/(gamma-1)) * (pc/pe)**(1/gamma) \
   / ( (gamma+1)/(gamma-1) * (1 - (pe/pc)**(1-1/gamma)) )**.5 
    
    Aexit = Athroat * expansion_ratio
    
    rexit = np.sqrt(Aexit/np.pi)
    
    dexit = 2 * rexit
    
    return Aexit, rexit, dexit

###############################################################################
    
def misc_char(rthroat,Aexit,Athroat,Lstar):    
   
    rn = 0.382 * rthroat
    
    rb = 1.5 * rthroat
    
    expansion_ratio = Aexit/Athroat
    
    Vol_chamber = Lstar * Athroat

    return rn, rb, expansion_ratio, Vol_chamber

###############################################################################
    
def nozzle_length(dthroat,aexit,athroat,alphad,adj_coeff,Vc,Ac,rc,rt,betad):  
    
    alpha = np.deg2rad(alphad)
    
    beta = np.deg2rad(betad)
    
    l = dthroat/2 * (np.sqrt(aexit/athroat)-1) * np.arctan(alpha)
    
    l_nozzle = l * adj_coeff
    
    l_cylinder = Vc/Ac - (rc-rt)*np.arctan(beta)
    
    return l_nozzle, l_cylinder

###############################################################################
    
def engine_calculator(MW,gamma,Tci,pei,pai,pci,area_ratio,thrusti,MR,\
                      contraction_ratio, Lstari, alphad, betad, adj_coeff):
    
    thrust_m, Tc_m, pc_m, pa_m, pe_m, Lstar_m = imp2metric(thrusti,Tci,pci,pai,pei,Lstari)
    
    R_m, cstar_m = constants(MW,gamma,Tc_m)
    
    Ct, Isp = performance(gamma,pe_m,pc_m,pa_m,area_ratio,R_m,cstar_m)
    
    mdot_m, mdot_lox_m, mdot_fuel_m = massflow(thrust_m,Isp,MR)
    
    Tthroat_m, pthroat_m = throat_properties(Tc_m,pc_m,gamma)
    
#    uexhaust_m, Ma_exit = exit_properties(gamma,R_m,Tc_m,pc_m,pa_m,pe_m,tstag_m)
    
    Athroat_m, rthroat_m, dthroat_m = throat_geometry(cstar_m,mdot_m,pc_m)
    
    Achamber_m, rchamber_m, dchamber_m = chamber_geometry(Athroat_m,contraction_ratio)
    
    Aexit_m, rexit_m, dexit_m = exit_geometry(gamma,Athroat_m,pe_m,pc_m)
    
    rn_m, rb_m, expansion_ratio, Vol_chamber_m = misc_char(rthroat_m,Aexit_m,Athroat_m,Lstar_m)
    
    l_nozzle_m, l_cyl_m = nozzle_length(dthroat_m,Aexit_m,Athroat_m,alphad,adj_coeff,\
                                    Vol_chamber_m,Achamber_m,rchamber_m,rthroat_m,betad)
    
    cstar, mdot, mdot_oxidizer, mdot_fuel, Tthroat, pthroat, \
    Athroat, rthroat, Achamber, rchamber, Aexit, rexit, rn, rb, Vol_chamber, \
    l_nozzle, l_cyl = metric2imp(cstar_m,mdot_m,mdot_lox_m,mdot_fuel_m,\
                                 Tthroat_m,pthroat_m, Athroat_m, \
                                 rthroat_m, Achamber_m, rchamber_m, Aexit_m, \
                                 rexit_m, rn_m, rb_m, Vol_chamber_m, \
                                 l_nozzle_m, l_cyl_m)  
           
           
    return cstar, Ct, Isp, mdot, mdot_oxidizer, mdot_fuel, Tthroat, \
           pthroat, Athroat, rthroat, Achamber, rchamber, \
           Aexit, rexit, rn, rb, expansion_ratio, Vol_chamber, l_nozzle, l_cyl
           
###############################################################################
           
def engine_excel_writer(MW,gamma,Tci,pei,pai,pci,area_ratio,thrusti,MR,\
                      contraction_ratio, Lstar, alphad, betad, adj_coeff, \
                      filename):
    
    cstar, Ct, Isp, mdot, mdot_oxidizer, mdot_fuel, Tthroat, \
    pthroat, Athroat, rthroat, Achamber, rchamber, \
    Aexit, rexit, rn, rb, expansion_ratio, Vol_chamber, l_nozzle, l_cyl = \
    \
    engine_calculator(MW,gamma,Tci,pei,pai,pci,area_ratio,thrusti,MR,\
                      contraction_ratio, Lstar, alphad, betad, adj_coeff)
    
    workbook = xlsxwriter.Workbook(filename)
###############################################################################    
    Constants = workbook.add_worksheet('Characteristics')
    
    Constants.write('A1','Characteristic Velocity')
    Constants.write('A2','Coefficient of Thrust')
    Constants.write('A3','Specific Impulse')
    
    Constants.write('B1',cstar)
    Constants.write('B2',Ct)
    Constants.write('B3',Isp)
    
    Constants.write('C1','ft/s')
    Constants.write('C2','Unitless')
    Constants.write('C3','s')
###############################################################################    
    Massflow = workbook.add_worksheet('Massflow')
    
    Massflow.write('A1','Total Massflow')
    Massflow.write('A2','Massflow of Oxidizer')
    Massflow.write('A3','Massflow of Fuel')
    
    Massflow.write('B1',mdot)
    Massflow.write('B2',mdot_oxidizer)
    Massflow.write('B3',mdot_fuel)
    
    Massflow.write('C1','lbm/s')
    Massflow.write('C2','lbm/s')
    Massflow.write('C3','lbm/s')
###############################################################################    
    Properties = workbook.add_worksheet('Various Properties')
    
    Properties.write('A1','Throat Temperature')
    Properties.write('A2','Throat Pressure')
    
    Properties.write('B1',Tthroat)
    Properties.write('B2',pthroat)
    
    Properties.write('C1','Rankine')
    Properties.write('C2','psi')
###############################################################################   
    Geometry = workbook.add_worksheet('Engine Geometry')
    
    Geometry.write('A1','Throat Area')
    Geometry.write('A2','Throat Radius')
    Geometry.write('A3','Chamber Area')
    Geometry.write('A4','Chamber Radius')
    Geometry.write('A5','Radius before throat')
    Geometry.write('A6','Radius after throat')
    Geometry.write('A7','Exit Area')
    Geometry.write('A8','Exit Radius')
    Geometry.write('A9','Total Chamber Volume')
    Geometry.write('A10','Nozzle Length')
    Geometry.write('A11','Chamber Cylinder Length')
    
    Geometry.write('B1',Athroat)
    Geometry.write('B2',rthroat)
    Geometry.write('B3',Achamber)
    Geometry.write('B4',rchamber)
    Geometry.write('B5',rb)
    Geometry.write('B6',rn)
    Geometry.write('B7',Aexit)
    Geometry.write('B8',rexit)
    Geometry.write('B9',Vol_chamber)
    Geometry.write('B10',l_nozzle)
    Geometry.write('B11',l_cyl)

    Geometry.write('C1','in^2')
    Geometry.write('C2','in')
    Geometry.write('C3','in^2')
    Geometry.write('C4','in')
    Geometry.write('C5','in')
    Geometry.write('C6','in')
    Geometry.write('C7','in^2')
    Geometry.write('C8','in')
    Geometry.write('C9','in^3')
    Geometry.write('C10','in')
    Geometry.write('C10','in')
    Geometry.write('C11','in')
    
###############################################################################
    workbook.close()    