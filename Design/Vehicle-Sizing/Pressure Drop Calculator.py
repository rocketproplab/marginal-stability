import math
from scipy.optimize import fsolve 

def barlows(pressure,outer_diameter,stress,fos):
    t = (pressure*outer_diameter*fos)/(2*stress)

#Thickness values picked from easily purchased 304 stainless steel tubing from onlinemetals.com
    if t < 0.020:
        thickness = 0.020
    elif t > 0.020 and t < 0.028:
        thickness = 0.028
    elif t > 0.028 and t < 0.035: #max thickness flaring tool can handle is 0.035 inches, so it is treated as a hard cap
        thickness = 0.035
    else:
        raise ValueError
    
    inner_diameter = outer_diameter - 2*thickness
    area = inner_diameter**2 * math.pi/4
    
    return thickness, inner_diameter, area

class lox():
    def __init__(self):
        self.density = 0.0412212402 #lb/in^3 - https://en.wikipedia.org/wiki/Liquid_oxygen
        self.kv = 0.0024336 #kinematic viscosity, in^2/s - https://www.engineeringtoolbox.com/oxygen-O2-dynamic-kinematic-viscosity-temperature-pressure-d_2081.html?vA=-183&degree=C#
        self.sg = 1.14 #specific gravity @ 20C, 1 atm - http://www.airproducts.com/~/media/files/pdf/company/safetygram-6.pdf

class rp1():
    def __init__(self):
        self.density = 0.0292631065 #lb/in^3 - https://en.wikipedia.org/wiki/RP-1
        self.kv = 0.00333250667 #kinematic viscosity, in^2/s - https://kinetics.nist.gov/RealFuels/macccr/macccr2008/Bruno2.pdf
        self.sg = 0.82 #specific gravity assuming density of RP-1 is 820 kg/m^3 - https://kinetics.nist.gov/RealFuels/macccr/macccr2008/Bruno2.pdf
        
class helium():
    def __init__(self):
        self.density = ((pressure*0.068046)*4/(.0821*298))*3.61273e-5 #lb/in^3, ideal gas law SI -> imperial
        self.kv = (0.00001317059/12)/self.density #kinematic viscosity, in^2/s - https://www.engineeringtoolbox.com/gases-absolute-dynamic-viscosity-d_1888.html 20C
        self.sg = 0.138 #https://www.engineeringtoolbox.com/specific-gravities-gases-d_334.html
        
def fluid_properties(inputfluid):
    fluid = inputfluid
    density = fluid.density
    kv = fluid.kv
    sg = fluid.sg
    
    return density,kv,sg

def pressuredrop(massflow,density,area,kv,pipe_roughness,k_factor,pressure,pipe_length):
    
    gravity = 32.17405 #ft/s^2 - https://www.engineeringtoolbox.com/mass-weight-d_589.html

#Based on massflow input, calculate the velocity. Assumed large tank, so only velocity considered is that at point 2
#Assume massflow used by engine is a constant value along all points outside the tank
    velocity = massflow / (density*area)

    reynolds = velocity*inner_diameter/(kv) #https://www.engineeringtoolbox.com/reynolds-number-d_237.html - circular pipe

#friction factor calculation for laminar flow
#https://www.nuclear-power.net/nuclear-engineering/fluid-dynamics/major-head-loss-friction-loss/friction-factor-for-laminar-flow/
    if reynolds < 2000:
        f = 64/reynolds
    
#friction factor calculation for turbulent flow with Re < 100,000
#https://www.nuclear-power.net/nuclear-engineering/fluid-dynamics/major-head-loss-friction-loss/friction-factor-turbulent-flow-colebrook/
    if reynolds > 2000:
        def friction(f):
            return (-2*math.log10(2.51/(reynolds*math.sqrt(f))+((pipe_roughness/inner_diameter)/3.72)))**2 - 1/f

        f = fsolve(friction,0.01)

#####################################################################################################################

#Can't find another source that uses this equation for such a relatively low Reynold's number, uncomment if verified

##friction factor calculation for turbulent flow with Re > 100,000    
##https://www.pipeflowcalculations.com/pipe-valve-fitting-flow/flow-in-pipes.php#Bernoulli-theorem
#    if reynolds > 100000:
#        f = (2.0*math.log10(inner_diameter/pipe_roughness)+1.14)**(-2)

#####################################################################################################################    
    
    density_calc = density*0.031081*(12**3) #lb/in^3 -> slug/ft^3
    pipe_length_calc = pipe_length/12 #in -> ft
    velocity_calc = velocity/12 #in/s -> ft/s
    inner_diameter_calc = inner_diameter/12 #in -> ft
    
    head_loss_major = (f*pipe_length_calc*velocity_calc**2)/(2*inner_diameter_calc*gravity)
    head_loss_minor = (k_factor*velocity_calc**2)/(2*gravity)
    
    total_head_loss = head_loss_major + head_loss_minor
    
    return reynolds, f, total_head_loss, velocity

def pressure_drop_calc(length,pressure_initial,density,velocity,head_loss): 
    gravity = 32.17405
    density_calc = density*0.031081*(12**3) #lb/in^3 -> slug/ft^3
    length_calc = length/12 #in -> ft
    velocity_calc = velocity/12 #in/s -> ft/s
    
    pressure_change = head_loss * density_calc * gravity - density_calc * length_calc * gravity + 1/2 * velocity_calc**2 * density_calc
    psi_change = pressure_change/144
    
    return psi_change

def cv_calc(massflow,sg,delP,density):
    Q = massflow/density #in^3/s
    vol_flow = 0.004329 * Q *60 #gal/min
    Cv = vol_flow * math.sqrt(sg)/math.sqrt(delP)
    
    return Cv

pressure = 1000 #psi
outer_diameter = 0.5 #inches
stress = 31200 #yield strength of material, psi
fos = 1.5 #factor of safety
massflow = 20 #lb/s
pipe_roughness = 0.00059055118 #converted mm to inches - http://www.enggcyclopedia.com/2011/09/absolute-roughness/
k_factor = 6 #addition of total k-factors for valves, fitings, etc - https://www.lmnoeng.com/surface.htm
pipe_length = 60 #total length of pipe - inches

thickness, inner_diameter, area = barlows(pressure,outer_diameter,stress,fos)
density, kv,sg = fluid_properties(lox())
reynolds, f, head_loss,velocity = pressuredrop(massflow,density,area,kv,pipe_roughness,k_factor,pressure,pipe_length)
psi_change = pressure_drop_calc(pipe_length,pressure,density,velocity,head_loss)
cv = cv_calc(massflow,sg,psi_change,density)

print(cv)