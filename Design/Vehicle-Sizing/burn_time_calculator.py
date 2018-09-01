import numpy
import MS_TANK_SIZING_CODE as sizing

altitudeGoal = int(input("Enter altitude goal: "))
thrust = int(input("Enter thrust: "))
diameter = int(input("Enter rocket diameter: "))
g0 = 9.8 #m/s^2 (acceleration due to gravity)
specificImpulse = 260 #s
of_ratio = 2.33
drag = 300
#http://web.mit.edu/16.unified/www/SPRING/propulsion/notes/node103.html has rocket equations

def lbToKg(lbMass):
    return lbMass / 2.20462

def kgToLb(kgMass):
    return kgMass * 2.20462

def inToM(inches):
    return inches * 2.54 / 100

def fixedMass():
    fixedMass = 470.25 - sizing.wall_mass
    return lbToKg(fixedMass)

def dryMass(fixedMass):
    return fixedMass + lbToKg(sizing.wall_mass)

def wetMass(dryMass):
    return dryMass + lbToKg(sizing.m_lox) + lbToKg(sizing.m_rp1)

def massRatio(wetMass, dryMass):
    ratio = wetMass / dryMass
    return ratio

def burnTime(specificImpulse, massRatio, g0, altitudeGoal):
    rFrac = (massRatio*numpy.log(massRatio))/(massRatio - 1)
    v_e = specificImpulse * g0
    num1 = v_e * (numpy.log(massRatio)**2)
    denom1 = 2*g0*(rFrac - 1)
    num2 = altitudeGoal
    denom2 = v_e * (rFrac - 1)
    burnTime = (num1/denom1) - (num2/denom2)
    return burnTime

def rp1MassCalc(propMass):
    return propMass / (of_ratio + 1)

def loxMassCalc(propMass):
    return propMass * of_ratio / (of_ratio + 1)

def heightCalc(specificImpulse, massRatio, burnTime, g0):
    v_e = specificImpulse * g0
    hMax = ((v_e)**2 * (numpy.log(massRatio))**2 / (2 * g0)) - v_e * burnTime * ((massRatio/(massRatio - 1))*numpy.log(R) - 1)
    return hMax

sizing.changeDiam(diameter)
F = fixedMass()
mdot = thrust / (specificImpulse * g0)
D = dryMass(F)
W = wetMass(D)
R = massRatio(W, D)
dragAccel = 2*drag/(D + W)
print(dragAccel)
t = burnTime(specificImpulse, R, g0, altitudeGoal)
propMass = t * mdot

print("D1 = " + str(kgToLb(D)))
print("W1 = " + str(kgToLb(W)))
print("Old propMass = " + str(sizing.m_lox + sizing.m_rp1))
print("t1 = " + str(t))
print("height = " + str(heightCalc(specificImpulse, R, t, g0)))
print("mdot = " + str(kgToLb(mdot)))
print("new propMass = " + str(kgToLb(propMass)))

for count in range(20):
    rp1Mass = rp1MassCalc(propMass)
    loxMass = loxMassCalc(propMass)
    sizing.changeRP1(kgToLb(rp1Mass))
    sizing.changeLOx(kgToLb(loxMass))
    newD = dryMass(F)
    newW = wetMass(newD)
    newR = massRatio(newW, newD)
    dragAccel = 2*drag/(newD + newW)
    newT = burnTime(specificImpulse, newR, g0, altitudeGoal)
    propMass = newT * mdot
    print("new propMass = " + str(kgToLb(propMass)))
    print("new D = " + str(kgToLb(newD)))
    print("new W = " + str(kgToLb(newW)))
    print("new T = " + str(newT))
    print("height = " + str(heightCalc(specificImpulse, newR, newT, g0)))
    print("mdot = " + str(kgToLb(mdot)))
