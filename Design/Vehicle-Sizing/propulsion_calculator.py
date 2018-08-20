import math
import Pressure Drop Calculator as pdc

#=============Global variables================
g = 32.17405 #Gravity constant in ft/s^2
atmosphericPressure = 14.696 #psi
efficiency = 0.45 #Efficiency of our engine (subject to possibly change?)
specificImpulse = 260 #Ratio of thrust produced to weight flow of propellants in seconds
exitArea = #exit area of the nozzle/engine

def openFile(filename):
    inputData = []
    with open(filename) as myFile:
        for line in myFile:
            inputData.append([int(x) for x in line.split()])
    return inputData


#chamberPressure = free stream pressure/chamber pressure
#mDot = mass flow rate
#exitV
def calculateThrust(mDot,equivalentVelocity):
    #The equation for thrust I am using is here https://www.grc.nasa.gov/WWW/K-12/airplane/specimp.html
    thrust = mDot * equivalentVelocity
    return thrust

def calculateMassFlowRate(chamberPressure,exitVelocity):
    #Solve for mDot from equivalent velocity here: https://www.grc.nasa.gov/WWW/K-12/airplane/specimp.html
    mDot = (exitArea * (atmosphericPressure - chamberPressure)) / (equivalentVelocity - exitVelocity)
    return mDot

def calculateEquivalentVelocity():
    #Equivalent velocity equations: https://www.grc.nasa.gov/WWW/K-12/airplane/specimp.html
    equivalentVelocity = specificImpulse * g * efficiency
    return equivalentVelocity

def calculateTotalPressures(pressureDrop,chamberPressure):
    propellantTankPressure = pressureDrop + chamberPressure
    return propellantTankPressure

def calculateExitVelocity():
    #Use the formula specified here: https://www.grc.nasa.gov/WWW/K-12/airplane/rktthsum.html
    #Or use another forumula?
    #WIP
    return exitVelocity

def main():
    grabFile = input("Hi there! What is the file name: ")
    chamberPressures = openFile(grabFile)
    Ve = calculateExitVelocity()
    Veq = calculateEquivalentVelocity()
    for value in chamberPressures:
        mDot = calculateMassFlowRate(value,Ve)
        thrust = calculateThrust(mDot,Veq)
        #Calculate pressure drop through pdc
        for
            pD = pdc.pressure_drop() #TODO: SUPPLY ARGUMENTS
            pressureDrop = pdc.pressure_drop_calc() #TODO: SUPPLY ARGUMENTS
            totalPressure = calculateTotalPressures(pressureDrop,value)
            #Use Neal's code to obtain pressure drops
            #WIP

        #TODO: WRITE VALUES TO OUTPUT FILE



main()
