import math
import Pressure_Drop_Calculator as pdc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

#=============Global variables================
g = 32.17405 #Gravity constant in ft/s^2
atmosphericPressure = 14.696 #psi
efficiency = 0.45 #Efficiency of our engine (subject to possibly change?)
specificImpulse = 260 #Ratio of thrust produced to weight flow of propellants in seconds
maxTankPressure = 750 #psi
#=============================================

def openFile(filename):
    inputData = []
    with open(filename) as myFile:
        for line in myFile:
            for i, value in enumerate(line.split()):
                inputData.append(int(value))
    return inputData

#def calculateThrust(mDot,equivalentVelocity):
    #The equation for thrust I am using is here https://www.grc.nasa.gov/WWW/K-12/airplane/specimp.html
    #thrust = mDot * equivalentVelocity
    #return thrust

def calculateMassFlowRate(thrust,equivalentVelocity):
    #Solve for mDot from equivalent velocity here: https://www.grc.nasa.gov/WWW/K-12/airplane/specimp.html
    #mDot = (exitArea * (atmosphericPressure - chamberPressure)) / (equivalentVelocity - exitVelocity)
    mDot = (thrust * g) / equivalentVelocity
    return mDot

def calculateEquivalentVelocity():
    #Equivalent velocity equations: https://www.grc.nasa.gov/WWW/K-12/airplane/specimp.html
    equivalentVelocity = specificImpulse * g * efficiency
    return equivalentVelocity

def calculateTotalPressures(pressureDrop):
    chamberPressure = maxTankPressure - pressureDrop
    return chamberPressure

def main():
    lox = pdc.lox()
    rp1 = pdc.rp1()
    helium = pdc.helium()

    #=========MANIPULATE THESE VALUES FOR PRESSURE DROP FUNCTION===========
    pressure = 500 #Max tank pressure in psi
    stress = 31200 #yield strength of material, psi
    fos = 1.5 #factor of safety
    pipe_roughness = 0.00059055118 #converted mm to inches - http://www.enggcyclopedia.com/2011/09/absolute-roughness/
    k_factor = 6 #addition of total k-factors for valves, fitings, etc - https://www.lmnoeng.com/surface.htm
    pipe_length = 60 #total length of pipe - inches
    fluidname = lox #change between lox,rp1,helium
    #==============================================================

    grabThrustFile = input("Hi there! What is the thrust file name: ")
    thrustPossibilities = openFile(grabThrustFile)

    Veq = calculateEquivalentVelocity()
    pipeSizes = [0.5,0.625,0.75,0.875,1]

    x1 = []
    y1 = []
    x2 = []
    y2 = []

    inches = " inches"

    plt.figure(1)
    for pipeSize in pipeSizes:
        for value in thrustPossibilities:
            mDot = calculateMassFlowRate(value,Veq)
            pD = pdc.exportfunc(pressure,pipeSize,stress,fos,fluidname,mDot,k_factor,pipe_roughness,pipe_length)
            chamberPressure = calculateTotalPressures(pD)
            x1.append(value)
            y1.append(chamberPressure)
            y2.append(pD)
        plt.subplot(2,1,1)
        plt.plot(x1,y1,label=str(pipeSize)+inches)
        plt.subplot(2,1,2)
        plt.plot(x1,y2,label=str(pipeSize)+inches)
        del x1[:]
        del y1[:]
        del y2[:]

    plt.subplot(2,1,1)
    plt.plot(x1,y1)
    plt.axis([0,5000,0,1000])
    plt.xlabel('Thrust')
    plt.ylabel('Chamber Pressure')
    plt.legend(loc='lower left')

    plt.subplot(2,1,2)
    plt.plot(x1,y2)
    plt.axis([0,5000,0,1000])
    plt.xlabel('Thrust')
    plt.ylabel('Pressure Drop')
    plt.legend(loc='upper left')

    plt.show()

main()
