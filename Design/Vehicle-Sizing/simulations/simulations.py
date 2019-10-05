""" simulations.py
Simulations is a python based flight simulation package
for rocket and missle trajectory analysis. """
import numpy as np
import csv
from scipy import interpolate as inter

__author__ = "Cameron Flannery"
__copyright__ = "Copyright 2018"
__license__ = "MIT"
__version__ = "0.0.1"
__status__ = "alpha"

class Rocket(object):
    """ Rocket is a simulation class for rocket simulations.

    This module performs calculations for the estimation of
    launch vehicle sizing and trajectory simulation for
    vertical launch vehicles.

    A list of assumptions, capabilities, and limitations
    will be added here as features are solidified. """

    def import_thrust_vec(self):
        unit_convert = unit()
        tlist = []
        alist = []
        with open('simulations/ThrustAlt.csv') as csvfile:
            tareader = csv.reader(csvfile)
            names = False
            for i in tareader:
                if not names:
                    names = True
                else:
                    tlist.append(float(i[0]))
                    alist.append(unit_convert.mToFt(float(i[1])*1000))
        return tlist, alist

    def get_thrust(self,alt,t_list,a_list):
        spline = inter.CubicSpline(a_list,t_list)
        return spline(alt)

    def __init__(self, initialConditions, engines, burntime, timestep=0.5):
        """ Initialization of the Rocket simulation class

        Args:
            initialConditions: expects dictionary containing
            initial conditions for launch vehicle.
             -> Required keywords:
                time,               # [s]
                velocity,           # [ft/s]
                flight_angle,       # [rad] vertical flight path angle
                flight_heading,     # [rad] flight path heading
                latitude,           # [rad]
                longitude,          # [rad]
                altitude,           # [ft]
                mass,               # [slug]
                thrust_sl,          # [lbf] sea-level thrust
                thrust_angle,       # [rad]
                lift_coefficient,   # [1]
                bank_angle          # [rad]
                reference_area      # [ft^2]

            engines:
            -> Required keywords:
                thrust_sl:          # [lbf]
                Isp:                # [s]
                Ae:                 # [ft^2]
                nengines:           # [1]

            burntime: length of burn in seconds
        Keyword Args:
            timestep: (optional), timestep in seconds. Default
                timestep is 1s
        Returns:
            0: Completed with no errors
        """

        requiredArgs = [
            'time',
            'velocity',
            'flight_angle',
            'flight_heading',
            'latitude',
            'longitude',
            'altitude',
            'mass',
            'lift_coefficient',
            'bank_angle'
            'thrust_sl',
            'thrust_angle',
            'Ae',
            'Isp',
            'reference_area',
        ]

        for arg in requiredArgs:
            if arg not in (initialConditions or engines):
                pass

        self.initialConditions = initialConditions
        self.engines = engines
        self.burntime = burntime
        self.timestep = timestep
        self.CONST()
        self.atm = atm()
        self.unit = unit()

    def run(self, stopTime=None, stopApogee=None):
        """ runs simulation

        Automatically ends simulation when the vehicle impacts
        the ground, or reaches a stable orbit.
        """

        # initialize arrays with values from initialConditions
        self.time               = [self.initialConditions['time']]
        self.velocity           = [self.initialConditions['velocity']]
        self.flight_angle       = [self.initialConditions['flight_angle']]
        self.flight_heading     = [self.initialConditions['flight_heading']]
        self.latitude           = [self.initialConditions['latitude']]
        self.longitude          = [self.initialConditions['longitude']]
        self.altitude           = [self.initialConditions['altitude']]
        self.mass               = [self.initialConditions['mass']]
        self.lift_coefficient   = [self.initialConditions['lift_coefficient']]
        self.bank_angle         = [self.initialConditions['bank_angle']]
        self.Cd                 = [self.calc_Cd(0)]
        self.drag               = [0]
        self.dynamic_pressure   = [0]
        self.rho                = [self.STDATM(self.initialConditions['altitude'])[1]]
        self.temp               = [self.STDATM(self.initialConditions['altitude'])[0]]
        self.M                  = [0]

        # initialize arrays with values from engines
        self.nengines           = self.engines['nengines']
        self.t_vec, self.a_vec  = self.import_thrust_vec()
        self.thrust             = [self.t_vec[0]]
        self.thrust_angle       = [self.engines['thrust_angle']]
        self.Ae                 = [self.engines['Ae']]
        self.Isp                = self.engines['Isp']
        self.mdot               = self.thrust[0]/(self.g0*self.Isp)

        # initialize additional values
        self.acceleration    = [0]
        self.R               = [self.Rearth]  # [m] initial distance to the center of the earth
        self.reference_area  = self.initialConditions['reference_area']

        self.runIter = 0  # iterator
        self.apogee = False

        while True:
            self.time.append(self.time[self.runIter] + self.timestep)
            self.R.append(self.Rearth + self.altitude[self.runIter])
            T, rho, sos = self.STDATM(self.altitude[self.runIter])  # Thermoproperties
            self.rho.append(rho)
            self.temp.append(T)

            self.M.append(self.velocity[self.runIter]/sos)
            Cd = self.calc_Cd(self.M[self.runIter])

            # calculate altitude, velocity, and acceleration
            self.altitude.append(self.altitude[self.runIter] + self.calc_dalt())
            self.velocity.append(self.velocity[self.runIter] + self.calc_deltaV())
            self.drag.append(self.calc_drag(self.velocity[self.runIter], rho, self.reference_area, Cd))
            self.acceleration.append(self.calc_accel())

            # Thrust

            if self.time[self.runIter] <= self.burntime:
                self.thrust.append(self.get_thrust(self.altitude[self.runIter + 1], self.t_vec, self.a_vec))
                print(self.thrust[self.runIter])
                self.mass.append(self.mass[self.runIter] - self.mdot*self.timestep)
            else:
                self.thrust.append(0)
                self.mass.append(self.mass[self.runIter])
            self.flight_heading.append(self.flight_heading[0])  # initial values until calcs added
            self.flight_angle.append(self.flight_angle[0])      # initial values until calcs added
            self.thrust_angle.append(self.thrust_angle[0])
            self.Cd.append(self.Cd[0])
            if self.velocity[self.runIter] < 0:
                self.apogee = True

            if not self.apogee:
                self.dynamic_pressure.append(self.maxQ(rho))

            # END CONDITIONS
            if (self.altitude[self.runIter] < 1000 and self.time[self.runIter] > self.burntime) or self.time[self.runIter] > 10000:
                break

            self.runIter += 1
        return (self.altitude, self.velocity, self.acceleration, self.mass, self.time, self.thrust, self.drag, self.dynamic_pressure, self.rho, self.temp, self.M)

    def calc_Cd(self, M):
        return .52

    def calc_drag(self, vel, rho, S, Cd):
        return 1/2*rho*vel**2*S*Cd

    def calc_thrust(self, thrust_sl=None, Ae=None, pe=None, pa=None):
        """ calc_thrust determines the thrust """

        return thrust_sl + (pe-pa)*Ae

    def calc_accel(self, thrust=None, thrust_angle=None,
                   drag=None, mass=None, g0=None, R=None,
                   flight_heading=None, timestep=None):
        """ calc_accel is a method of Rocket

        This method is typically used to update the acceleration values during
        simulation runs.

        Note:
            All arguments are optional. If no arguments are thrown,
            the method will return the calculated delta velocity
            based on the current timestep values.

        Args:
            thrust:         # [N]
            thrust_angle:   # [rad] angle
            drag:           # [N]
            mass:           # [kg]
            flight_heading:  # [rad] angle
            R:              # [m] radius from center of the Earth
            timestep:        # [s] timestep

        Returns:
            accel: acceleration value of the rocket
        """

        i = self.runIter
        if not thrust:
            thrust = self.thrust[i]
        if not thrust_angle:
            thrust_angle = self.thrust_angle[i]
        if not drag:
            drag = self.drag[i]
        if not mass:
            mass = self.mass[i]
        if not R:
            R = self.R[i]
        if not flight_heading:
            flight_heading = self.flight_heading[i]
        if not timestep:
            timestep = self.timestep

        return self.dVdt(thrust, thrust_angle, drag, mass, R, flight_heading)

    def calc_deltaV(self, thrust=None, thrust_angle=None,
                    drag=None, mass=None, g0=None, R=None,
                    flight_heading=None, timestep=None):
        """ calc_velocity is a method of Rocket

        This method is typically used to update the velocity values during
        simulation runs.

        Note:
            All arguments are optional. If no arguments are thrown,
            the method will return the calculated delta velocity
            based on the current timestep values.

        Args:
            thrust:         # [N]
            thrust_angle:   # [rad] angle
            drag:           # [N]
            mass:           # [kg]
            flight_heading: # [rad] angle
            R:              # [m] radius from center of the Earth
            timestep:       # [s] timestep

        Returns:
            dV: change in velocity of the rocket
        """

        i = self.runIter
        if not thrust:
            thrust = self.thrust[i]
        if not thrust_angle:
            thrust_angle = self.thrust_angle[i]
        if not drag:
            drag = self.drag[i]
        if not mass:
            mass = self.mass[i]
        if not R:
            R = self.R[i]
        if not flight_heading:
            flight_heading = self.flight_heading[i]
        if not timestep:
            timestep = self.timestep

        dV = self.dVdt(thrust, thrust_angle, drag, mass, R, flight_heading)*timestep
        return dV

    def calc_dalt(self, velocity=None, flight_heading=None, timestep=None):
        """ calc_dalt is a method of Rocket

        This method is typically used to update the altitude values during
        simulation runs.

        Note:
            All arguments are optional. If no arguments are thrown,
            the method will return the calculated delta altitude
            based on the current timestep values.

        Args:
            velocity:       # [m/s] velocity at current timestep
            flight_heading: # [m] radius from center of the Earth
            timestep:       # [s] timestep

        Returns:
            dalt: change in velocity of the rocket
        """

        i = self.runIter
        if not velocity:
            velocity = self.velocity[i]
        if not flight_heading:
            flight_heading = self.flight_heading[i]
        if not timestep:
            timestep = self.timestep

        dalt = velocity*np.sin(flight_heading)*timestep
        return dalt

    def dVdt(self, thrust, thrust_angle, drag, mass, R, flight_heading):
        """ dVdt is a method of Rocket

        dVdt calculates the acceleration of the rocket at the current step
        """
        dVdt = (thrust*np.cos(thrust_angle)-drag)/mass - self.g0*(self.Rearth/R)**2*np.sin(flight_heading)
        return dVdt

    def maxQ (self,density):
        i = self.runIter
        currentVelocity = self.velocity[i]
        return (1/2) * density * (currentVelocity ** 2)

    def CONST(self):
        """ Define useful constants as instance variables """
        self.g0 = 32.17          # gravity constant [ft/s^2]
        self.R_air = 1716.49        # gas constant [ft*lbf/(slug*R)]
        self.gamma_air = 1.4    # ratio of specific heats
        self.Rearth = 2.0902 * 10**7   # [ft]

    #Standard atmosphere interpolation:
    def STDATM(self,alt):
        if alt < 250000:
            return (self.atm.get_temp(alt),self.atm.get_rho(alt),self.atm.get_sos(alt))
        else:
            return (0,0,940)

    # standard atmosphere model (SI units)
    '''def STDATM(self, altitude):
        layer = -1.0            # gradient layer
        gradient = -self.unit.kToR(0.0019812)
        altitude_base = 0.0
        temperature_base = 518.67
        density_base = 0.00237717

        if altitude > 36089:
            layer = 1.0       # isothermal layer
            altitude_base = 36089
            temperature_base = self.unit.kToR(216.65)
            density_base = 0.000706208
        if altitude > 65617:
            layer = -1.0      # gradient layer
            gradient = self.unit.kToR(0.0003048)
            altitude_base = 65617
            temperature_base = 216.65 * 9/5
            density_base = 0.000170834
        if altitude > 104987:
            layer = -1.0       # gradient layer
            gradient = self.unit.kToR(0.00085344)
            altitude_base = 104987
            temperature_base = 228.65*9/5
            density_base = 0.0000256636
        if altitude > 154199:
            layer = 1.0      # isothermal layer
            altitude_base = 154199
            temperature_base = 270.65*9/5
            density_base = 1.47056878 * (10**(-6))
        if altitude > 167323:
            layer = -1.0       # gradient layer
            gradinet = -self.unit.kToR(0.00085344)
            altitude_base = 167323
            temperature_base = 270.65*9/5
            density_base = 0.00000277025
        if altitude > 232940:
            layer = 1.0      # isothermal layer
            altitude_base = 232940
            temperature_base = 214.65*9/5
            density_base = 1.24603 * (10**(-7))
        if layer < 0.0:
            temperature = temperature_base + gradient*(altitude - altitude_base)
            power = -1.0*(self.g0/gradient/self.R_air + 1.0)
            density = density_base*(temperature/temperature_base)**power
        else:
            temperature = temperature_base
            power = -1.0*self.g0*(altitude - altitude_base)/self.R_air/temperature
            density = density_base*np.exp(power)
        sos = np.sqrt(self.gamma_air*self.R_air*temperature)

        meAlt = self.unit.ftToM(altitude)
        geAlt = self.atm.geopotentialAlt(meAlt)
        meDensity = self.atm.getDensity(geAlt)
        meTemp = self.atm.getTempK(geAlt)
        temperature = self.unit.kToR(meTemp)
        density = self.unit.kgm3ToSlugFt3(meDensity)
        meSos = np.sqrt(self.gamma_air*self.R_air*temperature)
        sos = self.unit.mToFt(meSos)

        return (temperature, density, sos)'''


def test_Rocket():
    burntime = 50  # s
    nengines = 1
    thrust_sl = 24000
    Isp = 400
    g0 = 9.81
    mdot = nengines*thrust_sl/(g0*Isp)
    twratio = 50  # estimated thrust 2 weight ratio
    mstructure = 300  # kg
    mpropulsion = thrust_sl/(twratio*g0)
    mpropellant = mdot*burntime
    mass = mpropulsion + mpropellant + mstructure
    initialConditions = {
        'time': 0,
        'velocity': 0,
        'flight_angle': 0,
        'flight_heading': np.deg2rad(90),
        'latitude': 0,
        'longitude': 0,
        'altitude': 0,
        'mass': mass,
        'lift_coefficient': 0,
        'bank_angle': 0,
        'reference_area' : 0
    }
    engines = {
        'thrust_sl': thrust_sl,
        'thrust_angle': 0,
        'Isp': Isp,
        'Ae': 0.25,
        'nengines': nengines
    }
    itsatest = Rocket(initialConditions,engines,burntime)
    altitude, velocity, acceleration, mass, time, thrust = itsatest.run()

    return 0

class atm(object):
    def __init__(self):
        self.altlist = []
        self.tlist = []
        self.rholist = []
        self.soslist = []
        with open('simulations/atm.csv') as atmcsv:
            atmreader = csv.reader(atmcsv)
            names = False
            for i in atmreader:
                if not names:
                    names = True
                else:
                    self.altlist.append(float(i[0]))
                    self.tlist.append(float(i[1]))
                    self.rholist.append(float(i[2]))
                    self.soslist.append(float(i[3]))
        self.tspline = inter.CubicSpline(self.altlist,self.tlist)
        self.rhospline = inter.CubicSpline(self.altlist,self.rholist)
        self.sosspline = inter.CubicSpline(self.altlist,self.soslist)

    def get_temp(self,alt):
        return self.tspline(alt)

    def get_rho(self,alt):
        return self.rhospline(alt)

    def get_sos(self,alt):
        return self.tspline(alt)

'''class atm(object):
    def __init__(self):
        self.molWeight = 28.9644
        self.g = 9.8
        self.R = 8.314

    def getTempK(self, geAlt):
        temp = 288.15
        if 0 < geAlt < 11000:
            temp = 288.15 - 0.0065 * geAlt
        elif geAlt < 20000:
            temp = 216.65
        elif geAlt < 32000:
            temp = 216.65 + 0.001 * geAlt
        elif geAlt < 47000:
            temp = 228.65 + 0.0028 * geAlt
        elif geAlt < 51000:
            temp = 270.65
        elif geAlt < 71000:
            temp = 270.65 - 0.0028 * geAlt
        elif geAlt >= 71000:
            temp = 214.65 - 0.002 * geAlt
        return temp

    def getTempGrad(self, geAlt):
        if geAlt < 11000:
            return 0.0065
        elif geAlt < 20000:
            return 0
        elif geAlt < 32000:
            return -0.001
        elif geAlt < 47000:
            return -0.0028
        elif geAlt < 51000:
            return 0
        elif geAlt < 71000:
            return 0.0028
        elif geAlt >= 71000:
            return 0.002

    #Using this corrected altitude instead of the actual altitude allows us to
    #treat gravity as constant. It is a fairly small difference at low
    #altitudes, since gravity doesn't change very quickly
    def geopotentialAlt(self,alt):
        return alt

    def getBasePress(self, geAlt):
        if geAlt < 11000:
            return 101325
        elif geAlt < 20000:
            return 22632.1
        elif geAlt < 32000:
            return 5474.89
        elif geAlt < 47000:
            return 868.019
        elif geAlt < 51000:
            return 110.906
        elif geAlt < 71000:
            return 66.9389
        else:
            return 3.95642

    def getBaseAlt(self, geAlt):
        if geAlt < 11000:
            return 0
        elif geAlt < 20000:
            return 11000
        elif geAlt < 32000:
            return 20000
        elif geAlt < 47000:
            return 32000
        elif geAlt < 51000:
            return 47000
        elif geAlt < 71000:
            return 51000
        else:
            return 71000

    def getPressure(self, geAlt):
        temp = self.getTempK(geAlt)
        baseAlt = self.getBaseAlt(geAlt)
        baseTemp = self.getTempK(baseAlt)
        basePress = self.getBasePress(geAlt)
        tempGrad = self.getTempGrad(geAlt)
        if tempGrad == 0:
            num = -self.g * self.molWeight * (geAlt - baseAlt)
            den = self.R * baseTemp
            return basePress * np.exp(num/den)
        else:
            num = -self.g * self.molWeight
            den = self.R * tempGrad
            return basePress * (temp / baseTemp)**(num/den)

    def getDensity(self, geAlt):
        temp = self.getTempK(geAlt)
        pressure = self.getPressure(geAlt)
        return pressure * self.molWeight /(self.R * temp)'''

class unit(object):
    def kgToLb(self, kg):
        return kg * 2.20462

    def lbToKg(self, lb):
        return lb / 2.20462

    def lbfToN(self, lbf):
        return lbf * 4.44822

    def ftToM(self, ft):
        return ft / 3.28084

    def mToFt(self, m):
        return m * 3.28084

    def kmToMi(self, km):
        return km / 1.60934

    def paToPsi(self, pa):
        return pa / 6894.76

    def kgm3ToLbft3(self, kgm3):
        return kgm3 / 16.0185

    def nToLbf(self, n):
        return n / 4.44822

    def lbmToSlug(self, lbm):
        return lbm / 32.174049

    def psfToPsi(self, psf):
        return psf / 144

    def kToR(self, k):
        return k * 9/5

    def kgm3ToSlugFt3(self, kgm3):
        return kgm3 * 0.00194032

    def inToMm(self, inch):
        return inch*25.4

    def mmToIn(self, mm):
        return mm/25.4

if __name__ == '__main__':
    test_Rocket()
