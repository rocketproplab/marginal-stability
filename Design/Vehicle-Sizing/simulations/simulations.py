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
                coef_drag           # [dimensionless]

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
            'coef_drag'
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
        self.Cd                 = self.initialConditions['coef_drag']
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
        self.R               = [self.Rearth + self.initialConditions['altitude']]  # [ft] initial distance to the center of the earth
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

            # calculate altitude, velocity, and acceleration
            self.altitude.append(self.altitude[self.runIter] + self.calc_dalt())
            self.velocity.append(self.velocity[self.runIter] + self.calc_deltaV())
            self.drag.append(self.calc_drag(self.velocity[self.runIter], rho, self.reference_area, self.Cd))
            #print(str(self.Cd) + str(self.drag[self.runIter]))
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
            if self.velocity[self.runIter] < 0:
                self.apogee = True

            if not self.apogee:
                self.dynamic_pressure.append(self.maxQ(rho))

            # END CONDITIONS
            if (self.altitude[self.runIter] < 1000 and self.time[self.runIter] > self.burntime) or self.time[self.runIter] > 10000:
                break

            self.runIter += 1
        return (self.altitude, self.velocity, self.acceleration, self.mass, self.time, self.thrust, self.drag, self.dynamic_pressure, self.rho, self.temp, self.M)

    def calc_drag(self, vel, rho, S, Cd):
        return (1/2)*rho*(vel**2)*S*Cd

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
            thrust:         # [lbf]
            thrust_angle:   # [rad] angle
            drag:           # [lbf]
            mass:           # [slug]
            flight_heading:  # [rad] angle
            R:              # [ft] radius from center of the Earth
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
            thrust:         # [lbf]
            thrust_angle:   # [rad] angle
            drag:           # [lbf]
            mass:           # [slug]
            flight_heading: # [rad] angle
            R:              # [ft] radius from center of the Earth
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
            velocity:       # [ft/s] velocity at current timestep
            flight_heading: # [ft] radius from center of the Earth
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
            return (360,0,940)

def test_Rocket():
    burntime = 50  # s
    nengines = 1
    thrust_sl = 24000
    Isp = 400
    g0 = 9.81
    mdot = nengines*thrust_sl/(g0*Isp)
    twratio = 50  # estimated thrust 2 weight ratio
    mstructure = 300  # slug
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
