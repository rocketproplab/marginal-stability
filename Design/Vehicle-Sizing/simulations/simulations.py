""" simulations.py
Simulations is a python based flight simulation package
for rocket and missle trajectory analysis. """
import numpy as np

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

    def __init__(self, initialConditions, engines, burntime, timestep=1):
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

        # initialize arrays with values from engines
        self.nengines           = self.engines['nengines']
        self.thrust             = [self.engines['thrust_sl']*self.nengines]
        self.thrust_angle       = [self.engines['thrust_angle']]
        self.Ae                 = [self.engines['Ae']]
        self.Isp                = self.engines['Isp']
        self.mdot               = self.thrust[0]/(self.g0*self.Isp)

        # initialize additional values
        self.acceleration    = [0]
        self.R               = [self.Rearth]  # [m] initial distance to the center of the earth
        self.reference_area  = self.initialConditions['reference_area']

        self.runIter = 0  # iterator
        while True:
            self.time.append(self.time[self.runIter] + self.timestep)
            self.R.append(self.Rearth + self.altitude[self.runIter])
            T, rho, sos = self.STDATM(self.altitude[self.runIter])  # Thermoproperties
            self.rho.append(rho)
            self.temp.append(T)

            M = self.velocity[self.runIter]/sos
            Cd = self.calc_Cd(M)

            # calculate altitude, velocity, and acceleration
            self.altitude.append(self.altitude[self.runIter] + self.calc_dalt())
            self.velocity.append(self.velocity[self.runIter] + self.calc_deltaV())
            self.drag.append(self.calc_drag(self.velocity[self.runIter], rho, self.reference_area, Cd))
            self.acceleration.append(self.calc_accel())

            # Thrust
            if self.time[self.runIter] <= self.burntime:
                self.thrust.append(self.thrust[0])
                self.mass.append(self.mass[self.runIter] - self.mdot*self.timestep)
            else:
                self.thrust.append(0)
                self.mass.append(self.mass[self.runIter])
            self.flight_heading.append(self.flight_heading[0])  # initial values until calcs added
            self.flight_angle.append(self.flight_angle[0])      # initial values until calcs added
            self.thrust_angle.append(self.thrust_angle[0])
            self.Cd.append(self.Cd[0])
            if self.runIter <= self.burntime * 2:
                self.dynamic_pressure.append(self.maxQ(rho))

            # END CONDITIONS
            if (self.altitude[self.runIter] < 1000 and self.time[self.runIter] > self.burntime) or self.time[self.runIter] > 10000:
                break

            self.runIter += 1
        return (self.altitude, self.velocity, self.acceleration, self.mass, self.time, self.thrust, self.drag, self.dynamic_pressure, self.rho, self.temp)

    def calc_Cd(self, M):
        return .75

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

    # standard atmosphere model (SI units)
    def STDATM(self, altitude):
        layer = -1.0            # gradient layer
        gradient = -0.00357
        altitude_base = 0.0
        temperature_base = 518.688
        density_base = 0.00237689241

        if altitude > 36089.239:
            layer = 1.0       # isothermal layer
            altitude_base = 36089.239
            temperature_base = 389.988
            density_base = 0.000707828857
        if altitude > 82020.997:
            layer = -1.0      # gradient layer
            gradient = 0.00165
            altitude_base = 82020.997
            temperature_base = 389.988
            density_base = 7.88546183 * (10**(-5))
        if altitude > 154199.48:
            layer = 1.0       # isothermal layer
            altitude_base = 154199.48
            temperature_base = 508.788
            density_base = 2.86391281 * (10**(-6))
        if altitude > 173884.51:
            layer = -1.0      # gradient layer
            gradient = -0.00246
            altitude_base = 173884.51
            temperature_base = 508.788
            density_base = 1.47056878 * (10**(-6))
        if altitude > 259186.35:
            layer = 1.0       # isothermal layer
            altitude_base = 259186.35
            temperature_base = 298.188
            density_base = 4.34631754 * (10**(-8))
        if altitude > 295275.59:
            layer = -1.0      # gradient layer
            gradient = 0.00219
            altitude_base = 295275.59
            temperature_base = 298.188
            density_base = 4.50154317 * (10**(-9))
        if layer < 0.0:
            temperature = temperature_base + gradient*(altitude - altitude_base)
            power = -1.0*(self.g0/gradient/self.R_air + 1.0)
            density = density_base*(temperature/temperature_base)**power
        else:
            temperature = temperature_base
            power = -1.0*self.g0*(altitude - altitude_base)/self.R_air/temperature
            density = density_base*np.exp(power)
        sos = np.sqrt(self.gamma_air*self.R_air*temperature)

        return (temperature, density, sos)


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
        'bank_angle': 0
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


if __name__ == '__main__':
    test_Rocket()
