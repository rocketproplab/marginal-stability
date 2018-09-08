''' Beginning of Marginal Stability Tank Sizing Code
    This code will take inputs from a .txt file containing all relevant
    vehicle sizing data and information and compute important estimates
    such as tank height, volume, thickness, and mass.  These values
    provided will only be estimates; however, the code will take into
    account as many variables and parameters as possible to ensure
    the most accurate result.
    Work in Progress

'''
# Establish needed variables and parameters for calculations
# This will all eventually be handled by a .txt file, when I learn how to do that :)
import math # gives me all the available constants I need

class tankSize:
    def __init__(self, rp1Mass, loxMass, diameter, pressure, factorOfSafety, pressurantDiameter, initialPressure, finalPressure):
        self.rp1Mass = rp1Mass
        self.loxMass = loxMass
        self.diameter = diameter
        self.pressure = pressure
        self.factorOfSafety = factorOfSafety
        self.pressurantDiameter = pressurantDiameter
        self.initialPressure = initialPressure
        self.finalPressure = finalPressure

        pi = math.pi
        rho_lox = 0.04122124 # density of LOx as given in Vehicle Sizing spreadsheet V1, lbm/in^3
        rho_rp1 = 0.0292631 # density of RP-1, MS Vehicle Sizing V1
        rho_alum = 0.0975 # density of Aluminum 6061-T6, the likely choice for tanks
        ys_alum = 40030.488 # yield strength of 6061-T6


        # Begin importing text file stuff...yay!
        #file = open('input.txt', 'r')
        # this opens the file 'input.txt' from whatever directory the program is
        # saved to--make sure to download it or move it to the same directory!

        #prop_diam_text = file.readlines(1) # reads first line of the input file
        #for line in prop_diam_text:
        #    prop_diam_stuff = line.split() # saves split results of first line
        #diam = float(prop_diam_stuff[1]) # saves proper variables from that array, converted to float
        diam = diameter

        #prop_press_text = file.readlines(2) # continues to read down the input file
        #for line in prop_press_text:
        #    prop_press_stuff = line.split()
        #    print(prop_press_stuff)
        #lox_press = float(prop_press_stuff[1])
        lox_press = pressure
        rp1_press = lox_press # can be better worked into the code later, but should be the same
        ''' Quick note here: I'd already written the code such that it would allow for two different
        tank pressure inputs before adding the input text file capability, and I really didn't want
        to go back and rewrite all my code at the risk of missing a variable somewhere or skewing an
        equation.  Leaving it like this also makes it a lot easier down the road if we decide we'd like
        to keep the two tanks at different pressures for some reason.  It also makes it a lot easier to
        tell which tanks you're doing calculations for if they have different variable names, rather than
        just a generic "tank pressure".'''

        #fos_text = file.readlines(3)
        #for line in fos_text:
        #    fos_stuff = line.split()
        #fos = float(fos_stuff[1])
        fos = factorOfSafety

        #press_diam_text = file.readlines(4)
        #for line in press_diam_text:
        #    press_diam_stuff = line.split()
        #press_tank_diam = float(press_diam_stuff[1])
        press_tank_diam = pressurantDiameter

        #press_ipress_text = file.readlines(5)
        #for line in press_ipress_text:
        #    press_ipress_stuff = line.split()
        #initial_press = float(press_ipress_stuff[1])
        initial_press = initialPressure

        #press_fpress_text = file.readlines(6)
        #for line in press_fpress_text:
        #    press_fpress_stuff = line.split()
        #final_press = float(press_fpress_stuff[1])
        final_press = finalPressure

        #file.close() # closes the input file once everything's read in

        k = 0.725158 # Roundness factor for 2:1 elliptical heads, as given by [insert link here]
        E = 0.85 # joint efficiency of any welded joints, might change or not be used
        m_lox = self.loxMass # Mass of LOx required, lbm
        m_rp1 = self.rp1Mass # Mass of RP-1 required, lbm
        of_ratio = 2.33 # oxidizer:fuel ratio, as given currently in MSVS V1
        #lox_press = 600 # Max. Expected Operating Pressure of LOx, psi
        #rp1_press = 600 # MEOP of RP-1, psi
        tubediam = 0.83 # OD for transfer tube, in (for subtracting volume)

        vol_lox = m_lox / rho_lox # computes the volume of the propellants to be used, in^3
        vol_rp1 = m_rp1 / rho_rp1

        # Barlow's equation for tank wall thickness
        t_lox = (lox_press * diam) / (2 * ys_alum / fos)
        t_rp1 = (lox_press * diam) / (2 * ys_alum / fos)

        # Tank head thickness equations (taken from [insert link here])
        th_lox = (lox_press * diam * k) / ((2 * E * ys_alum / fos) - (0.2 * lox_press))
        th_rp1 = (rp1_press * diam * k) / ((2 * E * ys_alum / fos) - (0.2 * lox_press))

        # Inner volume of tank heads (neutrium eqns)
        vh_lox = (diam - (2 * th_lox))**3 * 0.5 * (pi / 24) * 2

        vh_rp1 = (diam - (2 * th_rp1))**3 * 0.5 * (pi / 24) * 2

        # Tank straight-wall height requirements
        sh_lox = (vol_lox - (2 * vh_lox)) / (pi * ((diam - 2 * t_lox)/2)**2)
        sh_rp1 = vol_rp1 / (pi * (((diam - (2 * t_rp1))/2)**2 - (tubediam / 2)**2)) # factors out internal plumbing

        # Volume of LOx tank
        vol_loxtank = (2 * vh_lox) + (pi * ((diam - (2 * t_lox))/2)**2 * sh_lox)
        vol_rp1tank = (pi * ((diam - (2 * t_rp1))/2)**2 * sh_rp1)

        # Tank head heights (for CAD purposes, mostly)
        hi_lox = (diam - 2 * th_lox) / 4 # inner height
        ho_lox = hi_lox + th_lox # outer height
        hi_rp1 = (diam - 2 * th_rp1) / 4
        ho_rp1 = hi_rp1 + th_rp1

        # Tank heights, separate and together
        lox_height = 2 * ho_lox + sh_lox
        rp1_height = ho_rp1 + sh_rp1
        total_height = ho_rp1 + ho_lox + sh_rp1 + sh_lox
        '''  This isn't just LOx height + RP-1 height because the tanks share a common bulkhead.
            When installed, the bottom of the LOx tank actually sits inside the shoulder of the
            RP-1 tank, so the overall length of the tanks sittingg in the rocket will be the
            height of the upper LOx bulkhead + the straight-wall LOx height + the RP-1
            straight-wall height (little/no bulkhead space in between) + the height of the
            lower RP-1 bulkhead
            Later versions of this code should acknowledge that there is, in fact, a bit of
            space between the two tanks, especially if they are double-walled.  This can be
            updated as the design is finalized
        '''


        # Numbers for CAD ease and tank designing specifications
        id_lox = diam - 2 * t_lox # inside diameter of the tanks, in
        id_rp1 = diam - 2 * t_rp1
        KRi_lox = 0.17 * id_lox # knuckle radii and crown radii of ellipses, inside and out, in
        CRi_lox = 0.9 * id_lox
        KRo_lox = 0.17 * diam
        CRo_lox = 0.9 * diam
        KRi_rp1 = 0.17 * id_rp1
        CRi_rp1 = 0.9 * id_rp1
        KRo_rp1 = 0.17 * diam
        CRo_rp1 = 0.9 * diam

        # Tank metal volume and mass estimates
        loxheadmetal = ((diam**3) * 0.5 * (pi / 24) * 2) - vh_lox
        rp1headmetal = ((diam**3) * 0.5 * (pi/24) * 2) - vh_rp1
        metalwalllox = (pi * (diam / 2)**2 * sh_lox) - (pi * (((diam - 2 * t_lox)/2)**2) * sh_lox)
        metalwallrp1 = (pi * (diam / 2)**2 * sh_rp1) - (pi * (((diam - 2 * t_rp1)/2)**2) * sh_rp1)
        metalvolume = 2*loxheadmetal + metalwalllox + metalwallrp1 + 2 * rp1headmetal

        self.tank_mass = rho_alum * metalvolume
        head_mass = (2 * loxheadmetal + 2 * rp1headmetal) * rho_alum
        wall_mass = (metalwalllox + metalwallrp1) * rho_alum

        file2 = open('output.txt', 'w') # creates generic output file, can be renamed later
        # Print results to output file rather than command line
        file2.write('The propellant tank diameter is ' + str(diam) + ' in\n')
        file2.write('The LOx tank pressure is ' + str(lox_press) + ' psi\n')
        file2.write('The RP-1 tank pressure is ' + str(rp1_press) + ' psi\n')
        file2.write('The total mass is ' + str(round(self.tank_mass, 3)) + ' lbs\n')
        file2.write('The total tank height is ' + str(round(total_height, 3)) + ' in\n')
        file2.write('The LOx tank is ' + str(round(lox_height, 3)) + ' in high\n')
        file2.write('The RP-1 tank is ' + str(round(rp1_height, 3)) + ' in high\n')
        '''Another quick note: I floated everything to 3 digits here just for simplicity's sake.
        It won't hurt my feelings if you need to float the decimals further or just take
        rounding off altogether, if precision becomes more important than 3 decimals down
        the road'''


        ''' Moving on to the Pressurant Tank Sizing Section
        This part of the code will take in other values, probably from the same .txt
        file, and compute the needed mass of helium, as well as the volume of helium,
        and size a tank appropriately
        '''

        # Prepare other variables from propellant tank sizing
        #press_tank_diam = 6 # has to be user-defined, need a diameter to help find volume
        prop_vol = vol_lox + vol_rp1 # gives the total volume of propellants being ejected
        prop_press = lox_press # since they should both be the same anyways
        gamma = 5 / 3; # the ratio of specific heats for the gas equations
        #initial_press = 4500 # pressure inside helium tank on the pad
        #final_press = 50 # pressure inside tank after burnout
        press_temp = 536.67 # gives temperature of pressurant in *R
        R = 10.73159*12**3 # gas constant in imperial, in^3-psi / *R-lbmol
        molar_mass = 4.003 # also for gas equations, lbm/lbmol

        # Calculate the pressurant mass from equations given by AIAA report
        press_mass = ((prop_press * prop_vol) / (R * press_temp)) * (gamma / (1-(final_press/initial_press))) * molar_mass

        # Calculate volume of pressurant from ideal gas eqns
        press_vol = (press_mass * (1/molar_mass) * R * press_temp) / initial_press

        # Use the volume and diameter to find tank straight-wall height, then add height of spherical end caps
        tank_height = ((press_vol - (4/3) * pi * (press_tank_diam / 2)**3) / (pi * (press_tank_diam / 2)**2)) + press_tank_diam

        # ...and print the results!
        file2.write('The helium tank diameter is ' + str(press_tank_diam) + ' in\n')
        file2.write('The volume of helium needed is ' + str(round(press_vol, 3)) + ' in^3\n')
        file2.write('The mass of helium needed is ' + str(round(press_mass, 3)) + ' lbm\n')
        file2.write('The pressurant tank is ' + str(round(tank_height, 3)) + ' in high\n')

        file2.close()

example = tankSize(115.5001155,269.1152691,10,750,4,8,4500,50)
