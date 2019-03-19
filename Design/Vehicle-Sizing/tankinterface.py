import numpy as np
import tanksim

# Begin importing text file stuff...yay!
file = open('input.txt', 'r')
# this opens the file 'input.txt' from whatever directory the program is
# saved to--make sure to download it or move it to the same directory!

prop_diam_text = file.readlines(1) # reads first line of the input file
for line in prop_diam_text:
    prop_diam_stuff = line.split() # saves split results of first line
diam = float(prop_diam_stuff[1]) # saves proper variables from that array, converted to float

prop_press_text = file.readlines(2) # continues to read down the input file
for line in prop_press_text:
    prop_press_stuff = line.split()
lox_press = float(prop_press_stuff[1])
rp1_press = lox_press # can be better worked into the code later, but should be the same
''' Quick note here: I'd already written the code such that it would allow for two different
tank pressure inputs before adding the input text file capability, and I really didn't want
to go back and rewrite all my code at the risk of missing a variable somewhere or skewing an
equation.  Leaving it like this also makes it a lot easier down the road if we decide we'd like
to keep the two tanks at different pressures for some reason.  It also makes it a lot easier to
tell which tanks you're doing calculations for if they have different variable names, rather than
just a generic "tank pressure".'''

fos_text = file.readlines(3)
for line in fos_text:
    fos_stuff = line.split()
fos = float(fos_stuff[1])

press_diam_text = file.readlines(4)
for line in press_diam_text:
    press_diam_stuff = line.split()
press_tank_diam = float(press_diam_stuff[1])

press_ipress_text = file.readlines(5)
for line in press_ipress_text:
    press_ipress_stuff = line.split()
initial_press = float(press_ipress_stuff[1])

press_fpress_text = file.readlines(6)
for line in press_fpress_text:
    press_fpress_stuff = line.split()
final_press = float(press_fpress_stuff[1])

file.close() # closes the input file once everything's read in

testtank = tanksim.tank(diam,lox_press,fos,press_tank_diam,initial_press,final_press,796.618,357.228)
tankSpecs = testtank.tankCalc()
pressSpecs = testtank.pressTankCalc()

tank_mass = tankSpecs['tank_mass']
total_height = tankSpecs['total_height']
lox_height = tankSpecs['lox_height']
rp1_height = tankSpecs['rp1_height']

press_vol = pressSpecs['press_vol']
press_mass = pressSpecs['press_mass']
tank_height = pressSpecs['tank_height']

print('The propellant tank diameter is ' + str(diam) + ' in')
print('The LOx tank pressure is ' + str(lox_press) + ' psi')
print('The RP-1 tank pressure is ' + str(rp1_press) + ' psi')
print('The total mass is ' + str(round(tank_mass, 3)) + ' lbs')
print('The total tank height is ' + str(round(total_height, 3)) + ' in')
print('The LOx tank is ' + str(round(lox_height, 3)) + ' in high')
print('The RP-1 tank is ' + str(round(rp1_height, 3)) + ' in high')
print(tank_mass * lox_height / total_height)
print(tank_mass * rp1_height / total_height)

print('The helium tank diameter is ' + str(press_tank_diam) + ' in')
print('The volume of helium needed is ' + str(round(press_vol, 3)) + ' in^3')
print('The mass of helium needed is ' + str(round(press_mass, 3)) + ' lbm')
print('The pressurant tank is ' + str(round(tank_height, 3)) + ' in high')
