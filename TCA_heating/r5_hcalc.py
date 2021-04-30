from H_calc import *
import csv
import scipy.interpolate as inter
import matplotlib.pyplot as plt
import numpy as np

locvec = []
with open('rev5_props.csv') as file:
    csv_reader = csv.reader(file, delimiter=',')
    rownum = 0
    for row in csv_reader:
        rownum += 1
        if rownum > 1:
            locvec.append((float(row[0]),float(row[1]),float(row[2])))

sortedvec = sorted(locvec,key=lambda coor: coor[1])
shortvec = []
ylist = []
for loc in sortedvec:
    if loc[1] not in ylist:
        shortvec.append(loc)
    ylist.append(loc[1])
yvec = [tup[1] for tup in shortvec]
rvec = [(tup[0]**2 + tup[2]**2)**(1/2) for tup in shortvec]
r = inter.interp1d(yvec,rvec,fill_value='extrapolate')

roughness = 1.98*10**(-6)

def T(x):
    return  1353.133 + (3469.594 - 1353.133)/(1 + (x/8.790652)**20.48152)**0.2413815

def x(y):
    return (-y+0.228)*100/2.54

def Ty(y):
    return T(x(y))

def region_avg(func,start,end):
    space = np.linspace(start,end,1000)
    vals = [func(i) for i in space]
    return round(np.average(vals))

mdot = 8.05 #For MS

'''avgr = region_avg(r,-0.107,0.226)
print(avgr)
avga = area(avgr)
avgmv = mv(avgr,mdot)
print(avgmv)
avgT = region_avg(Ty,-0.107,0.226)
print(avgT)
avgmuval = avgmu(avgT,gases,ms_weights)
print(avgmuval)
avgre = Reynolds(avgmv,avgr,avgmuval)
ed = roughness/(2*avgr)
print(avgre)
print(ed)'''

f = 0.009 
# from moody diagram https://www.researchgate.net/profile/Marco_Ferro3/publication/296700902/figure/fig1/AS:335833400922112@1457080330658/Moodys-diagram-depicting-the-friction-factor-in-function-of-Reynolds-number.png
# using commented code above to get inputs

def hcalc(y):
    return avgh(r(y),Ty(y),gases,ms_weights,mdot,f)

regionlist_r5 =   [(-.107,-0.081),
                (-0.081,-0.056),
                (-0.056,-0.030),
                (-0.030,-0.010),
                (-0.010,0.000),
                (0.000,0.026),
                (0.026,0.042),
                (0.042,0.064),
                (0.064,0.100),
                (0.100,0.138),
                (0.138,0.187),
                (0.187,0.215),
                (0.215,0.221),
                (0.221,0.226)]

regionlist_r6 =   [(-.107,-0.094),
                (-0.094,-0.076),
                (-0.076,-0.041),
                (-0.041,-0.008),
                (-0.008,0.011),
                (0.011,0.026),
                (0.026,0.041),
                (0.041,0.072),
                (0.072,0.097),
                (0.097,0.130),
                (0.130,0.173),
                (0.173,0.214),
                (0.214,0.219),
                (0.219,0.226)]
regionlist_r6.reverse()

regionlist_r7 =   [(-.107,-0.075),
                (-0.075,-0.035),
                (-0.035,-0.005),
                (-0.005,0.001),
                (0.001,0.027),
                (0.027,0.045),
                (0.045,0.074),
                (0.074,0.113),
                (0.113,0.157),
                (0.157,0.214),
                (0.214,0.221),
                (0.221,0.225),
                (0.225,0.226)]
regionlist_r7.reverse()

regionlist_r8 = [(.226,.225),
                (.225,.221),
                (.221,.218),
                (.218,.216),
                (.216,.162),
                (.162,.077),
                (.077,.045),
                (.045,.026),
                (.026,.004),
                (.004,-.005),
                (-.005,-.090),
                (-.090,-.109),
                (-.109,-.112),
                (-.112,-.115)]


regionlist_r8_2 = [(0.226,    0.226),
(0.226,   0.222),
(0.222,   0.218),
(0.218,   0.216),
(0.216,   0.116),
(0.116,   0.075),
(0.075,   0.042),
(0.042,   0.026),
(0.026,   0.001),
(0.001,   -0.006),
(-0.006,  -0.109),
(-0.109,  -0.112),
(-0.112,  -0.115)]

regionlist_r8_3 = [(0.228,    0.22),
(0.22,    0.216),
(0.216,   0.079),
(0.079,   0.042),
(0.042,   0.027),
(0.027,   0.001),
(0.001,   -0.006),
(-0.006,  -0.069),
(-0.069,  -0.108),
(-0.108,  -0.112),
(-0.112,  -0.114),
(-0.114,  -0.115)
]

for region in regionlist_r6:
    print(region_avg(Ty,region[1],region[0]))

print()

for region in regionlist_r6:
    print(region_avg(hcalc,region[1],region[0]))

'''lst = np.linspace(-.111,.228,1000)
h = [hcalc(i) for i in lst]
plt.plot(lst,h)
plt.show()'''