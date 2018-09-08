import numpy as np

class tank:
    def __init__(self,diam,press,fos,pressDiam,initial_press,final_press,m_lox,m_rp1):
        self.pi = np.pi
        self.rho_lox = 0.04122124 # density of LOx as given in Vehicle Sizing spreadsheet V1, lbm/in^3
        self.rho_rp1 = 0.0292631 # density of RP-1, MS Vehicle Sizing V1
        self.rho_alum = 0.0975 # density of Aluminum 6061-T6, the likely choice for tanks
        self.ys_alum = 40030.488 # yield strength of 6061-T6
        self.diam = diam
        self.lox_press = press
        self.rp1_press = press
        self.fos = fos
        self.press_tank_diam = pressDiam
        self.initial_press = initial_press
        self.final_press = final_press
        self.k = 0.725158 # Roundness factor for 2:1 elliptical heads, as given by [insert link here]
        self.Ef = 0.85 # joint efficiency of any welded joints, might change or not be used
        self.m_lox = m_lox
        self.m_rp1 = m_rp1
        self.tubediam = 0.83 # OD for transfer tube, in (for subtracting volume)
        self.lox_rectangle_len = 1.855
        self.rp1_rectangle_len = 1.05
        self.gamma = 5 / 3; # the ratio of specific heats for the gas equations
        self.press_temp = 536.67 # gives temperature of pressurant in *R
        self.gasConst = 10.73159*12**3 # gas constant in imperial, in^3-psi / *R-lbmol
        self.molar_mass = 4.003 # also for gas equations, lbm/lbmol

    def volume(self,mass,density):
        return mass/density

    def tankThickness(self,pressure,diameter,ys,fos):
        return (pressure * diameter) / (2 * ys / fos)

    def headThickness(self,pressure,diameter,k,Ef,ys,fos):
        return (pressure * diameter * k) / ((2 * Ef * ys / fos) - (0.2 * pressure))

    def mountRectangleVolume(self,rectangleLen):
        indVol = 0.5**2 * rectangleLen
        totVol = indVol * 6
        return totVol

    def headInnerVolume(self,diameter,headThickness,rectangleLen):
        totalVolume = (diameter - (2 * headThickness))**3 * 0.5 * (self.pi / 24) * 2
        return totalVolume - self.mountRectangleVolume(rectangleLen)

    def LOxTankWallHeight(self,loxVolume,diameter,headVolume,tankThickness):
        return (loxVolume - (2 * headVolume)) / (self.pi * ((diameter - 2 * tankThickness)/2)**2)

    def RP1TankWallHeight(self,rp1Volume,diameter,tankThickness):
        return (rp1Volume) / (self.pi * (((diameter - (2 * tankThickness))/2)**2 - (self.tubediam / 2)**2))

    def LOxTankVolume(self,headVolume,diameter,tankThickness,wallHeight):
        return (2 * headVolume) + (self.pi * ((diameter - (2 * tankThickness))/2)**2 * wallHeight)

    def RP1TankVolume(self,diameter,tankThickness,wallHeight):
        return (self.pi * ((diameter - (2 * tankThickness))/2)**2 * wallHeight)

    def intHeadHeight(self,diameter,headThickness):
        return (diameter - 2 * headThickness) / 4

    def extHeadHeight(self,intHeadHeight,headThickness):
        return intHeadHeight + headThickness

    def loxTankHeight(self,extHeadHeight,wallHeight):
        return 2 * extHeadHeight + wallHeight

    def rp1TankHeight(self,extHeadHeight,wallHeight):
        return extHeadHeight + wallHeight

    def totalTankHeight(self,rp1ExtHeadHeight,loxExtHeadHeight,rp1TankWallHeight,loxTankWallHeight):
        return rp1ExtHeadHeight + loxExtHeadHeight + rp1TankWallHeight + loxTankWallHeight

    def CADSpecs(self,diameter,rp1TankThickness,loxTankThickness):
        id_lox = diameter - 2 * loxTankThickness # inside diameter of the tanks, in
        id_rp1 = diameter - 2 * rp1TankThickness
        KRi_lox = 0.17 * id_lox # knuckle radii and crown radii of ellipses, inside and out, in
        CRi_lox = 0.9 * id_lox
        KRo_lox = 0.17 * diameter
        CRo_lox = 0.9 * diameter
        KRi_rp1 = 0.17 * id_rp1
        CRi_rp1 = 0.9 * id_rp1
        KRo_rp1 = 0.17 * diameter
        CRo_rp1 = 0.9 * diameter
        return {'id_lox': id_lox,
                'id_rp1': id_rp1,
                'KRi_lox': KRi_lox,
                'CRi_lox': CRi_lox,
                'KRo_lox': KRo_lox,
                'CRo_lox': CRo_lox,
                'KRi_rp1': KRi_rp1,
                'CRi_rp1': CRi_rp1,
                'KRo_rp1': KRo_rp1,
                'CRo_rp1': CRo_rp1}

    def headMetalVolume(self,diameter,headInnerVolume):
        return ((diameter**3) * 0.5 * (self.pi / 24) * 2) - headInnerVolume

    def wallMetalVolume(self,diameter,wallHeight,wallThickness):
        return (self.pi * (diameter / 2)**2 * wallHeight) - (self.pi * (((diameter - 2 * wallThickness)/2)**2) * wallHeight)

    def totalMetalVolume(self,loxHeadMetalVolume,rp1HeadMetalVolume,loxWallMetalVolume,rp1WallMetalVolume,loxMountRectangleVolume,rp1MountRectangleVolume):
        return 2*(loxHeadMetalVolume + rp1HeadMetalVolume) + loxWallMetalVolume + rp1WallMetalVolume + loxMountRectangleVolume + rp1MountRectangleVolume

    def totalPropVol(self, loxMass, rp1Mass, loxDen, rp1Den):
        return self.volume(loxMass,loxDen) + self.volume(rp1Mass,rp1Den)

    def pressMass(self,propPress,propVol,finalPress,initialPress):
        return ((propPress * propVol) / (self.gasConst * self.press_temp)) * (self.gamma / (1-(finalPress/initialPress))) * self.molar_mass

    def pressVol(self,pressMass,initialPress):
        return (pressMass * (1/self.molar_mass) * self.gasConst * self.press_temp) / initialPress

    def pressTankHeight(self,pressVol,pressTankDiam):
        return ((pressVol - (4/3) * self.pi * (pressTankDiam / 2)**3) / (self.pi * (pressTankDiam / 2)**2)) + pressTankDiam

    def pressTankCalc(self):
        prop_vol = self.totalPropVol(self.m_lox,self.m_rp1,self.rho_lox,self.rho_rp1)
        press_mass = self.pressMass(self.lox_press,prop_vol,self.final_press,self.initial_press)
        press_vol = self.pressVol(press_mass,self.initial_press)
        tank_height = self.pressTankHeight(press_vol,self.press_tank_diam)
        return {'press_mass': press_mass,
                'press_vol': press_vol,
                'tank_height': tank_height}

    def tankCalc(self):
        vol_lox = self.volume(self.m_lox,self.rho_lox)
        vol_rp1 = self.volume(self.m_rp1,self.rho_rp1)
        t_lox = self.tankThickness(self.lox_press,self.diam,self.ys_alum,self.fos)
        t_rp1 = self.tankThickness(self.rp1_press,self.diam,self.ys_alum,self.fos)
        th_lox = self.headThickness(self.lox_press,self.diam,self.k,self.Ef,self.ys_alum,self.fos)
        th_rp1 = self.headThickness(self.rp1_press,self.diam,self.k,self.Ef,self.ys_alum,self.fos)
        vh_lox = self.headInnerVolume(self.diam,th_lox,self.lox_rectangle_len)
        vh_rp1 = self.headInnerVolume(self.diam,th_rp1,self.rp1_rectangle_len)
        sh_lox = self.LOxTankWallHeight(vol_lox,self.diam,vh_lox,t_lox)
        sh_rp1 = self.RP1TankWallHeight(vol_rp1,self.diam,t_rp1)
        vol_loxtank = self.LOxTankVolume(vh_lox,self.diam,t_lox,sh_lox)
        vol_rp1tank = self.RP1TankVolume(self.diam,t_rp1,sh_rp1)
        hi_lox = self.intHeadHeight(self.diam,th_lox)
        ho_lox = self.extHeadHeight(hi_lox,th_lox)
        hi_rp1 = self.intHeadHeight(self.diam,th_rp1)
        ho_rp1 = self.extHeadHeight(hi_rp1,th_rp1)
        lox_height = self.loxTankHeight(ho_lox,sh_lox)
        rp1_height = self.rp1TankHeight(ho_rp1,sh_rp1)
        total_height = self.totalTankHeight(ho_rp1,ho_lox,sh_rp1,sh_lox)
        CAD_specs = self.CADSpecs(self.diam,t_rp1,t_lox)
        loxheadmetal = self.headMetalVolume(self.diam,vh_lox)
        rp1headmetal = self.headMetalVolume(self.diam,vh_rp1)
        metalwalllox = self.wallMetalVolume(self.diam,sh_lox,t_lox)
        metalwallrp1 = self.wallMetalVolume(self.diam,sh_rp1,t_rp1)
        lox_total_rect_volume = self.mountRectangleVolume(self.lox_rectangle_len)
        rp1_total_rect_volume = self.mountRectangleVolume(self.rp1_rectangle_len)
        metalvolume = self.totalMetalVolume(loxheadmetal,rp1headmetal,metalwalllox,metalwallrp1,lox_total_rect_volume,rp1_total_rect_volume)
        tank_mass = metalvolume * self.rho_alum
        return {'vol_lox': vol_lox,
                'vol_rp1': vol_rp1,
                't_lox': t_lox,
                't_rp1': t_rp1,
                'th_lox': th_lox,
                'th_rp1': th_rp1,
                'vh_lox': vh_lox,
                'vh_rp1': vh_rp1,
                'sh_lox': sh_lox,
                'sh_rp1': sh_rp1,
                'vol_loxtank': vol_loxtank,
                'vol_rp1tank': vol_rp1tank,
                'hi_lox': hi_lox,
                'ho_lox': ho_lox,
                'hi_rp1': hi_rp1,
                'ho_rp1': ho_rp1,
                'lox_height': lox_height,
                'rp1_height': rp1_height,
                'total_height': total_height,
                'CAD_specs': CAD_specs,
                'loxheadmetal': loxheadmetal,
                'rp1headmetal': rp1headmetal,
                'metalwalllox': metalwalllox,
                'metalwallrp1': metalwallrp1,
                'lox_total_rect_volume': lox_total_rect_volume,
                'rp1_total_rect_volume': rp1_total_rect_volume,
                'metalvolume': metalvolume,
                'tank_mass': tank_mass}
