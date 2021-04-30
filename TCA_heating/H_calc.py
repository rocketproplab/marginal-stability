import numpy as np
import matplotlib.pyplot as plt

def area(r):
    return np.pi * (r**2)

def Reynolds(mv, r, mu):
    return mv*r/mu

def mv(r,mdot):
    a = area(r)
    return mdot/a

def Prandtl(c,mu,k):
    return c * mu / k

def Nusselt(r,mdot,c,mu,k,f):
    mr = mv(r,mdot)
    Re = Reynolds(mr,r,mu)
    Pr = Prandtl(c,mu,k)
    num = (f/8) * (Re - 1000) * Pr
    den = 1 + 12.7*np.sqrt(f/8) * (Pr**(2/3) - 1)
    return num/den

def h(r,mdot,c,mu,k,f):
    nu = Nusselt(r,mdot,c,mu,k,f)
    return nu*k/r

class gas():
    def __init__(self,m0,t0,m1,t1,c,k0):
        self.m0 = m0
        self.m1 = m1
        self.t0 = t0
        self.t1 = t1
        self.c = c
        self.k0 = k0
        self.fk = k0/(m0 * c)
        self.s = (t0*m1*np.sqrt(t0/t1)/m0 - t0)/(1 - (m1/m0)*np.sqrt(t0/t1)**3)

    def mu(self,T):
        return self.m0*((T/self.t0)**(3/2))*(self.t0+self.s)/(T + self.s)

    def k(self,T):
        m = self.mu(T)
        return m*self.fk*self.c

co2 = gas(41.47*10**(-6),1050,83.1*10**(-6),3225,1524,0.0737)
co = gas(29.24*10**(-6),600,78.2*10**(-6),3225,1334,0.0457)
h2 = gas(18.77*10**(-6),1000,35.35*10**(-6),3225,18790,0.4934)
h2o = gas(37*10**(-6),973,14.835*10**(-6),440,3133,0.0923)

gases = [co2,co,h2,h2o]
ms_weights = np.array([0.118,0.394,0.158,0.330])
phx_weights = np.array([0.0782,0.2632,0.1977,0.4609])


def avgmu(T,gaslist,weightlist):
    mus = np.array([gas.mu(T) for gas in gaslist])
    return np.sum(mus*weightlist)

def avgk(T,gaslist,weightlist):
    ks = np.array([gas.k(T) for gas in gaslist])
    return np.sum(ks*weightlist)

def avgc(gaslist,weightlist):
    cs = np.array([gas.c for gas in gaslist])
    return np.sum(cs*weightlist)

def avgh(r,T,gases,weights,mdot,f):
    c = avgc(gases,weights)
    mu = avgmu(T,gases,weights)
    k = avgk(T,gases,weights)
    return h(r,mdot,c,mu,k,f)

#print(avgh(0.0313436,2921,gases,ms_weights,8.05307894))

''' PHX hot gas extrapolations
print('CO2 ',co2.k(3225))
print('CO ',co.k(3225))
print('H2 ',h2.k(3225))
print('H2O ',h2o.k(3225))
print()
print('H2O mu',h2o.mu(3225))
print()
print('mu ',avgmu(3225,gases,phx_weights))
print('cp ',avgc(gases,phx_weights))
print('k ',avgk(3225,gases,phx_weights))
'''