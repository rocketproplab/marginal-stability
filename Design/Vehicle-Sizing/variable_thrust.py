import numpy as np

def natural_cubic_spline(xi, yi, xtarget):
    np.testing.assert_approx_equal(len(xi), len(yi), err_msg='lengths are not equal')
    diff = (xi[j+1] - xi[j] for j in range(len(xi) - 1))
    for i in diff:
        if i <= 0:
            raise Exception('x values are not ordered')

    n = len(xi)
    M = np.zeros(n)
    rhs = np.zeros(n,1)
    hi = diff
    for i in range(1,n-1):
        M[i,i-1] = hi[i-1]
        M[i,i] = 2*(hi[i-1]+hi[i])
        M[i,i+1] = hi[i]
        rhs[i] = 3*(((yi[i+1] -yi[i])/(xi[i+1] - xi[i])) - ((yi[i] - yi[i-1])/(xi[i] - xi[i-1])))

    M[1,1] = 1
    M[1,2] = 0
    M[n,n] = 1
    M[n,n-1] = 0
    rhs[1] = 0
    rhs[n] = 0

    a = np.zeros(n-1)
    b = np.zeros(n-1)
    d = np.zeros(n-1)
    C = np.linalg.solve(M,rhs)
    for i in range(n-1):
        a[i] = yi[i]
        d[i] = (c[i+1] - c[i])/(3*hi[i])
        b[i] = (yi[i+1] - yi[i]) / (hi[i]) - (c[i+1]+2*c[i])*hi[i]/3


    slist = []
    for i in range(n-1):
        f = lambda x : a[i] + b[i]*(x - xi[i]) + c[i]*(x - xi[i])**2 + d*(x - xi[i])**3
        slist.append(f)
    S = np.array(slist)

    for i in range(n-1):
        if xtarget >= xi[i] and xtarget <= xi[i+1]:
            ytarget = S[i](xtarget)
        if xtarget > xi[-1]:
            ytarget = S[i](xi[-1])
    return ytarget

def Thrust_Chamber_Design(Thrust,Impulse):
    Pc = 650                # psi
    Pc_atm = Pc/14.6959
    P_Launch = 0.849                   # atm
    Pe = .5*P_Launch               # atm ("Summerfield Criterion")
    P_amb = P_Exit                     # atm
    Burn_Time = Impulse/Thrust          # sec

    dPe = [34, 33, 34, 33, 33, 33, 34, 33, 34, 34, 34, 34, 34, 34, 34, 34, 34,
        33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33,
        33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33]
    hPe = [1, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17,
        18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
        30, 31, 32, 33, 33, 34, 35, 36, 37, 38, 39, 39,
        40, 41, 42, 43, 43, 44, 45, 45, 46, 47, 47]

    Rmin = 2.20                          # O/F            lowest value of mixture ratio
    P = Pc_atm - 25                      # atm            atm of pressure ofer 25 atm
    dPc_R = 25/49                        # atm/pixel      x-axis increments
    n_pix_Pc_R = P//dPc_R          # pixels         # of pixels past 25 atm on x-axis
    if n_pix_Pc_R < 1:
        n_pix_Pc_R = 1
    dR = .45/341                         # units/pixel    y-axis increments
    n_dPe = dPe[n_pix_Pc_R - 1]              # pixels         dPe for desired chamber pressure
    n_pix_Pe = n_dPe*(1-1/.9*(Pe-.1))
    Mixture_Ratio = Rmin + dR*(n_pix_Pe + hPe[n_pix_Pc_R - 1])

    dRT = [30, 30, 30, 30, 30, 30, 31, 31, 32, 32, 32, 32,    # number of pixels between Ratio curves (2.2 and 2.3)
       33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 34, 34,
       34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 35, 35,
       34, 35, 35, 35, 36, 35, 36, 36, 35, 36, 36, 36]

    hRT = [10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,     # number of pixels between Tmin and bottom curve
       22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 31, 32,
       33, 34, 35, 35, 36, 37, 37, 38, 39, 40, 40, 41,
       42, 42, 43, 44, 44, 45, 45, 46, 47, 47, 48, 48]

    Tmin = 3350                            # K              lowest value of temperature
    dPc_T = 25/48                          # atm/pixel      x-axis increments
    n_pix_Pc_T = P//dPc_T            # pixels         # of pixels past 25 atm on x-axis
    if n_pix_Pc_T < 1:
        n_pix_Pc_T = 1
    dT = 600/341                           # K/pixel        y-axis increments
    n_dRT = dRT[n_pix_Pc_T]                # pixels         dRT for desired chamber pressure
    n_pix_RT = n_dRT/.1*(R-2.2)            # pixels         # of pixels above bottom curve

    T_Chamber = Tmin + dT*(n_pix_RT + hRT[n_pix_Pc_T])  # Kelvin

    dRM = [74, 74, 74, 74, 74, 75, 75, 74, 75, 75, 75, 75,    # number of pixels between Ratio curves (2.2 and 2.3)
       75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75,
       75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 76, 76,
       76, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77]

    hRM = [32, 33, 33, 34, 35, 35, 36, 37, 37, 38, 38, 39,     # # of pixels between Mmin and bottom curve
       40, 40, 41, 41, 42, 43, 43, 44, 44, 45, 45, 46,
       46, 47, 47, 47, 48, 48, 49, 49, 49, 50, 50, 50,
       50, 51, 51, 51, 52, 52, 52, 53, 53, 53, 54, 54]

    Mmin = 20.8                          # amu            lowest value of weight
    dPc_M = 25/48                        # atm/pixel      x-axis increments
    n_pix_Pc_M = P//dPc_M          # pixels         # of pixels past 25 atm on x-axis
    if n_pix_Pc_M < 1:
        n_pix_Pc_M = 1
    dM = 2/341                            # amu/pixel      y-axis increments
    n_dRM = dRM[n_pix_Pc_M]               # pixels         dRT for desired chamber pressure
    n_pix_RM = n_dRM/.1*(R-2.2)           # pixels         # of pixels above bottom curve

    Gas_Molecular_Weight = Mmin + dM*(n_pix_RM + hRM[n_pix_Pc_M])  # amu

    dRy = [46, 46.5, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47,    # number of pixels between Ratio curves (2.2 and 2.3)
    47, 47, 47, 47, 47.5, 48, 48, 48, 48, 48, 48, 48,
    48, 48, 48.5, 48.5, 48, 48.5, 49, 49, 49, 49, 49, 49,
    49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49]

    hRy = [58, 56.5, 55, 54, 53, 52, 50.5, 49, 48, 47, 46, 45.5,     # number of pixels between ymin and bottom curve
        45, 44, 43, 42, 41, 40, 39, 38, 37, 36.5, 36, 35,
        34, 33, 32, 31.5, 31, 30, 29, 28, 27.5, 27, 26, 25.5,
        25, 24, 23.5, 23, 22.5, 22, 21, 20.5, 20, 19.5, 19]

    ymin = 1.219                         # amu            lowest value of specific heat ratio
    dPc_y = 25/47                        # atm/pixel      x-axis increments
    n_pix_Pc_y = P//dPc_y          # pixels         # of pixels past 25 atm on x-axis
    if n_pix_Pc_y < 1:
        n_pix_Pc_y = 1

    dy = .021/341                         # amu/pixel      y-axis increments
    n_dRy = dRy[n_pix_Pc_y]               # pixels         dRy for desired chamber pressure
    n_pix_Ry = n_dRy*(1-1/.1*(R-2.2))     # pixels         # of pixels above bottom curve

    y = ymin + dy*(n_pix_Ry + hRy[n_pix_Pc_y])  # Cp/Cv

    g = 32.1740*12                     # in/s^2
    R = 1545.348963*12                 # in-lbf/lb-mol-R
    FS = 1.75                          # Safety Factor
    Stress = 115000                    # psi High temperature yield stress of inconel walling (http://www.specialmetals.com/assets/smc/documents/inconel_alloy_718.pdf)
    L_star = 45                        # in

    P_Chamber = P_Chamber_psi/14.6959      # atm
    T_Chamber_R = 1.8*T_Chamber            # Rankine

    lamda = .983
    e = ( ((2/(y+1))**(1/(y-1))) * (P_Chamber/P_Exit)**(1/y) ) / np.sqrt(((y+1)/(y-1)) * (1-(P_Exit/P_Chamber)**((y-1)/y)) )
    C_star = lamda*np.sqrt( g*y*R*T_Chamber_R/Gas_Molecular_Weight )/( y*np.sqrt( (2/(y+1))**((y+1)/(y-1)) ) )
    Cf = lamda*(np.sqrt( ((2*y**2)/(y-1)) * ((2/(y+1))**((y+1)/(y-1))) * (1-(P_Exit/P_Chamber)**((y-1)/y) )) + E*((P_Exit-P_amb)/P_Chamber))
    Isp = lamda*C_star*Cf/g
    M = Thrust/Isp;                # lbs/s
    M_RP1 = M/(1+Mixture_Ratio);   # lbs/s
    M_LOX = M - M_RP1;             # lbs/s

    Propellant_Mass = M*Burn_Time # lbs
    LOx_Mass = M_LOX*Burn_Time    # lbs
    RP1_Mass = M_RP1*Burn_Time    # lbs