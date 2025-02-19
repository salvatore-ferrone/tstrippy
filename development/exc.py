import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from temp import temp 


def king_density_function_of_W(W):
    sqrtW = np.sqrt(W)
    return np.exp(W)*custom_erf(sqrtW) - (2.0/np.sqrt(np.pi))*sqrtW*(1.0 + 2.0*W/3.0 )
    # return np.exp(W)*special.erf(sqrtW) - (2.0/np.sqrt(np.pi))*sqrtW*(1.0 + 2.0*W/3.0 )

def king_scale_radius(W0):
    rho0 = king_density_function_of_W(W0)
    return np.sqrt(9/ (4*np.pi*rho0))

def king_break_radius(W0):
    r0=king_scale_radius(W0)
    if W0 < 5.0:
        rbreak = r0/100
    else:
        rbreak = r0
    rbreak = r0
    return rbreak

def king_ode_in_r(r, y):
    W = y[0]
    dwdr = y[1]
    rho = king_density_function_of_W(W)
    if r <= 0: # numerical gaurd
        d2wdr2 = 0.0
    else:
        d2wdr2 = -4*np.pi*rho - 2*dwdr/r 
    return dwdr, d2wdr2

def king_ode_in_w(t, y):
    """ this one is weird, w is the independent variable, r is the dependent variable, and d/dr(INV(dw/dr)) is the second ode """ 
    W = t 
    r= y[0]
    dwdr = y[1]
    rho = king_density_function_of_W(W)
    term1 =-(1.0/dwdr)
    term2 = (4*np.pi*rho + 2.0*dwdr/r)
    return 1.0/dwdr, term1*term2

# Using a numerical approximation for the error function
def custom_erf(x):
    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # Save the sign of x
    sign = 1 if x >= 0 else -1
    x = abs(x)

    # A&S formula 7.1.26
    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * np.exp(-x * x)
    return sign * y

def solve_king_density_profile(W0,npoints=1001):
    # initialize 
    r,W,dWdr = np.zeros(npoints),np.zeros(npoints),np.zeros(npoints)
    rbreak=king_break_radius(W0)
    r[:npoints//2] = np.linspace(0,rbreak,npoints//2)
    limits = [0,rbreak]
    initial_conditions = [W0,0]
    sol = integrate.solve_ivp(king_ode_in_r, limits, initial_conditions, t_eval=r[:npoints//2],method="RK45",)
    # extract the first half of the solution 
    W[:npoints//2] = sol.y[0]
    dWdr[:npoints//2] = sol.y[1]
    Wmid = sol.y[0,-1]
    # now, solve for r given W
    W[npoints//2-1:] = np.linspace(Wmid,0,npoints-npoints//2 + 1)
    limits = [Wmid,0]
    initial_conditions = [rbreak,sol.y[1,-1]]
    sol = integrate.solve_ivp(king_ode_in_w, limits, initial_conditions, t_eval=W[npoints//2:],method="RK45",)
    r[npoints//2 :] = sol.y[0]
    dWdr[npoints//2  :] = sol.y[1]
    return r, np.array([W, dWdr])


if __name__=="__main__":
    npoints=1000
    nprint=3
    W0=11
    # r,w,dwdr=temp.solve_king_potential_profile(W0,npoints)
    temp.initialize_king_potential_profile(W0,npoints)
    r = temp.r_.copy()
    w = temp.w_.copy()
    dwdr = temp.dwdr_.copy()

    r_break = r[npoints//2]
    fig,axis=plt.subplots()
    axis.plot(r,w,'.')
    axis.stem([r_break],[w[npoints//2]],'r')
    # axis.plot(r,dwdr,'.')
    axis.set_xlabel('r')
    axis.set_title('King Potential Profile')
    plt.show()
    fig.savefig('king_potential_profile.png')