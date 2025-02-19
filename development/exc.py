import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy import integrate
from temp import temp 
import galpy.df


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


def leapfrogstep(x,y,z,vx,vy,vz,dt):
    """Leapfrog step"""
    # kick 
    ax0,ay0,az0,phi0=temp.king_unscaled([W0],x,y,z)
    # drift 
    xf = x + 0.5*vx*dt
    yf = y + 0.5*vy*dt
    zf = z + 0.5*vz*dt
    # kick
    axf,ayf,azf,phi=temp.king_unscaled([W0],xf,yf,zf)
    # update the velocities
    vxf = vx + 0.5*(ax0+axf)*dt
    vyf = vy + 0.5*(ay0+ayf)*dt
    vzf = vz + 0.5*(az0+azf)*dt
    
    # update the velocities
    return xf,yf,zf,vxf,vyf,vzf

if __name__=="__main__":
    npoints=1000
    nprint=3
    W0=11
    # r,w,dwdr=temp.solve_king_potential_profile(W0,npoints)
    temp.initialize_king_potential_profile(W0,npoints)
    r = temp.r_.copy()
    w = temp.w_.copy()
    dwdr = temp.dwdr_.copy()
    r_tidal = temp.tidal_radius_
    # sample from bovy
    myking = galpy.df.kingdf(W0)
    R,vR,vT,z,vz,phi=myking.sample(n=100,return_orbit=False)
    xp,yp,zp = R*np.cos(phi),R*np.sin(phi),z
    vxp,vyp,vzp = vR*np.cos(phi) - vT*np.sin(phi),vR*np.sin(phi) + vT*np.cos(phi),vz
    ax,ay,az,phi=temp.king_unscaled([W0],xp,yp,zp)    
    mag = np.sqrt(ax**2 + ay**2 + az**2)
    speed = np.sqrt(vxp**2 + vyp**2 + vzp**2)
    norm = mpl.colors.LogNorm(vmin=mag.min(),vmax=mag.max())
    norm_speed = mpl.colors.Normalize(vmin=0,vmax=speed.max())
    cmap = plt.cm.plasma
    cmap_speed = plt.cm.plasma_r
    colors = cmap(norm(mag))
    colors = cmap_speed(norm_speed(speed))
    fig,axis=plt.subplots(1,1,figsize=(5,6))
    dt=1e-3
    NSTEP = 1000
    NSKIP = 10
    c=0
    for i in range(NSTEP):
        xp,yp,zp,vxp,vyp,vzp = leapfrogstep(xp,yp,zp,vxp,vyp,vzp,dt)
        ax,ay,az,phi=temp.king_unscaled([W0],xp,yp,zp)
        mag = np.sqrt(ax**2 + ay**2 + az**2)
        speed = np.sqrt(vxp**2 + vyp**2 + vzp**2)

        if i%NSKIP == 0:
            axis.cla()
            axis.set_aspect("equal")
            axis.set_xlim(-r_tidal/2,r_tidal/2)
            axis.set_ylim(-r_tidal/2,r_tidal/2)
            axis.set_xlabel("x [kpc]")
            axis.set_ylabel("y [kpc]")
            axis.set_title("King Model with W0={}".format(W0))
            axis.quiver(xp,yp,ax/mag,ay/mag,color=colors)
            axis.quiver(xp,yp,vxp/speed,vyp/speed,color="k",alpha=0.5)
            fig.tight_layout()
            plt.show()
            fig.savefig("frame-{:d}.png".format(c))
            c+=1