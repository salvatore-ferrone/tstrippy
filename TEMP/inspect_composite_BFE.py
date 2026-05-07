import gravitymini
import numpy as np 
import matplotlib.pyplot as plt
import datetime 


def get_one_comp(coordPHI,coordQ,comp):
    grav = gravitymini.gravity
    grav.cleargravity()
    grav.addgravitycomponent(*comp)
    grav.finalizegravity()
    phi = grav.evaluategravitypotential(*coordPHI)
    starttime = datetime.datetime.now()
    ax,ay,az=grav.evaluategravityforces(*coordQ)
    endtime = datetime.datetime.now()
    comptime = endtime-starttime
    return ax,az,phi,comptime

def get_two_comp(coordPHI,coordQ,comp1,comp2):
    grav = gravitymini.gravity
    grav.cleargravity()
    grav.addgravitycomponent(*comp1)
    grav.addgravitycomponent(*comp2)
    grav.finalizegravity()
    phi = grav.evaluategravitypotential(*coordPHI)
    starttime = datetime.datetime.now()
    ax,_,az=grav.evaluategravityforces(*coordQ)
    endtime = datetime.datetime.now()
    comptime = endtime-starttime
    return ax,az,phi,comptime


## set the parameters 
rho_h_table = 11.4  # [Msun pc^-3]/1000
rho0_halo = (rho_h_table / 1000.0) * 1e9  # Msun / kpc^3
r0 = 14.7; rt = 1e3; q = 0.5; gamma = 1; beta=3
plummer                 =   [1e10, r0]
exponentialoblatehalo   =   [rho0_halo/2, r0/2, q]
ibata2024halo           =   [rho0_halo, r0, rt, q, gamma, beta]
hernquist               =   [1e10, r0]

# set the components
components = ["plummer", [1e12, r0]], ['ibata2024halo',ibata2024halo],  ['exponentialoblatehalo',exponentialoblatehalo], ['hernquist',hernquist]

### set the grid points
npoints = 110
x = np.linspace(-2*r0,2*r0,npoints); z = np.linspace(-2*r0,2*r0,npoints)
xq = np.linspace(-2*r0,2*r0,npoints//5); zq = np.linspace(-2*r0,2*r0,npoints//5)
X,Z = np.meshgrid(x,z,indexing="xy")
Xq,Zq = np.meshgrid(xq,zq,indexing="xy")
xq=Xq.flatten() ; zq=Zq.flatten() 
coordsQ=(xq,np.zeros_like(xq),zq)
coordsPHI = (X.flatten(),np.zeros_like(X.flatten()),Z.flatten())

# package the single components
singlecomponents = {}
for i in range(len(components)):
    key = "{:d}".format(i)
    singlecomponents[key] = {}
    singlecomponents[key]['component'] = components
    singlecomponents[key]["title"] = "{:.9s}".format(components[i][0])
    ax,az,PHI,comptime = get_one_comp(coordsPHI,coordsQ,components[i])
    singlecomponents[key]['ax'] = ax
    singlecomponents[key]['az'] = az
    singlecomponents[key]['PHI'] = np.reshape(PHI,coordsPHI[0].shape)

# do all permutations
permutations = {}
for i in range(len(components)):
    for j in range(len(components)):
        key = "{:d}{:d}".format(i,j)
        grav = gravitymini.gravity
        grav.cleargravity()
        grav.addgravitycomponent(*components[i])
        grav.addgravitycomponent(*components[j])
        grav.finalizegravity()
        phi = grav.evaluategravitypotential(*coordsPHI)
        ax,_,az=grav.evaluategravityforces(*coordsQ)
        ax_c,_,az_c=grav.evaluategravityforcecomponents(*coordsQ)   
        ax_cs = np.sum(ax_c,axis=0); az_cs = np.sum(az_c,axis=0)
        permutations[key] = {}
        permutations[key]['title'] = "{:.4s}/{:.4s}".format(components[i][0],components[j][0])
        permutations[key]['comp1'] = components[i]
        permutations[key]['comp2'] = components[j]
        permutations[key]['ax'] = ax 
        permutations[key]['az'] = az
        permutations[key]['phi'] = np.reshape(PHI,X.shape)
        permutations[key]["ax_c"] = ax_c
        permutations[key]["az_c"] = az_c
        permutations[key]["ax_cs"] = ax_cs
        permutations[key]["ax_cs"] = az_cs

# inspect the differences for order independence
RMS = lambda arg1, arg2 : np.sqrt(np.mean((arg1.flatten()-arg2.flatten())**2) )
# now store all of the RMS values 
for i in range(len(components)):
    for j in range(len(components)):
        key = "{:d}{:d}".format(i,j)
        key_reverse = "{:d}{:d}".format(j,i)
        axrms=RMS(permutations[key]['ax'],permutations[key_reverse]['ax'])
        azrms=RMS(permutations[key]['az'],permutations[key_reverse]['az'])
        phirms=RMS(permutations[key]['phi'],permutations[key_reverse]['phi'])
        permutations[key]["axrms"]=axrms
        permutations[key]["azrms"]=azrms
        permutations[key]["phirms"]=phirms


## plot them 
fig,axis = plt.subplots(4,4,figsize=(8.26,8),sharex=True,sharey=True)
for i in range(len(components)):
    for j in range(len(components)):
        key = "{:d}{:d}".format(i,j)
        textcolor="k"
        axis[i,j].pcolormesh(X,Z,permutations[key]['phi'],shading="gouraud")
        axis[i,j].quiver(xq,zq,permutations[key]['ax'],permutations[key]['az'])
        axis[i,j].set(aspect="equal",title=permutations[key]['title'])
        if permutations[key]["axrms"] > 0: textcolor="red"
        else: textcolor="black"
        axis[i,j].text( 0.05, 0.05, "RMS: ax {:.1e}".format(permutations[key]["axrms"]),transform=axis[i,j].transAxes,color=textcolor)
        if permutations[key]["azrms"] > 0: textcolor="red"
        else: textcolor="black"
        axis[i,j].text( 0.05, 0.15, "RMS: ax {:.1e}".format(permutations[key]["azrms"]),transform=axis[i,j].transAxes,color=textcolor)
        if permutations[key]["phirms"] > 0: textcolor="red"
        else: textcolor="black"
        axis[i,j].text( 0.05, 0.85, "RMS: phi {:.1e}".format(permutations[key]["phirms"]),transform=axis[i,j].transAxes,color=textcolor)
fig.tight_layout()
fig.savefig("combinations_composite_potential.png",dpi=300)


# three different ways for checking piece wise sum: they should all be equivalent. 
#   1. sum two independent gravity objects
#   2. composite potential: evaluategravityforces
#   3. composite potential: evaluategravityforcecomponents
# additionally, the forces from the independent gravity objects should be the same as the evaluategravityforcecomponents

for i in range(len(components)):
    for j in range(len(components)):
        keyi = "{:d}".format(i)
        keyj = "{:d}".format(j)
        key = "{:d}{:d}".format(i,j)
        key_reverse = "{:d}{:d}".format(j,i)
        print("Check if summing over evaluategravityforcecomponents is the same as evaluategravityforces")
        print("     RMS(Fnet, F1+F2)")
        print("     checking for", permutations[key]['title'])
        print("     ", RMS(permutations[key]['ax'],permutations[key]['ax_cs']))
        if not np.isclose(RMS(permutations[key]['ax'],permutations[key]['ax_cs']),0,): print("FAILED")
        else: print("SUCCESS")
        print("     checking for", permutations[key_reverse]['title'])
        print("     ", RMS(permutations[key_reverse]['ax'],permutations[key_reverse]['ax_cs']))
        if not np.isclose(RMS(permutations[key_reverse]['ax'],permutations[key_reverse]['ax_cs']),0,): print("FAILED")
        else: print("SUCCESS")
        print("")
        print("check if the component forces in gravity are the same as independently")
        print("     doing the order: ", permutations[key]["title"])
        print("    ", singlecomponents[keyi]['title'], "RMS",  RMS(singlecomponents[keyi]['ax'], permutations[key]["ax_c"][0]))
        if not np.isclose(RMS(singlecomponents[keyi]['ax'], permutations[key]["ax_c"][0]), 0): print("FAILED")
        else: print("SUCCESS")
        print("    ",singlecomponents[keyj]['title'], "RMS", RMS(singlecomponents[keyj]['ax'], permutations[key]["ax_c"][1]))
        if not np.isclose(RMS(singlecomponents[keyj]['ax'], permutations[key]["ax_c"][1]), 0): print("FAILED")
        else: print("SUCCESS")
        print("     checking reverse order: ", permutations[key_reverse]["title"])
        print("    ", singlecomponents[keyj]['title'], "RMS",  RMS(singlecomponents[keyj]['ax'], permutations[key_reverse]["ax_c"][0]))
        if not np.isclose(RMS(singlecomponents[keyj]['ax'], permutations[key_reverse]["ax_c"][0]),0):print("FAILED")
        else: print("SUCCESS")
        print("    ", singlecomponents[keyi]['title'], "RMS", RMS(singlecomponents[keyi]['ax'], permutations[key_reverse]["ax_c"][1]))
        if not np.isclose(RMS(singlecomponents[keyi]['ax'], permutations[key_reverse]["ax_c"][1]),0 ):print("FAILED")
        else: print("SUCCESS")
        print("")
        print("checking if single component sum is the same as sum of evaluategravityforcecomponents")
        axsinglesum= singlecomponents[keyi]["ax"] + singlecomponents[keyj]["ax"] 
        print("    ", permutations[key]["title"],RMS(axsinglesum,permutations[key]["ax_cs"]))
        if not np.isclose(RMS(axsinglesum,permutations[key]["ax_cs"]),0 ):print("FAILED")
        else: print("SUCCESS")

        print("    ", permutations[key_reverse]["title"],RMS(axsinglesum,permutations[key_reverse]["ax_cs"]))
        if not np.isclose(RMS(axsinglesum,permutations[key_reverse]["ax_cs"]),0 ):print("FAILED")
        else: print("SUCCESS")


        print("")
        print("checking to see if sum of single gravity module is the same as evaluategravityforces")
        print("   ", permutations[key]["title"],RMS(axsinglesum,permutations[key]["ax"]))
        if not np.isclose(RMS(axsinglesum,permutations[key]["ax"]),0 ):print("FAILED")
        else: print("SUCCESS")
        print("   ", permutations[key_reverse]["title"],RMS(axsinglesum,permutations[key_reverse]["ax"]))
        if not np.isclose(RMS(axsinglesum,permutations[key_reverse]["ax"]),0 ): print("FAILED")
        else: print("SUCCESS")