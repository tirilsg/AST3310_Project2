import numpy as np
import matplotlib.pyplot as plt

def cross_section_mod(R, L, F_C, ax=None, sanity=False, savefig=False, name="name"):
    R_sun      = 6.96e8         # m
    r_range    = 1.2#*np.max(R)  # normalized 
    core_limit = 0.995 * np.max(L)
    
    star_zone = np.int32(np.where(L>core_limit,0.5,-0.5) * np.where(F_C>0,3,1) + 2)
    colors=['b','c','yellow','r']
    if ax is None:
        fig, ax = plt.subplots(figsize=(6,6))
    else:
        fig = ax.figure

    ax.set_xlim(-r_range, r_range)
    ax.set_ylim(-r_range, r_range)
    ax.set_aspect('equal')
    star_zone_prev = -1
    for k in range(0,len(R)-1):
        if star_zone[k]!=star_zone_prev: # only plot a new *filled* circle if a new star_zone
            star_zone_prev = star_zone[k]
            circle = plt.Circle((0,0),R[k]/R[0],fc=colors[star_zone[k]],fill=True,ec=None)
            ax.add_artist(circle)
    circle_white=plt.Circle((0,0),R[-1]/R_sun,fc='white',fill=True,lw=0)
    ax.add_artist(circle_white)
    circle_red    = plt.Circle((2*r_range, 2*r_range), 0.1*r_range, color=colors[3], fill=True)
    circle_yellow = plt.Circle((2*r_range, 2*r_range), 0.1*r_range, color=colors[2], fill=True)
    circle_cyan   = plt.Circle((2*r_range, 2*r_range), 0.1*r_range, color=colors[1], fill=True)
    circle_blue   = plt.Circle((2*r_range, 2*r_range), 0.1*r_range, color=colors[0], fill=True)
    ax.legend([circle_red,circle_yellow,circle_cyan,circle_blue],\
            [f'Outer Convection Zone ',f'Outer Radiation Zone',\
                f'Inner Radiation Zone', f'Inner Convection Zone '], loc="lower right",fontsize=14)
    return star_zone
