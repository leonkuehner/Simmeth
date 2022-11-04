#!/usr/bin/env python3

import numpy as np

import ex_2_1

"""""
First of all we define all constants and functions which we're about to use
"""""
# CONSTANTS=====================================

m = 2.0  # Mass in [kg]

g = 9.81  # Gravitational acceleration [m/s^2]

Gammas = np.array([0, 0.1])  # Friction coefficients [kg/s]

V_W = np.array([[0, 0], [-50, 0]])  # Wind velocity fields [m/s]
vw = np.array([0., 0.])
Trajectories = []
x = np.array([0.0, 0.0])
v = np.array([50.0, 50.0])
dt = 0.1

"""""
Secondly we define all functions which we're going to use. 
-The force field function 
-The run function (this accepts a certain set of parameters and calculates a corresponding trajectory  
-The main function 
"""""


def force(mass, gravity, v, gamma,
          v_0):  # v_0 is the wind field vector and v is the current speed vector of the particle
    F_g = np.array([0, -mass * gravity])
    F_fric = -gamma * (v - v_0)
    F_tot = F_g + F_fric
    return F_tot


def run(x0, v0, dt, mass, gravity, gamma, v_0):
    #trajectory = [x]
    x = x0.copy()       # deep copy of initial position
    v = v0.copy()       # deep copy of initial velocity
    traj = [x.copy()]   # list of positions (trajectory)
    while x[1] >= 0:    # check if x coordinate is below ground
        x += v * dt
        v += force(mass, gravity, v, gamma,v_0) / m * dt
        traj.append(x.copy())
    return traj

"""""
Inside the Main function we examine the effects of friction and wind speed.
In the second part we look for the wind speed for which the canon ball is hitting the ground at the initial position
"""""
def main():
#============ PART 3.2a
    for gamma,vww in [(0., 0.),(.1, 0.),(.1,-50.)]:
        vw = np.array([vww, 0.])
        ctraj = np.array(run(x, v, dt, m, g, gamma, vw))
        plt.plot(ctraj[:, 0], ctraj[:, 1], '-',label=f'$\gamma= {gamma}$, $v_w= {vww}$ m/s')
#label=r'$\gamma = \num{{ {} }}$; $v_w = \SI{{ {} }}{{\m\per\s}}$'.format(gamma, int(vw[0]))
 #   plt.figure()
    plt.legend(loc='upper right')
 #   plt.grid()
    plt.axis([0,550,0,150])
    plt.title('Trajectory of a cannonball')
    plt.xlabel('distance x [m]')
    plt.ylabel('height h [m]')
 #   plt.show()
    plt.savefig("./cannonball_examples.png", bboxes_inches="tight")
 #   plt.close()
#============ PART 3.2b
    amp = .1  # Amplifier
# We calculate the hitting position of the canon ball eg the first component of last array
    hit = run(x, v, dt, m, g, gamma, vw)[-1][0]
# As long as the hitting x-position is not greater the 10^-6 m, the wind speed will be increased
    while abs(hit) >= 1e-6:
        hit = run(x, v, dt, m, g, gamma, vw)[-1][0]
        vw[0] -= amp * hit
    print(vw)

    # Trajectories for different wind speeds
    plt.figure()
    plt.xlabel(r'distance [m]')
    plt.ylabel(r'height [m]')

    for vww in np.linspace(0, -200, 9):
        vw = np.array([vww, 0.])
        ctraj = np.array(run(x, v, dt, m, g, gamma, vw))
        plt.plot(ctraj[:, 0], ctraj[:, 1], '-',label = f"$v_w $= {vww} m/s ")
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.axis([-10,480,0,120])
    plt.savefig("./cannonball_vary.png", bboxes_inches="tight")
#    plt.savefig("./Gcannonball_vary.pgf", bboxes_inches="tight")
    plt.close()


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    main()
