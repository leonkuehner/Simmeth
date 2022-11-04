#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True

# l o a d i n i t i a l p o s i t i o n s and masses from f i l e
data = np.load('solar_system.npz')  # liest datei aus
names = data['names']
x_init = data['x_init']
v_init = data['v_init']
m = data['m']
g = data['g']
names = ['Sun', 'Earth', 'Moon','Mars','Venus','Jupiter']

dt = .0001
year = 31556952
steps = 10000
# in Astronomischen Einheiten und relativ zur Erdmasse

def force(r_ij, m_i, m_j, g):
    F = np.array([0,0])
    F = -g*m_i*m_j*r_ij/(np.sqrt(r_ij[0]**2+r_ij[1]**2)**3)
    return F
# [-g*m_i*m_j*r_ij[0]/(np.sqrt(r_ij[0]**2+r_ij[1]**2)**3),-g*m_i*m_j*r_ij[1]/(np.sqrt(r_ij[0]**2+r_ij[1]**2)**3)]
            

def step_euler(x, v, dt, mass, g, forces):
    x = x_init.copy()                   # deep copy of initial positions
    v = v_init.copy()                   # deep copy of initial velocities
    traj = np.zeros((steps,2,6))
    traj[0] = x_init
    for i in range(steps):
        x += v * dt
        v += forces(x, m,g) / m * dt
        traj[i] = x
    return traj
def forces(x, m, g):
    Force = np.zeros((2,6))
    for j in range(6):
        for i in [q for q in range(6) if q !=j]:
            rij = x[:,j]-x[:,i]
            Force[:,j] += force(rij,m[j],m[i],g)
    return Force
   # return 
#traj = np.zeros((steps,2,6))
#traj[0] = x_init
#print(traj)
traj = step_euler(x_init, v_init, dt, m, g, forces)
print(traj[:,1,2])
xerde = traj[:,0,2]
yerde = traj[:,1,2]
#for i in range(6):
#    plt.plot(traj[:,0,i]-traj[:,0,1],traj[:,1,i]-traj[:,1,1],label=names[i])
plt.plot(traj[:,0,2]-traj[:,0,1],traj[:,1,2]-traj[:,1,1],label=names[2])
#plt.plot(0,0,'o', color='royalblue')
plt.legend()
plt.xlabel(r'$x$ [au]')
plt.ylabel(r'$y$ [au]')
plt.savefig('3_2moon.pdf')
plt.show()

