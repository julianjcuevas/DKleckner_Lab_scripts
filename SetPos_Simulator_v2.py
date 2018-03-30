import time
import os, sys
from hoomd import *
from hoomd import md
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import hoomd
from scipy import signal as si
import math


context.initialize('--mode=cpu')

Num_colloids = 100
spawn_radius = 75

cuttoff_radius = 30
Amp = 10
freq = 0.5
phi = 0

eps = 7
sig = 2

temperature = 0.3
delta_t = 0.005

count = 0
max_run = 1e6;
run_sum = 0;
ran_run_max = 5e4;
###########################################################################################

hour = time.strftime("%H")
minute = time.strftime("%M")
day = time.strftime("%d")
month = time.strftime("%m")
year = time.strftime("%Y")
date = month + "/" + day + "/" + year
label = month + '_' + day + '_' + year + '_' + hour + '_' + minute

###########################################################################################
xy = '/home/labuser/hoomd-blue-bind2/build-phase/Colloid_Coordinates' + 'Num_of_Colloids_' + str(Num_colloids) + '_Initial_Radius_' + str(spawn_radius) + '.txt'

output_dir = os.path.expanduser('Simulation_' + label)
if not os.path.exists(output_dir): os.makedirs(output_dir)

path = "/home/labuser/hoomd-blue-bind2/build-phase/" + "Simulation_" + label

if os.path.exists(output_dir):
    os.chdir(path)
###########################################################################################

def ob_potential(A, k, r, phi):
    theta = k * r
    V_r = A*((np.pi/2) - (np.cos((theta + phi))/theta) - (np.sin((theta + phi))/(theta**2)))
    return V_r
'''    
def ob_force(A, k, r, phi):
    theta = k * r
    F_r = -A*np.sin((theta + phi))*((2/theta**3)-(1/theta))
    return F_r
'''

def temperature_well(V_r, r, temperature):
    rel_max = si.argrelextrema(V_r, np.greater)[0] #argrelextrema extracts the indices in form of a tuple
    rel_min = si.argrelextrema(V_r, np.less)[0]
    rel_diff = V_r[rel_max[1]] - V_r[rel_min[1]]
    temp_vs_ob = temperature/V_r[rel_max[1]]
    const_temp_line = V_r[rel_min[1]] + (rel_diff * temp_vs_ob)
    print(V_r[rel_min[1]])
    return const_temp_line
###########################################################################################


r = np.linspace(1, cuttoff_radius, num = (cuttoff_radius*10))
Vr = ob_potential(Amp, freq, r, phi)
T = temperature_well(Vr, r, temperature)
T = np.repeat(T, len(r))

plt.plot(r, Vr, color='m')
plt.plot(r, T, color = 'cyan')
plt.xlabel('Radial distance, r', fontsize = 'large')
plt.ylabel('V(r)', fontsize = 'large')
plt.suptitle('Optical Binding Potential', fontsize = 'large')
plt.savefig('Pplot.png', dpi = 192)

###########################################################################################

conditions = open('Conditions.txt', 'wt')

conditions.write(str(Amp) + '\n')
conditions.write(str(freq) + '\n')
conditions.write(str(phi) + '\n')
conditions.write(str(temperature) + '\n')
conditions.write(str(delta_t)+'\n')
conditions.write(str(cuttoff_radius) + '\n')
conditions.write(str(ran_run_max) + '\n')
conditions.write(str(eps) + '\n')
conditions.write(str(sig) + '\n')
conditions.write(str(spawn_radius))

conditions.close()

###########################################################################################

initial = open('Simulation_Initial_Conditions.txt', 'wt')

initial.write("The initial conditions for this simulation are: \n")
initial.write("Optical Binding Force: \n")
initial.write("Amplitude: " + str(Amp) + "\n")
initial.write("Frequency: " + str(freq) + "\n")
initial.write("Phase: " + str(phi) + "\n")
initial.write("Simulation Conditions: \n")
initial.write("Temperature: " + str(temperature) + "\n")
initial.write("Delta t: " + str(delta_t) + "\n")
initial.write("Cutoff Radius: " + str(cuttoff_radius) + "\n")
initial.write("Simulation Duration: \n")
initial.write("Brownian Motion, no optical binding: " + str(ran_run_max) + "time steps \n")
initial.write('WCA Epsilon: ' + str(eps) + '\n')
initial.write('WCA Sigma: ' + str(sig) + '\n')
initial.write("Total time: " + str(max_run) + "\n")
initial.write("Time Stamp: " + str(label))

initial.close()
###########################################################################################
'''
There are several initial lattices that can be initialized with HOOMD.
The square lattice, sq
The hexagonal lattice, hex
The cubic lattice, sc
The body centered cubic, bcc
The face centered cubic, fcc

We're sticking to the 2-D lattice structures. But, you are also able to create your own unique
unit cell and control other parameters.

'''
file = open(xy, 'r')

coordinates = np.loadtxt(file, dtype = float, delimiter = '\t')

snap = data.make_snapshot(N = Num_colloids, box = hoomd.data.boxdim(L = (20*spawn_radius), dimensions = 2))

for i in range(0, len(coordinates)):
    snap.particles.position[i] = coordinates[i, :]


system = init.read_snapshot(snap)

nl = md.nlist.cell()

wca = md.pair.lj(r_cut=2.5, nlist=nl)
wca.set_params(mode="shift")
wca.pair_coeff.set('A', 'A', epsilon=eps, sigma=sig, r_cut=((2**(1./6.))*2))

all = group.all();
md.integrate.mode_standard(dt=delta_t)
md.integrate.brownian(group=all, kT=temperature, seed=987)

###########################################################################################

data_dir = os.path.expanduser('Simulation_Text_Files')
if not os.path.exists(data_dir): os.makedirs(data_dir)

path2 = "/home/labuser/hoomd-blue-bind2/build-phase/" + "Simulation_"  + label + "/Simulation_Text_Files"

if os.path.exists(data_dir):
    os.chdir(path2)

###########################################################################################

txt_file = open('coordinates%06d.txt'%count, 'wt')

run(0)
snapshot =system.take_snapshot(all=True)
pos = snapshot.particles.position
x_pos = pos[...,0].ravel()
y_pos = pos[...,1].ravel()
z_pos = pos[...,2].ravel()
    
x = x_pos.astype(str)
y = y_pos.astype(str)
z = z_pos.astype(str)
    
for i in range(len(x)):
    txt_file.write(x[i] + "\t" + y[i] + "\t" + z[i] + "\n")
    
txt_file.close()

############################################################################################

while run_sum <= ran_run_max:
    txt_file = open('coordinates%06d.txt'%count, 'wt')
    run(1e3)
    run_sum += 1e3
    count += 1
    snapshot =system.take_snapshot(all=True)
    pos = snapshot.particles.position
    x_pos = pos[...,0].ravel()
    y_pos = pos[...,1].ravel()
    z_pos = pos[...,2].ravel()
    
    x = x_pos.astype(str)
    y = y_pos.astype(str)
    z = z_pos.astype(str)
    
    for i in range(len(x)):
        txt_file.write(x[i] + "\t" + y[i] + "\t" + z[i] + "\n")
    
    txt_file.close()

##############################################################################################

nl = md.nlist.cell()
ob = md.pair.OBind(r_cut=cuttoff_radius, nlist=nl)
ob.pair_coeff.set('A', 'A', amplitude= Amp, frequency=freq, phase = phi)
ob2 = md.pair.OBind(r_cut = cuttoff_radius, nlist=nl)
ob2.pair_coeff.set('A', 'A', amplitude = Amp, frequency = (freq/math.sqrt(2)), phase = (phi + math.pi/4))

##############################################################################################

while run_sum <= max_run:
    #txt_file = open('coordinates' + str(count) + '.txt', 'wt')
    txt_file = open('coordinates%06d.txt'%count, 'wt')
    run(1e3)
    run_sum += 1e3
    count += 1
    snapshot =system.take_snapshot(all=True)
    pos = snapshot.particles.position
    x_pos = pos[...,0].ravel()
    y_pos = pos[...,1].ravel()
    z_pos = pos[...,2].ravel()
    
    x = x_pos.astype(str)
    y = y_pos.astype(str)
    z = z_pos.astype(str)
    
    for i in range(len(x)):
        txt_file.write(x[i] + "\t" + y[i] + "\t" + z[i] + "\n")
    
    txt_file.close()