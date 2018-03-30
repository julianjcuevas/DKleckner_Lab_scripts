import random
import os
import numpy as np

#############################################################################
'''
Main body of script referenced from stack overflow
'''

Num_colloids = 20 #set the number of particles you wish to simulate
spawn_radius = 50 # the maximum radius in which particles will be distributed
Colloid_min_distance = 3 # the minimum distance between particles when particles are generated
range_X = (-spawn_radius, spawn_radius)
range_Y = (-spawn_radius, spawn_radius)

label = 'Num_of_Colloids_' + str(Num_colloids) + '_Initial_Radius_' + str(spawn_radius) #sets the label to distinguish between files.

deltas = set() #creates a python object that will allow for detection and removal of duplicate postions, as well as allow
               # for setting the exclusions for position generation

for x in range(-Colloid_min_distance, Colloid_min_distance+1):
    for y in range(-Colloid_min_distance, Colloid_min_distance+1):
        if x**2 + y**2 <= Colloid_min_distance**2:
            deltas.add((x,y))
            
random_points = []
excluded_points = set()

i = 0

while i < Num_colloids:
    x = random.randrange(*range_X)
    y = random.randrange(*range_Y)
    if (x,y) in excluded_points: continue
    random_points.append((x,y))
    i += 1
    excluded_points.update((x+dx, y+dy) for (dx,dy) in deltas)
    
coord = np.asarray(random_points)
zeros = np.zeros((Num_colloids,1))
coordinates = np.append(coord,zeros, axis = 1)

X = coordinates[:,0]
Y = coordinates[:,1]
Z = coordinates[:,2]

txt = open('Colloid_Coordinates' + label + '.txt', 'wt')

for i in range(len(X)):
    txt.write(str(X[i]) + "\t" + str(Y[i]) + "\t" + str(Z[i]) + "\n")

txt.close()