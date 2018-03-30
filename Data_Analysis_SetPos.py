import os
import numpy as np
import pylab as P
from ob_forces import *

rootdir = os.getcwd()
path = 'Simulation_Text_Files'
fn_prefix = 'coordinates'

output_dir = os.path.expanduser('Data')
if not os.path.exists(output_dir): os.makedirs(output_dir)



texts = rootdir + '/' + path


############################################################################################

conditions = np.loadtxt('Conditions.txt', dtype = float, delimiter = '\t')

Amp = conditions[0]
freq = conditions[1]
phase = conditions[2]
temperature = conditions[3]
delta_t = conditions[4]
r_cut = conditions[5]
random_run = conditions[6]
eps = conditions[7]
sig = conditions[8]

############################################################################################
count = 0
time_step = 1e3
fig = P.figure(figsize=(9, 5), dpi=192)

force = OB()
'''

for root, dirs, files in os.walk(path):
    for file in files:
        print(file)
        #if not file.startswith(fn_prefix): continue
        
        #name = fn_prefix + str(count) + '.txt'
        f = open(file, 'r')
        positions = np.loadtxt(f, dtype=float, delimiter = '\t')
        
        X = positions[:,0]
        Y = positions[:,1]
        
        ax = fig.add_subplot(111)
        
        circs = []
        
        
        for n in range(0, len(X)):
            circs.append(P.Circle((X[n], Y[n]), np.pi, facecolor='blue', edgecolor='k', lw=1))
            
        
        for circ in circs:
            ax.add_artist(circ)
        
        P.savefig(output_dir, dpi=192)
        
        count += 1
'''
if os.path.exists(path):
    os.chdir(path)

for filename in sorted(os.listdir(texts)):
    if not filename.startswith(fn_prefix): continue
    
    
    f = open('coordinates%06d.txt'%count, 'r')
    X = np.loadtxt(f, dtype=float, delimiter = '\t')
        
    t = time_step*delta_t*count
    
    if (time_step*count) <= random_run:
        V_on = False
    else:
        V_on = True
    
        #fig = P.figure(figsize=(8, 4), dpi=192)
        
    
        
    ax1 = fig.add_axes([.05, .1, .4, .8])
    ax2 = fig.add_axes([.55, .1, .4, .8])
        
    circs = []
    circs2 = []
        
    for n in range(len(X)):
            # P.plot(Xt[:, n, 0], Xt[:, n, 1], 'w')
        circs.append(P.Circle((X[n, 0], X[n, 1]), np.pi, facecolor='blue', edgecolor='k', lw=1))
        circs2.append(P.Circle((X[n, 0], X[n, 1]), np.pi, facecolor='none', edgecolor='k', lw=0.5))
        
        # ax = P.gca()
        
        
    xr = 70
    yr = xr
    nr = 500
        
    y, x = np.ogrid[yr:-yr:nr*1j, -xr:xr:nr*1j]
        
    step = len(X)
        
    if V_on:
        V = 0
            # for i in range(count2*10, (count2+1)*10):
       # R = freq*X + phase
        for i in range(0, len(X), step):
            i0 = i
            V += force.compute_2D_grid_V(X[i0:i0+step, :2], x, y)
            V = Amp*V   
        print(V)
    else:
        V = 0*x*y
    
    
    #img = ax2.imshow(V, extent=[-xr, xr, -yr, yr], cmap='inferno', clim=(.75, (-0.75)))
    img = ax2.imshow(V, extent=[-xr, xr, -yr, yr], cmap='inferno')
    
    
    for circ in circs:
        ax1.add_artist(circ)
    
    for circ in circs2:
        ax2.add_artist(circ)
        
    for ax in (ax1, ax2):
        ax.set_aspect(1)
        
    ax1.set_xlim(-225, 225)
    ax1.set_ylim(-225, 225)
        
    ax2.set_xlim(-xr, xr)
    ax2.set_ylim(-yr, yr)
        
    ax1.set_xlabel('x-position')
    ax1.set_ylabel('y-position')
        
    P.xlabel('x-position')
    P.ylabel('y-position')
        #ax1.ylabel('y')
        
    ax1.add_artist(P.Rectangle((-xr, -yr), xr*2, yr*2, facecolor='none', edgecolor='red', lw=1))
        
    fig.colorbar(img)
        #fig.clf()
    P.title('t=%6.2f' % t, family='monospace')
        #P.savefig(ofn, dpi=192)
        #fig.clf()
    
        #P.close(fig)
   
    
    data = os.path.join(output_dir, '%08d.png'%count)
    print('t=%7.3f => %s' % (t, data))
    if os.getcwd() is not rootdir:
        os.chdir(rootdir)
        P.savefig(data, dpi=192)
        fig.clf()
        os.chdir(texts)
        
    count += 1