import numpy as np
import pylab as py
from scipy import misc
from scipy import ndimage as im
from scipy import optimize as opt
import os

''' The following script is used to take multiple beam spot images -- that were saved under similar conventions -- and 
used to extract the laser's parameters after incidence on a focusing lens.

'''

framedir = '.'
fn_prefix = "GLR_unfocused_"

zs = []
widths = []

#parameters of optics and laser
wavelength = 532e-9
focal_length = 50e-3

scale_factor = 7.8125e-7 #m per pixel
min_reff = 10000000 #outrageous number

def gaussian(x, y, A, B, x0, y0, sigma_x, sigma_y, theta):
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    Gauss = A + B * np.exp( - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) + c*((y-y0)**2)))
    return Gauss

def crop(pic):
    pic[0,:4] = 0 #eliminates the time stamp present in the top corner of image

    h, w = pic.shape #finds the height and width of image 

    picy = h//2 #turns height and width to integers
    picx = w//2

    pixel = 200
    crop = pic[picy - pixel : picy + pixel, picx - pixel : picx + pixel] #initial crops to account for dead pixel



    mp = im.measurements.maximum_position(crop * (crop != 255)) #finds the max position(intensity) of picture

    pixel2 = 50 #change according to the image
     #crops the image to minimize noise
    crop2 = crop[mp[0] - pixel2: mp[0] + pixel2, mp[1] - pixel2 : mp[1] + pixel2] #secondary crop
    return crop2

for filename in os.listdir(framedir):
#for filename in filter(lambda f: f.startswith(fn_prefix), os.listdir(framedir)):
    if not filename.startswith(fn_prefix): continue
    
    num = int(filename[len(fn_prefix):][:3])
    filename = os.path.join(framedir, filename)
    
    pic = misc.imread(filename)[...,0]
    crop2 = crop(pic)
    mp2 = im.measurements.maximum_position(crop2 * (crop2!= 255)) #finds the new max position after the secondary crop
    #establishes the numerical length of x and y by finding the length of a sliver of the image
    
    N = len(crop2[mp2[0], :]) #we must use the secondary crop of the image to get the right N value
    x = np.arange(N)
    y = np.arange(N)

    #creates the plot grid
    x, y = np.meshgrid(x, y)
        
    '''The raw data plots will now be constructed'''
    
    data = crop2

    '''Now we will attempt to establish the guess parameters for the least square fitting'''
    
    com = im.measurements.center_of_mass(crop2)
    x0 = com[1]
    y0 = com[0]
    
    #establishes the guess tuple for our distribution
    data_guess = (
        data.min(), #finds the minimum values of the image/array
        data.max() - data.min(), #the guess for the amplitude is the max value minus the min value(we are accounting for noise)
        x0, # I am assuming that there may be some stray light or brightspots that may vary the center of mass
        y0,
        N//5, #i am taking a standard deviation along an axis to as a guess for sigma x and y
        N//5,
        0 #since the raw data appears circular I am taking a theta of 0
    )
    
    popt, pcov = opt.leastsq(lambda p: (gaussian(x, y, *p) - data).ravel(), data_guess)
    
    a = (np.cos(popt[6])**2)/(2*popt[4]**2) + (np.sin(popt[6])**2)/(2*popt[5]**2) #used to determine the waist in x
    c = (np.sin(popt[6])**2)/(2*popt[4]**2) + (np.cos(popt[6])**2)/(2*popt[5]**2) #used to determine the waist in y
    
    w_a = popt[4]*scale_factor*2 #2 converts from sigma to w; sigma a
    w_b = popt[5]*scale_factor*2 #sigma b
    
    r_eff = np.sqrt(w_a * w_b) #effective radius
    
    
    widths.append(r_eff)
    zs.append(num/1000 * 25.4E-3)
    
    # if r_eff < min_reff:
    #     min_reff = r_eff
    #     z0 = zs[-1]
     
    
    print('Position:', (num/1000), 'in')
    print('Radius:', (r_eff*1E6),'um')
 
 
zs = np.array(zs) 
widths = np.array(widths)

i_min = np.argmin(widths)
min_reff = widths[i_min]
z0 = zs[i_min]
Z_R_f = (np.pi * min_reff**2)/wavelength #calculates the rayleigh range after the focusing lens

zs -= z0

def theoreticalwidth(z, wo, zr):
    return wo*np.sqrt(1+(z/(zr))**2)

z = np.linspace(zs[0], zs[len(zs)-1])

fig, ax1 = py.subplots()
theory = theoreticalwidth(z, min_reff, Z_R_f)
ax1.plot(z, theory, 'b-')
py.xlabel('Position Z (mm)')
py.ylabel('w(z) (um)')
py.suptitle('Beam Width Vs. Position (theoretical)')

ax2 = ax1.twiny()
ax2.plot(zs, widths, 'go')

py.show()
