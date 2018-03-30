# -*- coding: utf-8 -*-
"""
Created on Sat Apr  9 22:00:22 2016

@author: Julian
"""

import numpy as np
import pylab as py
from scipy import misc
from scipy import ndimage as im
from scipy import optimize as opt
import sys
import os

def gaussian(x, y, A, B, x0, y0, sigma_x, sigma_y, theta):
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    Gauss = A + B * np.exp( - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) + c*((y-y0)**2)))
    return Gauss


#imports the picture and does preliminary analysis
#pic = misc.imread('focused_diode_beamexpander1_50mm_4_18_2016.png')[...,0]
#pic = misc.imread('Unfocused_230_4_19_2016.png')[...,0]
#pic = misc.imread('focused_190_4_19_2016.png')[...,0]
#pic = misc.imread('GLR_unfocused_139_5_02_2016.png')[...,0]

filename = 'beamspot_glass_slide2_10_12_16.png'
name = os.path.splitext(filename)
name = name[0] + '.txt'

pic = misc.imread(filename)
if pic.ndim == 3:
    pic = pic[...,0]

#print (pic)
#sys.exit()
#pic = misc.imread('Unfocused_275_4_19_2016.png')[...,0]
#pic = misc.imread('Unfocused_177_4_20_2016.png')[...,0]
pic[0,:4] = 0 #eliminates the time stamp present in the top corner of image

h, w = pic.shape #finds the height and width of image 

picy = h//2 #turns height and width to integers
picx = w//2

pixel = 500
crop = pic[picy - pixel : picy + pixel, picx - pixel : picx + pixel] #initial crops to account for dead pixel



mp = im.measurements.maximum_position(crop * (crop != 255)) #finds the max position(intensity) of picture
#print(mp[1])
#sys.exit()
pixel2 = 50 #change according to the image
if pixel2 > mp[1]:
    pixel2 = mp[1]
 #crops the image to minimize noise
crop2 = crop[mp[0] - pixel2: mp[0] + pixel2, mp[1] - pixel2 : mp[1] + pixel2] #secondary crop



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

guess_gaussian = gaussian(x, y, *data_guess) #the values for the guess gaussian

popt, pcov = opt.leastsq(lambda p: (gaussian(x, y, *p) - data).ravel(), data_guess)

fit = gaussian(x, y, *popt)
print(popt)

py.figure()

py.subplot(121)
py.suptitle('2D Gaussian of Laser Beam Image Raw Data')
py.imshow(data, cmap =  'inferno')
py.colorbar()

py.subplot(122)
py.suptitle('2D Gaussian of Laser Beam Image with Fit')
py.imshow(data, cmap =  'hot')

GC = py.contour(x, y, guess_gaussian, 8, colors = 'pink')
GC.collections[0].set_label('Guess = pink')

FC = py.contour(x, y, fit, 8, colors = 'lime')
FC.collections[0].set_label('Fit = lime')
#py.legend()
py.show()

''' We will now convert our least square's fit values to actual physical values'''

txt_file = open(name, "wt")

print('Our least squares returned the following values in pixels for: ')
txt_file.write('Our least squares returned the following values in pixels for: \n');
txt_file.write('      A: %5.3f' % (popt[0]) + '\n');
print('      A: %5.3f' % (popt[0]))
txt_file.write('      B: %5.3f' % (popt[1]) + '\n');
print('      B: %5.3f' % (popt[1]))
txt_file.write('     x0: %5.3f' % (popt[2]) + '\n');
print('     x0: %5.3f' % (popt[2]))
txt_file.write('     y0: %5.3f' % (popt[3]) + '\n');
print('     y0: %5.3f' % (popt[3]))
txt_file.write('sigma_a: %5.3f' % (popt[4]) + '\n');
print('sigma_a: %5.3f' % (popt[4]))
txt_file.write('sigma_b: %5.3f' % (popt[5]) + '\n');
print('sigma_b: %5.3f' % (popt[5]))
txt_file.write('  theta: %5.3f' % (popt[6]) + '\n');
print('  theta: %5.3f' % (popt[6]))

#now we use our data to attain the physical waist of the beam.

#a = (np.cos(popt[6])**2)/(2*popt[4]**2) + (np.sin(popt[6])**2)/(2*popt[5]**2) #used to determine the waist in x

#c = (np.sin(popt[6])**2)/(2*popt[4]**2) + (np.cos(popt[6])**2)/(2*popt[5]**2) #used to determine the waist in y

#scale_factor = 7.8125e-7 #m per pixel
scale_factor = 7.874e-7 #m per pixel color camera

w_a = popt[4]*scale_factor*2 #2 converts from sigma to w; sigma a
w_b = popt[5]*scale_factor*2 #sigma b

r_eff = np.sqrt(w_a * w_b) #effective radius

#waist_x = (1/(np.sqrt(a)))*scale_factor
#waist_y = (1/(np.sqrt(c)))*scale_factor

print()
print('From our analysis, the beam waist for this beam is:')
txt_file.write('From our analysis, the beam waist for this beam is:\n');
# print('w_fx: ', waist_x, 'm')
# print('w_fy:', waist_y, 'm')
# print('r_eff: ', r_eff, 'm')
#print(' w_fx: %6.3f um' % (waist_x * 1E6))
#print(' w_fy: %6.3f um' % (waist_y * 1E6))
print('r_eff: %6.3f um' % (r_eff*1E6))
txt_file.write('r_eff: %6.3f um' % (r_eff*1E6) + '\n');
#the image being calculated is that after the beam has traveled through the focal lens
#so we must reverse calculate the actual waist from the focal waist, therefore for the best measurements, the beam must be focused, otherwise our calculations are not as approximate as possible

#waist_f = (waist_x + waist_y)/2
#w_f = waist_f

#parameters of optics and laser
wavelength = 532e-9
focal_length = 50e-3

w_o = (wavelength*focal_length)/(np.pi * r_eff) #calculates effective waist of the beam

Z_R_f = (np.pi * r_eff**2)/wavelength #calculates the rayleigh range after the focusing lens

Z_R = (np.pi * w_o**2)/wavelength #the distance the beam can travel before it starts to diverge

print()
print('Beam Waist w_o: %6.4f mm' % (w_o * 1E3))
txt_file.write('Beam Waist w_o: %6.4f mm' % (w_o * 1E3) +'\n');
#print('    Beam Width: %6.4f um' % (w_f * 1E6))
print()
print('Therefore, we calculate our Z_R is:')
txt_file.write('Therefore, we calculate our Z_R is:\n');
print('Z_R_f: %6.4f m or %6.3f in' % (Z_R_f, Z_R_f * 39.3701))
txt_file.write('Z_R_f: %6.4f m or %6.3f in' % (Z_R_f, Z_R_f * 39.3701)+'\n');
print('  Z_R: %6.4f m or %6.3f in' % (Z_R, Z_R * 39.3701))
txt_file.write('  Z_R: %6.4f m or %6.3f in' % (Z_R, Z_R * 39.3701)+'\n');

txt_file.close()
##pixel length of GS3_U3 in  pixels per mm is (10^(-3))/1280 or 7.8125e-7 m is 1 pixel
'''
available colors:
b is blue,
c is cyan,
g is green,
k is black,
m is magenta,
pink is pink
r is red
silver is silver
w is white
y is yellow
gold is gold
navy is navy
brown is brown
i believe that if you name the color it will just give you that color
'''