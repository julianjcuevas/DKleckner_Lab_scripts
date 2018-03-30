# DKleckner_Lab_scripts
This is a collection of scripts that I used during my undergraduate research at UC Merced

My undergraduate research involved the analysis of a laser beam spot image after incidence in a focusing lens to determine
the physical parameters of the beam. Continuous Wave Laser's generally carry a gaussian profile, so I used regression analysis
to fit a gaussian onto the image and determine the waist and astigmatism of the beam. The waist gives us insight into the
laser's long range behavior. 

Aside from fitting the gaussian onto the beam, I also conducted particle simulations, where I had to modify files from HOOMD
and implement our own inter-particle force into the software(I no longer have access to these files). From this, I wrote a
few scripts that ran the simulation and stored the data and initial parameters for the simulation in a streamlined fashion. 
The data could later be analyzed and converted into a video that overlayed our particle force potential in order to determine 
the interactions.
