# 
# HG.py
#
# Description:
#	python snippet that produces a Henyey-Greenstein 
#	phase matrix.
#
# Author: Patrick Stegmann
# Date: 2018-10-22
#       2020-12-27
#
import numpy as np


class HG:
	
	def __init__(self,asym):
		self.g = asym

	def Delta(self):
  		return (0.5*(10./9.*self.g + 250./729.*self.g**3))**2 \
			+ (1./3.*(4. - 25./27.*self.g**2.))**2

	def BigG(self):
  		return 5./9.*self.g \
			+ (0.5*(10./9.*self.g + 250./729.*self.g**3) \
					+ np.sqrt(self.Delta()))**(1./3.) 
				- (np.sqrt(self.Delta()) - 0.5*(10./9.*self.g \
							+ 250./729.*self.g**3.))**(1./3.)

	def phasefu(self,theta):
		Prefactor = 2./(2.+self.BigG()**2.) \
					* (1. - self.BigG()**2.)/(1.+self.BigG()**2. \
							- 2.*self.BigG()*np.cos(theta))**(3./2.)
		P11 = Prefactor*(1. + np.cos(theta)**2.)
		P12 = Prefactor*(-1. + np.cos(theta)**2.)
		P33 = Prefactor*np.cos(theta)
		return P11
