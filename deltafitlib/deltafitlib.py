'''

Description:
============
A simple Python code for applying the delta-fit
method of Hu et al., 2002 to a given input
phase function.
The computed output are the delta-fit
expansion coefficients and the fitting error.
The results are illustrated as plots.

Author: Patrick Stegmann 
Date:   2019-2020

'''

import numpy as np
from scipy.special import legendre

class DeltaFit:
	
	''' Specify cut-off angle '''
	cutoff = 3.0 # [degrees]

	def weight(self,angle):
		if angle < cutoff:
			return 0.
		else:
			return 1.

	
