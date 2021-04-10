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

class DeltaFit(object):
	
	def __init__(self, \
			cutoff \
			theta \
			pacc \
			NStreams):
		''' Specify cut-off angle '''
		self.cutoff = cutoff                # Cut-off angle [degrees]
		self.theta = theta                  # Scattering angle [degrees]
		self.pacc = pacc                    # Phase function values [-]
		self.nx = 2.*theta/theta[-1] - 1.
		self.NStreams = NStreams            # Number of streams [-]
		self.rror = np.zeros(len(NStreams)) # Fitting error [-]


	def weight(self,angle):
		if angle < self.cutoff:
			return 0.
		else:
			return 1.

	def fit(self):
		from scipy.special import legendre
		tt = 0
		for kk in self.NStreams:
			legvals = np.zeros((self.theta.size,kk))
			A = np.zeros((self.theta.size,kk))
			
			''' Compute the matrix '''
			for jj in range(kk):
				Pn = legendre(jj)
				for ii in range(self.theta.size):
					legvals[ii,jj] = Pn(self.nx[ii])
			
			for ii in range(self.theta.size):
				for jj in range(kk):
					A[ii,jj] = legvals[ii,jj] \
							   /self.pacc[ii] \
					*self.weight(self.theta[ii])
			''' Compute the least-squares fit
			  of the expansion coefficients '''
			b = np.ones(self.theta.size)
			''' U, S, V^T = np.linalg.svd(a) '''
			sol = np.linalg.lstsq(A,b,rcond=None)
			coef = sol[0]

			''' Computing delta-fit phase function '''
			dfit = np.polynomial.legendre.legval(self.nx,coef)
			for ii in range(10, self.theta.size):
				self.rror[tt] += (dfit[ii] - self.pacc[ii])**2
			tt += 1
			return dfit

	def plotPhaseFunction(self,theta,phfu,kk):
		import matplotlib.pyplot as plt
		plt.figure()
		plt.semilogy(theta,phfu,label=str(kk))
		plt.xlabel(r'Scattering Angle $\theta$ [deg]')
		plt.ylabel(r'Phase Function $P_{11}$ [-]')
		plt.show()



