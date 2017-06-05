# coding: utf-8
import scipy.special as ss
import scipy.interpolate
import numpy as np
import matplotlib.pyplot as plt

def ogata_hankel(lnum):
	#lnum = 1e4
	h = 1e-5
	tanh, sinh, pi = np.tanh, np.sinh, np.pi
	Jnu = ss.jv
	xi = ss.jn_zeros
	w_nul = 2 / (pi**2 * xi(2,lnum) * Jnu(3, pi*xi(2,lnum)))
	psi = lambda t: t*tanh((pi/2)*sinh(t))
	def fps(func,x,h):
	    return (-func(x+2*h)+8*func(x+h)-8*func(x-h)+func(x-2*h))/(12*h)
	def P(interp_points):
		func_x, func_y = k, k*pk
		return np.interp(interp_points, func_x, func_y)
	
	X = (pi/h) * psi(h*xi(2, lnum))
	J2X = Jnu(2, X)
	p2 = w_nul * J2X * fps(psi, h*xi(2,lnum), h)

	p1 = k * pk
	Sigma = p1 * p2
	Sigma = sum(Sigma)
	Sigma *= pi
	return Sigma
