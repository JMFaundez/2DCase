import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from scipy import optimize, integrate

def shape_1(xa,ya,x_10,s1,s2):

	A = 0.2969
	B = 0.1260
	C = 0.3516
	D = 0.2843
	E = 0.1015
	thick = 0.08

	# Airfoil y-coordinate
	profile = lambda x: 5*thick*(A*np.sqrt(x) - B*x - C*x**2 + D*x**3 - E*x**4)
	# Airfoil tangent
	tan = lambda x: 5*thick*(0.5*A/np.sqrt(x) - B - 2*C*x + 3*D*x**2 - 4*E*x**3)
	# Airfoil normal
	norm = lambda x: -1/tan(x)
	ec_norm = lambda x,xp: norm(xp)*(x-xp) + profile(xp)
	# Shape of input actuator for flat plate
	b2 = lambda x,y: np.exp(-(x - x_10)**2/(s1**2) - (y)**2/(s2**2))
	# Function to find the x-coordinate of the airfoil xp
	find = lambda xp: ya - norm(xp)*(xa-xp) - profile(xp)
	tol = 1E-15
	# Actually find the intersection with the profile
	x_profile = optimize.newton(find,xa,tol=tol)
	y_profile = profile(x_profile)
	# Distance from x_p to x_10
	ds = lambda x: np.sqrt(1+(tan(x))**2)
	S,error_S = integrate.quad(ds,x_10,x_profile) 
	#print(S, error_S)
	# Normal distance to airfoil
	Dy = np.sqrt((xa-x_profile)**2 + (ya-y_profile)**2)
	b = b2(x_10+S, Dy)
	#print(b)
	return b
