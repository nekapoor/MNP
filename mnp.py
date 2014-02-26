import math
import numpy as np
import scipy as sp
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random

def erf(x):
    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # Save the sign of x
    sign = 1
    if x.any() < 0:
        sign = -1
    x = abs(x)

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*np.exp(-x*x)

    return sign*y



def Af(r_star, t_star, pe):
	return np.exp(-1*r_star*np.sqrt(pe))*(1-erf(r_star/(np.sqrt(4*t_star)) - np.sqrt(pe*t_star)))	


def Bf(r_star, t_star, pe):
	return np.exp(r_star*np.sqrt(pe))*(1-erf(r_star/(np.sqrt(4*t_star)) + np.sqrt(pe*t_star)))


def Cf(r_star, t_star):
	return 1-erf(r_star/(np.sqrt(4*t_star)))


def exp_pet(pe, t_star):
	return np.exp(-pe*t_star)



def T_final(r_star, t_star, pe, gamma):

	AB_avg_one_minus_r = (Af(1-r_star, t_star, pe) + Bf(1-r_star, t_star, pe))/2
	AB_avg_one_plus_r = (Af(1+r_star, t_star, pe) + Bf(1+r_star, t_star, pe))/2
	C_one_minus_r = exp_pet(pe, t_star)*Cf(1-r_star, t_star)
	C_one_plus_r = exp_pet(pe, t_star)*Cf(1+r_star, t_star)
	AB_subt_one_minus_r = (Af(1-r_star, t_star, pe) - Bf(1-r_star, t_star, pe))/(2*np.sqrt(pe))
	AB_subt_one_plus_r = (Af(1+r_star, t_star, pe) - Bf(1+r_star, t_star, pe))/(2*np.sqrt(pe))

	seventh_term = exp_pet(pe, t_star)*(2*np.exp((-1*(1-r_star)**2)/(4*t_star))*np.sqrt(t_star/np.pi) - (1-r_star)*Cf(1-r_star, t_star))
	eighth_term = exp_pet(pe, t_star)*(2*np.exp((-1*(1+r_star)**2)/(4*t_star))*np.sqrt(t_star/np.pi) - (1+r_star)*Cf(1+r_star, t_star))

	initial_term = (gamma*(1-exp_pet(pe, t_star)))/(pe) - gamma/(2*r_star*pe)

	return initial_term*(AB_avg_one_minus_r - AB_avg_one_plus_r - C_one_minus_r + C_one_plus_r + AB_subt_one_minus_r - AB_subt_one_plus_r - seventh_term + eighth_term)


pp = 1
g = 1


#rs = np.arange(0,2,0.1).tolist()
#ts = np.arange(0.1,50,1).tolist()
#X,Y = np.meshgrid(rs, ts) # grid of point

x = linspace(0.001, 2, 10)
Z = T_final(x, 2, pp, g) # evaluation of the function 
figure()
plot(x, Z, 'r')
xlabel('r*')
ylabel('T*')
title('title')

# fig = plt.figure(figsize=(14,6))

# # `ax` is a 3D-aware axis instance because of the projection='3d' keyword argument to add_subplot

# # surface_plot with color grading and color bar
# ax = fig.add_subplot(1, 2, 2, projection='3d')
# p = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
# cb = fig.colorbar(p, shrink=0.5)
show()