import math
import numpy as np
import scipy as sp
from pylab import *
from scipy.special import erfc
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random

pp = 1
g = 1





def Af(r_star, t_star, pe):
	return np.exp(-1*r_star*np.sqrt(pe))*1*(erfc(r_star/(np.sqrt(4*t_star)).any() - np.sqrt(pe*t_star)))	


def Bf(r_star, t_star, pe):
	return np.exp(r_star*np.sqrt(pe))*1*(erfc(r_star/(np.sqrt(4*t_star)) + np.sqrt(pe*t_star)))


def Cf(r_star, t_star):
	return 1*(erfc(r_star/(np.sqrt(4*t_star))))


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

	initial_term = (gamma*(1-exp_pet(pe, t_star)))/(pe) - gamma/(2*r_star)*pe

	return initial_term*(AB_avg_one_minus_r - AB_avg_one_plus_r - C_one_minus_r + C_one_plus_r + AB_subt_one_minus_r - AB_subt_one_plus_r - seventh_term + eighth_term)


def T_steady_less_than_1(r_star, t_star, pe, gamma):
    first = gamma/pe
    second_first = (1+1/np.sqrt(pe))
    second_second = gamma*np.exp(-np.sqrt(pe))
    second_third = (np.exp(-r_star*np.sqrt(pe)) - np.exp(-r_star*np.sqrt(pe)))
    second_fourth = (1/(2*r_star*pe))
    second = second_first*second_second*second_third*second_fourth

    return first+second

def T_final_greater_than_1(r_star, t_star, pe, gamma):
    AB_avg_r_minus_one = (Af(r_star-1, t_star, pe) + Bf(r_star-1, t_star, pe))/2
    AB_avg_r_plus_one = (Af(r_star+1, t_star, pe) + Bf(r_star+1, t_star, pe))/2
    C_one_minus_r = exp_pet(pe, t_star)*Cf(r_star-1, t_star)
    C_one_plus_r = exp_pet(pe, t_star)*Cf(r_star+1, t_star)
    AB_subt_r_minus_one = (Af(r_star-1, t_star, pe) - Bf(r_star-1, t_star, pe))/(2*np.sqrt(pe))
    AB_subt_r_plus_one = (Af(r_star+1, t_star, pe) - Bf(r_star+1, t_star, pe))/(2*np.sqrt(pe))

    seventh_term = exp_pet(pe, t_star)*(2*np.exp((-1*(r_star-1)**2)/(4*t_star))*np.sqrt(t_star/np.pi) - (r_star-1)*Cf(r_star-1, t_star))
    eighth_term = exp_pet(pe, t_star)*(2*np.exp((-1*(r_star+1)**2)/(4*t_star))*np.sqrt(t_star/np.pi) - (r_star+1)*Cf(r_star+1, t_star))

    initial_term = gamma/((2*r_star)*pe)

    return initial_term*(AB_avg_r_minus_one + AB_avg_r_plus_one - C_one_minus_r - C_one_plus_r - AB_subt_r_minus_one + AB_subt_r_plus_one + seventh_term - eighth_term)

def T_steady_greater_than_1(r_star, t_star, pe, gamma):
    first = (1/(2*r_star*pe))
    second = np.exp(-r_star*np.sqrt(pe))
    exp_pos = np.exp(np.sqrt(pe))
    exp_neg = np.exp(-np.sqrt(pe))
    inv_sqrt_pe = 1/np.sqrt(pe)
    third = exp_pos + exp_neg - inv_sqrt_pe*(exp_pos - exp_neg)

    return first*second*third
  

#rs = np.arange(0,2,0.1).tolist()
#ts = np.arange(0.1,50,1).tolist()
#X,Y = np.meshgrid(rs, ts) # grid of point

x = linspace(0.001, 3, 6)
# x1 = linspace(2, 3, 2)
Z = T_final(x, 2, pp, g) # evaluation of the function 
Z_1 = T_steady_less_than_1(x, 2, pp, g)
Z_2 = T_final_greater_than_1(x, 2, pp, g)
Z_3 = T_steady_greater_than_1(x, 2, pp, g)
# Z_1 = T_final_greater_than_1(x1, 2, pp, g)
figure()
plot(x, Z+Z_1+Z_2+Z_3, 'r')
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