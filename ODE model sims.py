# -*- coding: utf-8 -*-
"""
Created on Mon Feb  3 01:12:23 2025

@author: theoj
"""

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import numpy as np
import scipy.integrate as integrate
from scipy.integrate import solve_ivp

#First-order ODE for temperature T(t)
#   T'(t) = c*exp(bT) - hT

#Function for T'(t)
def temp(t, T, c, b, h):
    T_prime = (c*np.exp(b*T)) - (h*T)
    return T_prime


#Solution for T(t) given values for constants c,b,h and initial condition T(0) = T0
def Temp_Sol(c,b,h,T0):
    sol = solve_ivp(fun=temp, y0=T0, t_span=[0,40], t_eval=np.linspace(0,40,500), args=(c,b,h,), dense_output=True)
    return [sol.t, sol.y]


#Plot temperature over time
def Temperature_Plot(c,b,h,T0):
    fig=plt.figure()
    plt.plot(Temp_Sol(c,b,h,T0)[0], Temp_Sol(c,b,h,T0)[1][0])
    #plt.title(f'λ = {int(round(h/(c*b), 0))}, T(0) = {T0[0]}')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature')
    plt.grid()
    return fig


#Define lambda = h/bc
#Case when c=b=h=1 i.e. lambda = 1
Temperature_Plot(c=1, b=1, h=1, T0=[1]) #Temperature 'explodes'

#Case when lambda = 2
Temperature_Plot(c=1, b=1, h=2, T0=[1]) #Temperature explodes but at a later time

#Case when lambda large
#Choose h=10, c=b=1 i.e. lambda = 10
Temperature_Plot(c=1, b=1, h=10, T0=[1]) #Temperature decays

#Lambda = e
Temperature_Plot(c=1, b=1, h=np.exp(1), T0=[1.01]) #Behaviour depends on initial condition

#Plot temperature over time (small lambda)
fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(Temp_Sol(1,1,1,[1])[0], Temp_Sol(1,1,1,[1])[1][0])
ax2.plot(Temp_Sol(1,1,2,[0.5])[0], Temp_Sol(1,1,2,[0.5])[1][0], color='green')
ax1.set_ylabel('Temperature')
ax2.set_ylabel('Temperature')
ax2.set_xlabel('Time (s)')
ax1.set_title('λ = 1, T(0) = 1')
ax2.set_title('λ = 2, T(0) = 0.5')
plt.tight_layout()

#Plot temperature over time (critical lambda)
fig, (ax1, ax2) = plt.subplots(2, figsize = (8,8))
#ax1.plot(Temp_Sol(1,1,np.exp(1),[1])[0], Temp_Sol(1,1,np.exp(1),[1])[1][0], color='green')
ax1.plot(Temp_Sol(1,1,np.exp(1),[0.75])[0], Temp_Sol(1,1,np.exp(1),[0.75])[1][0], color='blue')
ax2.plot(Temp_Sol(1,1,np.exp(1),[1.02])[0], Temp_Sol(1,1,np.exp(1),[1.02])[1][0], color='red')
ax1.set_ylabel('Temperature')
ax2.set_ylabel('Temperature')
ax2.set_xlabel('Time (s)')
ax1.set_title('T(0) = 0.75')
ax2.set_title('T(0) = 1.02')
plt.tight_layout()
ax1.grid()
ax2.grid()




