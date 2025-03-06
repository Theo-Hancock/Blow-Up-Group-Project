# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 09:43:26 2025

@author: theoj
"""

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

#Constants
number_of_points = int(101)
delta_x = 2/number_of_points
delta_t = 1/10000
D_T = 1             #Diffusion constant - temperature
D_a = 1             #Diffusion constant - reactant
D = D_a/D_T
gamma = 1           #Reaction constant - temperature
mu = 1              #Reaction constant - reactant
r_T = delta_t/delta_x**2            
r_a = D * (delta_t/delta_x**2)      #Make sure r_T and r_a are both < 0.5

#Define function for spatially non-uniform initial reactant concentration
#def f(x):
#    return 0.5*x**2 + 0.5

#Initial amount of reactant (constant case)
initial_a = 1

#Coupled system
#   dT/dt = d^2T/dx^2 + (gamma * a * e^T)
#   da/dt = (D_a/D_T * d^2a/dx^2) - (mu * a * e^T)


#Boundary conditions for temperature and reactant concentration at t=0
left_boundary_temp = 0
right_boundary_temp = 0
left_boundary_reactant = 1
right_boundary_reactant = 1

#Initial state for temperature at t=0
initial_temp = []
for s in range(0,number_of_points):
    initial_temp.append(0)
initial_temp[0] = left_boundary_temp
initial_temp[-1] = right_boundary_temp
initial_temp = tuple(initial_temp) #To prevent initial state from being modified

#Initial state for reactant at t=0
initial_reactant = []
for s in range(0, number_of_points):
    initial_reactant.append(initial_a)
#x = np.linspace(-1, 1, number_of_points)
#for s in range(0,number_of_points):
#    initial_reactant.append(f(x[s]))
initial_reactant = tuple(initial_reactant) #To prevent initial state from being modified

#Run this block (and set the boundary conditions if needed) to reset the system
T_state = list(initial_temp)
a_state = list(initial_reactant)
time = 0
plotsplease=[]
plotsplease.append((list(initial_reactant), list(initial_temp), time))
central_values=[[],[],[]]
central_values[0].append(initial_reactant[50])
central_values[1].append(initial_temp[50])
central_values[2].append(0)


#Update reactant concentration each timestep
def update_a_state():
    new_state = []
    new_state.append(left_boundary_reactant)
    for i in range(1,number_of_points-1):
        k = delta_t * mu * a_state[i]
        new_state.append(r_a*a_state[i-1] + r_a*a_state[i+1] + (1-2*r_a)*a_state[i] - k*math.exp(T_state[i]))
    new_state.append(right_boundary_reactant)
    for d in range(0,number_of_points):
        a_state[d] = new_state[d]
    return a_state
                 
#Update temperature each timestep                
def update_T_state():
    new_state = []
    new_state.append(left_boundary_temp)
    for i in range(1,number_of_points-1):
        k = delta_t * gamma * a_state[i]
        new_state.append(r_T*T_state[i-1] + r_T*T_state[i+1] + (1-2*r_T)*T_state[i] + k*math.exp(T_state[i]))
    new_state.append(right_boundary_temp)
    for d in range(0,number_of_points):
        T_state[d] = new_state[d]
    return T_state


#Combined function to update reactant and temperature, add plots
def update_state(n_steps, time_gap):
    for l in range(0,n_steps):
        global time
        update_a_state()
        update_T_state()
        a_entry = tuple(a_state)
        T_entry = tuple(T_state)
        time += 1
        if time % time_gap == 0:
            plotsplease.append((list(a_entry), list(T_entry), time))
        if time % (1/(1000*delta_t)) == 0:
            central_values[0].append(a_entry[50])
            central_values[1].append(T_entry[50])
            central_values[2].append(time*delta_t)
                
            
            
#Simulate the system
update_state(20000,10000)

#Save results (for making plots with subplots for different parameter values)
central_values_save = tuple(central_values)
plotsplease_save = tuple(plotsplease)


#Plot time evolution of system
x = np.linspace(-1, 1, number_of_points)


#Plot with two subplots for different time values
fig, (ax1, ax2) = plt.subplots(2, figsize = (8,8))
ax1.plot(x, plotsplease[0][0], label = 'Reactant a')
ax1.plot(x, plotsplease[0][1], label = 'Temperature u')
ax1.set_xlabel("x")
ax1.set_title(f'τ = {round(plotsplease[0][2]*delta_t, 2)}')
ax2.plot(x, plotsplease[1][0], label = 'Reactant a')
ax2.plot(x, plotsplease[1][1], label = 'Temperature u')
ax2.set_xlabel("x")
ax2.set_title(f'τ = {round(plotsplease[1][2]*delta_t, 2)}')
plt.tight_layout()
ax1.legend()
ax2.legend()
ax1.grid()
ax2.grid()


#Plot central temperature and reactant over time
plt.plot(central_values[2], central_values[0], label = 'Reactant a(0,τ)')
plt.plot(central_values[2], central_values[1], label = 'Temperature u(0,τ)')
plt.xlabel('τ')
plt.ylim(0,1)
plt.tight_layout()
plt.legend()
plt.grid()


#Plots for comparing parameter values
fig, (ax1, ax2) = plt.subplots(2, figsize = (6,8))
ax1.plot(central_values_save[2], central_values_save[0], label = 'Reactant a(0,τ)')
ax1.plot(central_values_save[2], central_values_save[1], label = 'Temperature u(0,τ)')
ax1.set_xlabel('τ')
ax1.set_xlim(0,5)
ax1.set_title('γ = 5, μ = 1')
ax2.plot(central_values[2], central_values[0], label = 'Reactant a(0,τ)')
ax2.plot(central_values[2], central_values[1], label = 'Temperature u(0,τ)')
ax2.set_xlabel('τ')
ax2.set_xlim(0,5)
ax2.set_title('γ = 5, μ = 5')
plt.tight_layout()
ax1.legend()
ax2.legend()
ax1.grid()
ax2.grid()


#Two subplots - space/time plot and central value plot
fig, (ax1, ax2) = plt.subplots(2, figsize = (6,8))
for m in range(0,len(plotsplease)):
    ax1.plot(x, plotsplease[m][0], label = f'Reactant - τ = {round(plotsplease[m][2]*delta_t, 2)}')
for m in range(0,len(plotsplease)):
    ax1.plot(x, plotsplease[m][1], label = f'Temperature - τ = {round(plotsplease[m][2]*delta_t, 2)}')
ax1.set_xlabel('τ')
ax2.plot(central_values[2], central_values[0], label = 'Reactant a(0,τ)')
ax2.plot(central_values[2], central_values[1], label = 'Temperature u(0,τ)')
ax2.set_xlabel('τ')
ax2.set_xlim(0,2)
plt.tight_layout()
ax1.legend(prop={'size':8})
ax2.legend()
ax1.grid()
ax2.grid()






        

    

