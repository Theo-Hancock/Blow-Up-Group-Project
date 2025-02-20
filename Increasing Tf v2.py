# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 09:47:37 2025

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
#alpha = float(8.21e-08) #Diffusion constant
#beta = float(9.34e-19)  #Pre-exponential term
#gamma = 0.072           #Exponential term
Lambda = 0.01            #Non-dimensionalised constant
r = delta_t/delta_x**2   #alpha*delta_t/delta_x**2, need r = 1/2 for system to be numerically stable
k = Lambda*delta_t       #beta*delta_t

#Tf = Qt, Q > 0 constant
Q = 5

#Boundary conditions at t=0
left_boundary_condition = 0
right_boundary_condition = 0

#Initial state for t=0
initial_state = []
for s in range(0,number_of_points):
    initial_state.append(0)
initial_state[0] = left_boundary_condition
initial_state[-1] = right_boundary_condition
initial_state = tuple(initial_state) #To prevent initial state from being modified

#Run this block and set the boundary conditions back to 0 to reset the system
state = list(initial_state)
counter = 0
plotsplease=[]
plotsplease.append((list(initial_state), counter))

#Function which creates a new state and replaces the old one. 
#Parameter time_gap controls how often to take 'snapshots' of system
#10k timesteps = 1 second (nondimensionalised time)
def update_state(time_gap):
    new_state = []
    new_state.append(left_boundary_condition)
    for i in range(1,number_of_points-1):
        if i == number_of_points-2:
            new_state.append(r*state[i-1] + r*right_boundary_condition + (1-2*r)*state[i] + k*math.exp(state[i]))
        else:
            new_state.append(r*state[i-1] + r*state[i+1] + (1-2*r)*state[i] + k*math.exp(state[i]))
    new_state.append(right_boundary_condition)
    if counter%time_gap == 0:
        plotsplease.append((new_state, counter))
    for d in range(0,number_of_points):
        state[d] = new_state[d]
    return state

#Simulate the system
for l in range(1,50001):
    counter += 1
    left_boundary_condition += Q*delta_t
    right_boundary_condition += Q*delta_t
    update_state(15000)
    
for l in range(1,8001):
    counter += 1
    left_boundary_condition += Q*delta_t
    right_boundary_condition += Q*delta_t
    update_state(58000)

for l in range(1,272):
    counter += 1
    left_boundary_condition += Q*delta_t
    right_boundary_condition += Q*delta_t
    update_state(100)


#Plot time evolution of system
x = np.linspace(-1, 1, number_of_points)

for m in range(0,len(plotsplease)):
    plt.plot(x,plotsplease[m][0], label = f'τ = {round(plotsplease[m][1]/10000, 2)}')
plt.plot(x, state, label = f'τ = {round(counter/10000,3)}')
plt.xlabel("x")
plt.ylabel("u")
plt.ylim(0,30)
plt.legend()
plt.show()


#Plot final state
plt.plot(x, state, label = f'τ = {counter}')
plt.ylim(0,30)
plt.legend()
plt.show()


#Time to blowup for different values of Q (and Lambda=0.01):
    #Q = 1:  58273 (central blowup)
