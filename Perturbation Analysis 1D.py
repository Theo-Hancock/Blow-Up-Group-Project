# -*- coding: utf-8 -*-
"""
Created on Sat Mar  8 01:54:41 2025

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
Lambda = 0.01            #Non-dimensionalised reaction rate constant
r = delta_t/delta_x**2   #Need r = 1/2 for system to be numerically stable
k = Lambda*delta_t

#Tf = Qt, Q > 0 constant
Q = 2.88

#Boundary conditions at t=0
left_boundary_condition = 0
right_boundary_condition = 0

#Initial state for t=0 (without perturbation)
initial_state = []
for s in range(0,number_of_points):
    initial_state.append(0)
initial_state[0] = left_boundary_condition
initial_state[-1] = right_boundary_condition
initial_state = tuple(initial_state) #To prevent initial state from being modified


#Add asymmetric perturbation
np.random.seed(42)
sigma = 0.25        #Standard deviation
p = 0.5             #Perturbation size
x = np.linspace(-1, 1, 101)
perturb = p * (1/(sigma*np.sqrt(2*np.pi))) * np.exp(-0.5 * ((x-0)/sigma)**2)  # Gaussian.


#Run this block and set the boundary conditions back to 0 to reset the system
state = list(initial_state + perturb)
time = 0
plotsplease=[]
#point_values = [[],[]]  #Stores temperature values for a set point over time
plotsplease.append((list(state), time))
#point_values[0].append(state[97])
#point_values[1].append(0)


#Function which creates a new state and replaces the old one. 
#10k timesteps = 1 second (nondimensionalised time)
def update_state():
    new_state = []
    new_state.append(left_boundary_condition)
    for i in range(1,number_of_points-1):
        new_state.append(r*state[i-1] + r*state[i+1] + (1-2*r)*state[i] + k*math.exp(state[i]))
    new_state.append(right_boundary_condition)
    for d in range(0,number_of_points):
        state[d] = new_state[d]
    return state

#Parameter time_gap controls how often to take 'snapshots' of system
#10k timesteps = 1 second (nondimensionalised time)
def run_solution(n, time_gap):
    for i in range(0,n):
        global time
        global left_boundary_condition
        global right_boundary_condition
        left_boundary_condition += Q*delta_t
        right_boundary_condition += Q*delta_t
        update_state()
        state_entry = tuple(state)
        time += 1
        if time % time_gap == 0:
            plotsplease.append((state_entry, time))
        #point_values[0].append(state_entry[97])
        #point_values[1].append(time*delta_t)
        
        
#Simulate the system
run_solution(25000,5000)
run_solution(40, 40)
run_solution(6,25046)

#Save results
state_save = (state, time)

#Plot time evolution of system
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize = (12,6))
for m in range(0,len(plotsplease)):
    ax1.plot(x, plotsplease[m][0], label = f'τ = {round(plotsplease[m][1]/10000, 4)}')
ax2.plot(x, state, label = f'τ = {round(time/10000, 4)}')
ax1.set_xlabel("x", fontsize=14)
ax1.set_ylabel("u", fontsize=14)
ax2.set_xlabel("x", fontsize=14)
ax2.set_ylabel("u", fontsize=14)
ax1.set_ylim(0,20)
ax1.legend()
ax2.legend(prop={'size':13})
plt.tight_layout()
plt.grid()


#p = 0.01, sigma = 0.1
#Q = 5      No effect
#Q = 7      No perturbation -> Blowup at t=1.2567, x = 0.82. Perturbation little effect
#Q = 8      No perturbation -> Blowup at t=1.1273, x = 0.84. Perturbation little effect
#Q = 9      No perturbation -> Blowup at t=1.0231, x = 0.86. Perturbation has effect
#Q = 10     No perturbation -> Blowup at t=0.9372, x = 0.86. Perturbation has effect
#Q = 15     No perturbation -> Blowup at t=0.664, x = 0.9. Perturbation has effect
#Q = 20     No perturbation -> Blowup at t=0.5170, x = 0.92. Perturbation has effect


#p = 0.1, sigma = 0.1
#Q = 5      No effect
#Q = 8      Little effect
#Q = 9      Has effect

#p = 0.5, sigma = 0.1
#Q = 5      No symmetry breaking but blowup happens slightly earlier
#Q = 8      Perturbation has some effect, blowup slightly earlier

#p = 0.1, sigma = 0.02
#Q = 5      No effect
#Q = 8      No effect
#Q = 9      Has effect

#p = 0.1, sigma = 0.25
#Q = 5      No effect
#Q = 8      Little effect
#Q = 9      Has effect

#p = 0.5, sigma = 0.25
#Q = 5      Has effect
#Q = 3      No effect
#Q = 4      No effect


#Attempt to merge two blowup points into one
#p = 0.5, sigma = 0.25      Works for Q=2.87, but not for Q = 2.88
#p = 1, sigma = 0.25        Same as above
#p = 1, sigma = 0.5         Same as above
#p = 2, sigma = 0.5         Same as above
#p = 2, sigma = 0.25        Same as above


#Plot blowup state for different parameter values
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize = (12,6))
ax1.plot(x, state_save[0], label = f'τ = {round(state_save[1]*delta_t, 4)}')
ax1.set_xlabel("x", fontsize=14)
ax1.set_ylabel("u", fontsize=14)
#ax1.set_ylim(0,20)
ax1.set_title('Q = 8', fontsize=16)
ax1.legend(prop={'size':13})
ax1.grid()

#Plot point temperature over time
ax2.plot(x, state, label = f'τ = {round(time*delta_t, 4)}')
ax2.set_xlabel('x', fontsize=14)
ax2.set_ylabel('u', fontsize=14)
#ax2.set_ylim(0,20)
ax2.set_title('Q = 9', fontsize=16)
ax2.legend(prop={'size':13}, loc = 'upper center')
ax2.grid()
plt.tight_layout()





     
        