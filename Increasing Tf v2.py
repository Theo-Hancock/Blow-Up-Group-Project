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
Q = 2.8701

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
        new_state.append(r*state[i-1] + r*state[i+1] + (1-2*r)*state[i] + k*math.exp(state[i]))
    new_state.append(right_boundary_condition)
    if counter%time_gap == 0:
        plotsplease.append((new_state, counter))
    for d in range(0,number_of_points):
        state[d] = new_state[d]
    return state

#Simulate the system
for l in range(1,24001):
    counter += 1
    left_boundary_condition += Q*delta_t
    right_boundary_condition += Q*delta_t
    update_state(8000)
    
for l in range(1,1001):
    counter += 1
    left_boundary_condition += Q*delta_t
    right_boundary_condition += Q*delta_t
    update_state(1000)
    
for l in range(1,41):
    counter += 1
    left_boundary_condition += Q*delta_t
    right_boundary_condition += Q*delta_t
    update_state(40)
    
for l in range(1,10):
    counter += 1
    left_boundary_condition += Q*delta_t
    right_boundary_condition += Q*delta_t
    update_state(8)

for l in range(1,3):
    counter += 1
    left_boundary_condition += Q*delta_t
    right_boundary_condition += Q*delta_t
    update_state(25051)  
    


#Plot time evolution of system
x = np.linspace(-1, 1, number_of_points)

for m in range(0,len(plotsplease)):
    plt.plot(x,plotsplease[m][0], label = f'τ = {round(plotsplease[m][1]/10000, 4)}')
#plt.plot(x, state, label = f'τ = {round((counter)/10000,4)}')
plt.xlabel("x")
plt.ylabel("u")
plt.ylim(0,40)
plt.legend()
plt.show()

#Save results for Q=3 to separate array
plotsplease_dupe = tuple(plotsplease)

#Plot final state
plt.plot(x, state, label = f'τ = {round((counter)/10000, 3)}')
#plt.ylim(0,30)
plt.legend()
plt.show()

#Save final state for Q = 2.86 to separate array
state_dupe = (tuple(state), counter-1)

#Time to blowup for different values of Q (and Lambda=0.01):
    #Q = 1:     58273 (central blowup)
    #Q = 5:     16414 (two blowups)
    #Q = 3:     24233 (two blowups)
    #Q = 2.8    25521 (central blowup)
    #Q = 2.9    24857 (two blowups)
    #Q = 2.87   25052 (two blowups)
    #Q = 2.88   24987 (two blowups)
    #Q = 2.86   25118 (central blowup)
    
#Bifurcation at some point between 2.86 and 2.87



#Create plot with two subplots to compare Q = 3 and Q = 2.88
fig, (ax1, ax2) = plt.subplots(2, figsize = (8,8))
for m in range(0,len(plotsplease_dupe)):
    ax1.plot(x,plotsplease_dupe[m][0], label = f'τ = {round(plotsplease_dupe[m][1]/10000, 3)}')
ax1.set_xlabel("x")
ax1.set_ylabel("u")
ax1.set_ylim(0,30)
for m in range(0,len(plotsplease)):
    ax2.plot(x,plotsplease[m][0], label = f'τ = {round(plotsplease[m][1]/10000, 4)}')
ax2.set_xlabel("x")
ax2.set_ylabel("u")
ax2.set_ylim(0,40)
ax1.set_title('Q = 3')
ax2.set_title('Q = 2.87')
plt.tight_layout()
ax1.legend()
ax2.legend()
ax1.grid()
ax2.grid()


#Create plot with two subplots to compare Q = 2.87 and Q = 2.88
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize = (16,6))
ax1.plot(x,state_dupe[0], label = f'τ = {round(state_dupe[1]/10000, 3)}')
ax1.set_xlabel("x")
ax1.set_ylabel("u")
ax2.plot(x, state, label = f'τ = {round((counter-1)/10000, 3)}')
ax2.set_xlabel("x")
ax2.set_ylabel("u")
ax1.set_title('Q = 2.86')
ax2.set_title('Q = 2.87')
plt.tight_layout()
ax1.legend()
ax2.legend()
ax1.grid()
ax2.grid()
