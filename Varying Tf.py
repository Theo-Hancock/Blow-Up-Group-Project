# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 14:29:24 2025

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
r = delta_t/delta_x**2   #alpha*delta_t/delta_x**2, need r < 1/2 for system to be numerically stable
k = Lambda*delta_t       #beta*delta_t

#Q > 0 constant, rate of increase of furnace temperature
Q = 1
#P > 0 constant, rate of decrease of furnace temperature
P = 1
#τ_f: time at which we reverse the furnace temperature (in timesteps)
t_f = 51950

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
time = 0
plotsplease=[]
point_values = [[],[]]  #Stores temperature values for a set point over time
plotsplease.append((list(initial_state), time))
point_values[0].append(initial_state[50])
point_values[1].append(0)

#Function which updates the state each timestep. 
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

#Function which simulates the system for n timesteps, updating the boundary conditions each step
#Parameter time_gap controls how often to take 'snapshots' of system
def run_solution(n, time_gap):
    for l in range(0,n):
        global time
        global left_boundary_condition
        global right_boundary_condition
        if time < t_f:
            left_boundary_condition += Q*delta_t
            right_boundary_condition += Q*delta_t
        elif time < (1+(Q/P))*t_f:
            left_boundary_condition -= P*delta_t
            right_boundary_condition -= P*delta_t
        else:
            left_boundary_condition = 0
            right_boundary_condition = 0
        update_state()
        state_entry = tuple(state)
        time += 1
        if time % time_gap == 0:
            plotsplease.append((state_entry, time))
        point_values[0].append(state_entry[50])
        point_values[1].append(time*delta_t)
        
        
        
#Simulate the system
run_solution(t_f,t_f)
run_solution(140000,140000)  
run_solution(15,23415)
run_solution(6000,6000)
run_solution(10000,10000)


#Plot time evolution of system
x = np.linspace(-1, 1, number_of_points)
#plt.plot(x,state)

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize = (20,8))
for m in range(1,len(plotsplease)):
    ax1.plot(x,plotsplease[m][0], label = f'τ = {round(plotsplease[m][1]/10000, 4)}')
ax1.set_xlabel("x")
ax1.set_ylabel("u")
#ax1.set_ylim(0,20)
ax1.legend(prop={'size':13})
ax1.grid()

#Plot point temperature over time
ax2.plot(point_values[1], point_values[0], color='blue')
ax2.axvline(x=t_f/10000, linestyle='dashed', label = f'τ_f = {round(t_f/10000, 4)}', color='orange')
ax2.set_xlabel('τ')
ax2.set_ylabel('u')
ax2.set_ylim(0,20)
ax2.legend(prop={'size':25})
ax2.grid()
plt.tight_layout()

#Q = P = 5        
#t_f = 1.2      No blowup
#t_f = 1.35     Blowup
#t_f = 1.3      No blowup
#t_f = 1.32     No blowup
#t_f = 1.34     Blowup (t=1.9934)
#t_f = 1.33     Blowup (t=2.3415)
#t_f = 1.325    No blowup
#t_f = 1.328    No blowup
#t_f = 1.329    No blowup (peak temp around t=2.1)

#Critical value between 1.329 and 1.33


#Q = 5, P = 10
#t_f = 1.5      Blowup (t=1.6980)
#t_f = 1.4      No blowup
#t_f = 1.45     Blowup (t=1.7503)
#t_f = 1.43     Blowup (t=1.7981)
#t_f = 1.41     Blowup (t=1.9342)
#t_f = 1.405    No blowup
#t_f = 1.407    Blowup (t=2.0294)
#t_f = 1.406    No blowup

#Critical value between 1.406 and 1.407

#Save central temperature values for t_f = 1.406 (Q=5, P=10)
point_values_save = tuple(point_values)


#Q = 1, P = 1
#t_f = 5.0      No blowup
#t_f = 5.5      Blowup
#t_f = 5.25     Blowup (t=6.1582)
#t_f = 5.125    No blowup
#t_f = 5.2      Blowup (t.6.6413)
#t_f = 5.16     No blowup
#t_f = 5.18     No blowup
#t_f = 5.19     No blowup
#t_f = 5.195    No blowup
#t_f = 5.198    Blowup (t=6.7318)
#t_f = 5.197    Blowup (t=6.8008)
#t_f = 5.196    Blowup

#Critical value between 5.195 and 5.196

             

#Plot central temperature over time for different P,Q values
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize = (20,8))
ax1.plot(point_values_save[1], point_values_save[0], color='blue')
ax1.axvline(x=1.406, linestyle='dashed', label = f'τ_f = {round(1.406, 4)}', color='orange')
ax1.set_xlabel("τ")
ax1.set_ylabel("u")
ax1.set_title('Q = 5, P = 10', fontsize=25)
ax1.set_ylim(0,10)
ax1.legend(prop={'size':25})
ax1.grid()
ax2.plot(point_values[1], point_values[0], color='blue')
ax2.axvline(x=t_f/10000, linestyle='dashed', label = f'τ_f = {round(t_f/10000, 4)}', color='orange')
ax2.set_xlabel('τ')
ax2.set_ylabel('u')
ax2.set_title('Q = 1, P = 1', fontsize=25)
ax2.set_ylim(0,10)
ax2.legend(prop={'size':25})
ax2.grid()
plt.tight_layout()
            
        

