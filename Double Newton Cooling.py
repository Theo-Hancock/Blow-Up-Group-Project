# -*- coding: utf-8 -*-
"""
Created on Thu Mar  6 15:27:13 2025

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
Lambda = 0.01            #Non-dimensionalised constant
r = delta_t/delta_x**2   #Need r < 1/2 for system to be numerically stable
K = Lambda*delta_t       #Not to be confused with Newton cooling constant k
k = 1                    #Newton cooling constant


#Q > 0 constant, rate of increase of furnace temperature
Q = 5
#τ_f: time at which we remove from furnace and apply double Robin boundary condition
t_f = 14500
#Note - assuming ambient temperature u_0 = 0 for Newton cooling


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


#Function which updates the state each timestep
def update_state():
    new_state = []
    for i in range(0, number_of_points):
        if i == 0:
            if time >= t_f:
                new_state.append(2*r*state[i+1] + (1 - 2*r - 2*r*k*delta_x)*state[i] + K*math.exp(state[i]))
            else:
                new_state.append(left_boundary_condition)
        elif i in range(1,number_of_points-1):
            new_state.append(r*state[i-1] + r*state[i+1] + (1-2*r)*state[i] + K*math.exp(state[i]))
        elif i == number_of_points-1:
            if time >= t_f:
                new_state.append(2*r*state[i-1] + (1 - 2*r - 2*r*k*delta_x)*state[i] + K*math.exp(state[i]))
            else:
                new_state.append(right_boundary_condition)
    for d in range(0,number_of_points):
        state[d] = new_state[d]
    return state


#Function which simulates the system for n timesteps, updating the boundary conditions each step
#Parameter time_gap controls when to take 'snapshots' of system
def run_solution(n, time_gap):
    for l in range(0,n):
        global time
        global left_boundary_condition
        global right_boundary_condition
        if time < t_f:
            left_boundary_condition += Q*delta_t
            right_boundary_condition += Q*delta_t
        update_state()
        state_entry = tuple(state)
        time += 1
        if time % time_gap == 0:
            plotsplease.append((state_entry, time))
        point_values[0].append(state_entry[50])
        point_values[1].append(time*delta_t)
        
        
run_solution(40000,40000)

#Save central temperature result
point_values_save = tuple(point_values)


#Plot time evolution of system
x = np.linspace(-1, 1, number_of_points)


#Central temperature over time for different parameter values
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize = (20,8))
ax1.plot(point_values_save[1], point_values_save[0], color='blue')
ax1.axvline(x=1.4499, linestyle='dashed', label = 'τ_f = 1.4499', color='orange')
ax1.set_xlabel('τ')
ax1.set_ylabel('u')
ax1.set_ylim(0,20)
ax1.legend(prop={'size':25})
ax1.grid()

ax2.plot(point_values[1], point_values[0], color='blue')
ax2.axvline(x=1.45, linestyle='dashed', label = 'τ_f = 1.45', color='orange')
ax2.set_xlabel('τ')
ax2.set_ylabel('u')
ax2.set_ylim(0,20)
ax2.legend(prop={'size':25}, loc='upper right')
ax2.grid()

plt.tight_layout()



#Temperature profile over time + central temperature over time
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize = (20,8))
for m in range(1,len(plotsplease)):
    ax1.plot(x,plotsplease[m][0], label = f'τ = {round(plotsplease[m][1]/10000, 4)}')
ax1.set_xlabel("x")
ax1.set_ylabel("u")
ax1.set_ylim(0,20)
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
        

#Results
#Q=5, k=1
#t_f = 1.5      Blowup (t=1.7539)
#t_f = 1.4      No blowup
#t_f = 1.45     Blowup (t=3.8305) !! look at central temperature plot
#t_f = 1.449    No blowup (also very intersting central temperature plot)
#t_f = 1.4495   No blowup
#t_f = 1.4498   No blowup
#t_f = 1.4499   No blowup (maybe use this value in the report alongside 1.45)

#Critical value between 1.4499 and 1.45


#Q=5, k=2
#t_f = 1.5      No blowup
#t_f = 1.6      Blowup (t=1.6562)
#t_f = 1.55     Blowup (t=1.7051)
#t_f = 1.52     Blowup (t=1.7874)
#t_f = 1.51     Blowup (t=1.8854)
#t_f = 1.505    Blowup (t=2.1054)
#t_f = 1.503    No blowup
#t_f = 1.504    No blowup

#Critical value between 1.504 and 1.505




