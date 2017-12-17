# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 10:44:10 2017

@author: Victor
"""

from stochastic_process import Parametres, Black_Scholes, Ornstein_Uhlenbeck, Lucia_Schwartz_2002
import matplotlib.pyplot as plt
import datetime 

figsize=(5,5)

number_simulation = 10

#BLACK_SCHOLES

param = Parametres(x_init=100, bs_mu=0.03, bs_sigma=0.015)

#Simulation de trajectoires pour method='analytical'
method='analytical'
start = datetime.datetime.now()
plt.figure(figsize=figsize)
plt.title('Black_Scholes, {} simulations'.format(number_simulation))
for i in range(number_simulation):
    traj = Black_Scholes(param=param).Trajectory(simulation_method=method)
    plt.plot(traj)
print('{} simulation (method=None) time = {}'.format(number_simulation, datetime.datetime.now()-start))
plt.show()

#Simulation de trajectoires pour method='euler_maruyama'
method='euler_maruyama'
start = datetime.datetime.now()
plt.figure(figsize=figsize)
plt.title('Black_Scholes, {} simulations'.format(number_simulation))
for i in range(number_simulation):
    traj = Black_Scholes(param=param).Trajectory(simulation_method=method)
    plt.plot(traj)
print('{} simulation (method=None) time = {}'.format(number_simulation, datetime.datetime.now()-start))
plt.show()

#ORNSTEIN_UHLENBECK

param = Parametres(x_init=3, ou_a=2, ou_b=3, ou_sigma=0.15)

#Simulation de trajectoires pour method='analytical'
method='analytical'
start = datetime.datetime.now()
plt.figure(figsize=figsize)
plt.title('Ornstein_Uhlenbeck, {} simulations'.format(number_simulation))
for i in range(number_simulation):
    traj = Ornstein_Uhlenbeck(param=param).Trajectory(simulation_method=method)
    plt.plot(traj)
print('{} simulation (method=None) time = {}'.format(number_simulation, datetime.datetime.now()-start))
plt.show()

#Simulation de trajectoires pour method='euler_maruyama'
method='euler_maruyama'
start = datetime.datetime.now()
plt.figure(figsize=figsize)
plt.title('Ornstein_Uhlenbeck, {} simulations'.format(number_simulation))
for i in range(number_simulation):
    traj = Ornstein_Uhlenbeck(param=param).Trajectory(simulation_method=method)
    plt.plot(traj)
print('{} simulation (method=None) time = {}'.format(number_simulation, datetime.datetime.now()-start))
plt.show()

#LUCIA_SCHWARTZ_2002

#Plot la saisonalité

param = Parametres(x_init=100, ou_a=5.84, ou_sigma=0.5, seas_a=4.86, seas_b=0.09, seas_c=0.306, seas_tau=0.836)

#Simulation de trajectoires pour method='analytical'
method='analytical'
start = datetime.datetime.now()
plt.figure(figsize=figsize)
plt.title('Lucia_Schwartz_2002, {} simulations'.format(number_simulation))
for i in range(number_simulation):
    traj = Lucia_Schwartz_2002(param=param).Trajectory(simulation_method=method)
    plt.plot(traj)
print('{} simulation (method=None) time = {}'.format(number_simulation, datetime.datetime.now()-start))
plt.show()

method='euler_maruyama'
start = datetime.datetime.now()
plt.figure(figsize=figsize)
plt.title('Lucia_Schwartz_2002, {} simulations'.format(number_simulation))
for i in range(number_simulation):
    traj = Lucia_Schwartz_2002(param=param).Trajectory(simulation_method=method)
    plt.plot(traj)
print('{} simulation (method=None) time = {}'.format(number_simulation, datetime.datetime.now()-start))
plt.show()
