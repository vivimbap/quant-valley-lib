# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 11:07:31 2017

@author: Victor
"""

from math import log, exp, cos, pi, sqrt
import numpy as np
from simulation_methods import Euler_Maruyama

class Parametres:
    
    def __init__(self, x_init=100, 
                 diff_sigma=(lambda t=0,x=0: 0), diff_mu=(lambda t=0,x=0: 0),
                 bs_sigma=0, bs_mu=0, 
                 ou_a=0, ou_b=0, ou_sigma=0,
                 seas_a=0, seas_b=0, seas_c=0, seas_tau=0,
                 demand=(lambda t=0: 0)):
        
        self.x_init = x_init
        self.diff_sigma = diff_sigma
        self.diff_mu = diff_mu
        self.bs_sigma = bs_sigma
        self.bs_mu = bs_mu
        self.ou_a = ou_a
        self.ou_b = ou_b
        self.ou_sigma = ou_sigma
        self.seas_a = seas_a
        self.seas_b = seas_b
        self.demand = demand
        self.seas_c = seas_c
        self.seas_tau = seas_tau

class Diffusion:
    
    def __init__(self, param=Parametres()):
        self.Param = param

#BLACK-SCHOLES
class Black_Scholes:
    
    def __init__(self, param=Parametres()):
        self.Param = param
    
    #simulation_method / str / None, 'euler_maruyama' /default = None 
    def Trajectory(self, simulation_method='analytical', n=1000, t=0, T=1):
        
        if simulation_method=='analytical':
            delta_time = (T-t)/float(n)
            return self.Param.x_init*np.cumprod(np.exp(np.insert((self.Param.bs_mu - 0.5*self.Param.bs_sigma**2)*delta_time + np.random.normal(loc=0, scale=np.sqrt(delta_time)*self.Param.bs_sigma, size=n-1), 0, 0)))
        
        if simulation_method=='euler_maruyama':
            return Euler_Maruyama(X_0=self.Param.x_init, t_0=t, T=T, b=(lambda t,x: self.Param.bs_mu*x), sigma=(lambda t,x: self.Param.bs_sigma*x), n=n).integrate()
        
        else:
            print('must pass a valid method')
            

class Ornstein_Uhlenbeck:
    
    def __init__(self, param=Parametres()):
        self.Param = param
    
    def Trajectory(self, simulation_method='analytical', n=1000, t=0, T=1):
        
        if simulation_method=='analytical':
            delta = (T-t)/float(n)
            dt = self.Param.x_init*np.exp(-self.Param.ou_a*np.linspace(t,T,n)) + self.Param.ou_b*(1-np.exp(-self.Param.ou_a*np.linspace(t,T,n)))
            std = (self.Param.ou_sigma/sqrt(2*self.Param.ou_a))*sqrt(1-exp(-2*self.Param.ou_a*delta))
            dw = np.cumsum(np.random.normal(scale=std, size=n))
            return dt + dw 
        
        if simulation_method=='euler_maruyama':
            return Euler_Maruyama(X_0=self.Param.x_init, t_0=t, T=T, b=(lambda t,x: self.Param.ou_a*(self.Param.ou_b-x)), sigma=(lambda t,x: self.Param.ou_sigma), n=n).integrate()
        
        else:
            print('must pass a valid method')
            
class Lucia_Schwartz_2002:
    
    def __init__(self, param=Parametres()):
        self.Param = param
        
    def Demand(self, t):
        return self.Param.demand(t)
    
    def Seasonality(self, t):
        return self.Param.seas_a + self.Param.seas_b*self.Demand(t) + self.Param.seas_c*cos(((t-self.Param.seas_tau)*2*pi))
        
    def Trajectory(self, simulation_method='analytical', n=730, t=0, T=2):
        
        if simulation_method=='analytical':
            y_init = log(self.Param.x_init/float(exp(self.Seasonality(t))))
            return np.exp(Ornstein_Uhlenbeck(Parametres(x_init=y_init, ou_a=self.Param.ou_a, ou_sigma=self.Param.ou_sigma)).Trajectory(simulation_method='analytical', n=n, t=t, T=T) + np.vectorize(self.Seasonality)(np.linspace(t,T,n)))

        if simulation_method=='euler_maruyama':
            y_init = log(self.Param.x_init/float(exp(self.Seasonality(t))))
            return np.exp(Euler_Maruyama(X_0=y_init, b=(lambda t,x: -self.Param.ou_a*x), sigma=(lambda t,x: self.Param.ou_sigma), t_0=t, T=T, n=n).integrate() + np.vectorize(self.Seasonality)(np.linspace(t,T,n)))
        
        else:
            print('must pass a valid method')
    
class Merton_Jump:
    
    def __init__(self, param=Parametres()):
        self.Param = param






