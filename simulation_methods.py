# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 11:06:44 2017

@author: Victor
"""

#change

import numpy as np

class Euler_Maruyama:
    
    def __init__(self, X_0=100, t_0=0, T=1, b=(lambda t,x : 0), sigma=(lambda t,x : 0), n=100):
        self.X_0=X_0
        self.t_0=t_0
        self.T=T
        self.b=b
        self.sigma=sigma
        self.n=n
        self.delta = (self.T-self.t_0)/float(n)
        self.list_t=np.linspace(self.t_0, self.T, num=self.n, endpoint=True, dtype=float)
        
    def integrate(self):
        results = np.zeros(self.n)
        results[0] = self.X_0
        for i in range(1,self.n):
            results[i] = results[i-1] + self.delta*self.b(self.list_t[i-1], results[i-1]) + self.sigma(self.list_t[i-1], results[i-1])*np.random.normal(scale=np.sqrt(self.delta))
        return results
