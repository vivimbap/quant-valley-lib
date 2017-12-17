# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 23:33:35 2017

@author: Victor
"""

import functions
from math import sqrt
import numpy as np

class Halton:
    
    def __init__(self, d, n=0):
        
        self.Dimension = d
        self.Size = n
        self.Values = np.empty((self.Size, self.Dimension))
        
        p = functions.first_primes(self.Dimension)
        
        for i in range(self.Size):
            for j in range(self.Dimension):
                
                decomp = functions.expansion(i+1, p[j])
                phi = 0
                pow_p = p[j]
                
                for k in range(len(decomp)):
                    phi += decomp[k]/pow_p
                    pow_p *= p[j]
                 
                self.Values[i,j] = phi
    
    def __getitem__(self, key):
        
        #If the key-th value of the Halton sequence has already been computed,
        #then it is returned.
        if (self.Size > key):
            return self.Values[key]
        #Else it computes the key-th value of the Halton sequence and returns 
        #that.
        else:
            p = functions.first_primes(self.Dimension)
            res = np.empty(self.Dimension)
            for j in range(self.Dimension):
                decomp = functions.expansion(key+1, p[j])
                phi = 0
                pow_p = p[j]
                for k in range(len(decomp)):
                    phi += decomp[k]/pow_p
                    pow_p *= p[j]
                res[j] = phi
            return res
        
    
    def resize(self, new_shape):
        return
        
class Kakutani:
    
    def __init__(self, d, n=0, x_init=None, y_init=None):
        
        self.Dimension = d
        self.Size = n
        
        p = functions.first_primes(self.Dimension)
        
        if (x_init==None):
            self.x = np.empty(d)
            for i in range(d):
                if (((i+1)!=3) and ((i+1)!=4)):
                    self.x[i] = 1/5
                else:
                    self.x[i] = (2*p[i]-1-sqrt((p[i]+2)**2+4*p[i]))/3
        else:
            self.x = x_init
        
        if (y_init==None):
            self.y = np.empty(d)
            for i in range(d):
                self.y[i] = 1/p[i]
                
        else:
            self.y = y_init
        
        self.Values = np.empty((self.Size, self.Dimension))
        self.Values[0] = self.x
        for i in range(1,self.Size):
            for j in range(self.Dimension):
                self.Values[i,j] = functions.integrate_expansion_01(functions.kakutani_adding_machine(self.Values[i-1,j],self.y[j],p[j]))
