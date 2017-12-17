# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 23:33:35 2017

@author: Victor
"""

from functions import decomposition, first_primes
import numpy as np

class Halton:
    
    def __init__(self, d, n=0):
        
        self.Dimension = d
        self.Size = n
        self.Values = np.empty((self.Size, self.Dimension))
        
        p = first_primes(self.Dimension)
        
        for i in range(self.Size):
            for j in range(self.Dimension):
                
                decomp = decomposition(i+1, p[j])
                phi = 0
                pow_p = p[j]
                
                for k in range(len(decomp)):
                    phi += decomp[k]/pow_p
                    pow_p *= p[j]
                 
                self.Values[i,j] = phi
                
                
                
                