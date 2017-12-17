# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 21:05:11 2017

@author: Victor

"""

from math import sqrt 
import numpy as np

#RETURN TRUE IF n IS PRIME, FALSE OTHERWISE
def is_prime(n):
    
    #Two is prime
    if n==2:
        return True
    
    #If it's even then it's not prime
    if n%2==0:
        return False 
    
    #Check wether any odd number between 3 and sqrt(n) evenly divides n
    maxi = int(sqrt(n))
    d = 3
    
    while (n%d != 0) and (d<=maxi):
        d+=2
    
    #If none did evenly divide n then it's prime
    if d>maxi:
        return True
    #Else it's not prime
    else:
        return False



#RETURNS A NUMPY ARRAY CONTAINING THE n FIRST PRIMES 
def first_primes(n):
    
    #n first primes will be stored in 'ans' array
    ans = np.empty(n, dtype=int)
    
    #2 is the first prime so we add it to 'ans'
    if n!=0:
        ans[0]=2
    
    #variable 'count' counts how many primes have been found so far
    count = 1
    #variable 'aux' contains the next integer that will be tested
    aux = 3
    
    while count<n:
        if is_prime(aux):
            ans[count]=aux
            count+=1
        aux+=2
    
    return ans



#RETURNS A NUMPY ARRAY CONTAINING ALL PRIMES LESSER OR EQUAL TO n
def all_primes_leq(n):
    
    #we use Python's list comprehension feature for this one
    return np.array([p for p in range(2,n+1) if is_prime(p)])



#RETURNS A NUMPY ARRAY CONTAINING THE p-ADIC DECOMPOSITION OF n
#ie IF n = a_0 + a_1*p + a_2*p^2 + ... + a_r*p^r , with a_r!=0
#THEN IT RETURNS [a_0, a_1, a_2, ..., a_r]
def expansion(n,p=10):
    
    decomp = np.array([])
    a = p-1
    
    while n>p-1:
        a = n//p
        b = n%p
        n=a
        decomp = np.append(decomp, b)
        
    decomp = np.append(decomp, a)
    return decomp

def integrate_expansion(x,p=10):
    res = 0
    for i in range(len(x)):
        res = res + x[i]*(p**i)
    return res

#RETURNS A NUMPY ARRAY CONTAINING THE p-ADIC DECOMPOSITION OF x 
#IE IF x = x_1*p^(-1) + x_2*p^(-2) + ...
#THEN IT RETURNS [x_1, ..., x_20]
def expansion_01(x,p=10,size=20):
    
    decomp = np.empty(size)
    pow_p=p
    
    for i in range(size):
        k=0
        while (k+1)*(1/pow_p) <= x:
            k+=1
        decomp[i] = k
        x = x - k*(1/pow_p)
        pow_p*=p
        
    return decomp

def integrate_expansion_01(x,p):
    res = 0
    for i in range(len(x)):
        res = res + x[i]*(p**(-(i+1)))
    return res
    
def kakutani_adding_machine(x,y,p,size=20):
    
    if (type(x)==float):
        x = expansion_01(x,p,size)
        
    if (type(y)==float):
        y = expansion_01(y,p,size)
    
    results = np.zeros(size)
    epsilon = np.array([])
    
    results[0] = (x[0] + y[0])%p
    
    if (x[0] + y[0]<p):
        epsilon = np.append(epsilon, 0)
        
    else:
        epsilon = np.append(epsilon, 1)
    
    for i in range(1,size):
        results[i] = (epsilon[i-1] + x[i] + y[i])%p
        if (epsilon[i-1]+x[i]+y[i]<p):
            epsilon = np.append(epsilon, 0)
        else:
            epsilon = np.append(epsilon, 1)
    
    return results

    
    

        
        
        
        
        
        
        
        
        
        
        