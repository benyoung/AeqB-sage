#!/usr/bin/env python

def polyfit(L):
    v = [[i**j for j in range(len(L))] for i in range(len(L))]
    c = Matrix(v).solve_right(vector(L))
    print sum([c[k]*x^k for k in range(len(L))]) 
