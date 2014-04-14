#!/usr/bin/env python

def polyfitter_vec(v):
    if len(v) == 1:
        return v;
    wp = polyfitter_vec([v[i+1]-v[i] for i in range(len(v)-1)])
    wpp = [wp[i]/(i+1) for i in range(len(v)-1)]
    wpp.insert(0, v[0])
    return wpp

def polyfit(v, n=var('n')):
    w = polyfitter_vec(v)
    acc = 0
    for i in range(len(w)):
        acc += w[i] * falling_factorial(n, i)
    if acc == 0:
        return 0;
    return acc.expand()


