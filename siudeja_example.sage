#!/usr/bin/env python

# n and l are as in your problem
# T is something like (abs(y)/abs(z))**2
# a,b,c,d,e,f are the coefficients of the recurrence
var("n,l, T, a,b,c,d,e,f")

# your function that you want to sum
def F(n,l):
    return binomial(2*n,l) * binomial(2*n-2*l, n-l) * T**l

# ratios of shifted terms in your function, simplified
def Q(a,b):
    return (F(n + a, l + b) / F(n, l)).factorial_simplify().factorial_simplify()

# This is the recurrence which your sum will satisfy; we will solve for
# a,b,c,d,e,f.  a,b,c,d,e,f will not depend on l.  Note Q(0,0) is 1.
lhs_long = a*Q(0,0) + b*Q(1,0) + c*Q(2,0) + d*Q(0,1) + e*Q(1,1) + f*Q(2,1)
lhs_poly = lhs_long().factor().numerator()
the_coeffs = lhs_poly.coefficients(l)

L = []
for C in the_coeffs:
    L.append([C[0].coeff(variable) for variable in [a,b,c,d,e,f]])
M = matrix(L)
soln = M.right_kernel().basis()[0]


# This will run for a long time.  When it is done you know the values of a,b,c,d,e,f
# above.  Now sum on l.  You can sum over all integers l, since your summand is zero
# outside the range [0,n].  Since a,b,c,d,e,f do not depend on l, then the result is
# the recurrence that your sum satisfies.
for x in soln:
    if x == 0:
        print "zero"
    else:
        print factor(x)
