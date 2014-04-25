def gosper_sum(t,n):
    # Takes an expression and the variable for which the expression takes as input, and returns
    # an expression for the resulting hypergeometric term that satisfies Gosper's algorithm
    r = (t(n+1)/t).full_simplify()
    h = var('h')
    R = PolynomialRing(QQ,[n,h])
    
    # Find the polynomials a,b,c such that gcd(a(n),b(n+h) = 1 for any positive integer h
    factorized = gosper_step2(r,n,h)
    a = R(factorized[0])
    b = R(factorized[1])
    c = R(factorized[2])
    
    # Find the polynomial x solve a*x(n+1) - b(n-1)*x = c
    x = gosper_step3(a,b,c,n)
    if x == "No polynomial solution exists":
        return "Could not find term that satisfies Gosper's algorithm."
    b = R(b)
    
    # Return the desired hypergeometric term
    return R(b(n-1,h))*x/c*t
