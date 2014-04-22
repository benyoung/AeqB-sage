def gosper_sum(t,n):
    # Takes an expression and the variable for which the expression takes as input, and returns
    # an expression for the resulting hypergeometric term that satisfies Gosper's algorithm
    r = (t(n+1)/t).full_simplify()
    h = var('h')
    R = PolynomialRing(QQ,[n,h])
    factorized = gosper_step2(r,n,h)
    a = R(factorized[0])
    b = R(factorized[1])
    c = R(factorized[2])
    x = gosper_step3(a,b,c,n)
    if x == "No polynomial solution exists":
        return "Could not find term that satisfies Gosper's algorithm."
    b = R(b)
    return R(b(n-1,h))*x/c*t
