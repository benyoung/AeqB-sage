def gosper_sum(t,n):
    # Takes an expression and the variable for which the expression takes as input, and returns
    # an expression for the resulting hypergeometric term that satisfies Gosper's algorithm
    r = (t(n=n+1)/t).full_simplify()
    P.<n,h> = PolynomialRing(QQ,'n,h')
    factorized = gosper_step2(r,n,h)
    a = P(factorized[0])
    b = P(factorized[1])
    c = P(factorized[2])
    x = gosper_step3(a,b,c,n)
    if x == "No polynomial solution exists":
        return "Could not find term that satisfies Gosper's algorithm."
    return b(n=n-1)*x/c*t
