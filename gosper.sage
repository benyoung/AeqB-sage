##########################################################################################
# DESCRIPTION: Follows the explanation of the complete Gosper's algorithm described in "A=B"
# by Marko Petkovsek, Herbert Wilf, and Doron Zeilberger.
#
# INPUT: Takes an expression that represents the hypergeometric term t(n) and the variable
# n for which t is defined over.
#
# OUTPUT: Returns an expression z_n for the resulting hypergeometric term that satisfies 
# Gosper's algorithm. I.E., sum(t(k),k,0,n-1) = z_n - z_0
#
# AUTHOR: Kevin Wilson, kwilson8@uoregon.edu
##########################################################################################

def gosper_sum(t,n):
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
