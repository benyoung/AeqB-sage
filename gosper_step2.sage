##########################################################################################
# DESCRIPTION: Follows the explanation of step 2 of Gosper's algorithm described in "A=B" by 
# Marko Petkovsek, Herbert Wilf, and Doron Zeilberger.
#
# INPUT: Takes as input a rational expression r defined in variables n and h, and the input      
# variables n and h.
#
# OUTPUT: Returns the polynomials a, b, and c such that r(n) = a(n)/b(n)*c(n+1)/c(n) and
# gcd(a(n),b(n+h)) = 1 for any non-negative integer h.
#
# AUTHOR: Kevin Wilson, kwilson8@uoregon.edu
##########################################################################################

def gosper_step2(r,n,h):
	R = PolynomialRing(QQ,[n,h])
	x = R(r.numerator())    
	y = R(r.denominator())
    
    # Calculate common roots of y(n+h) and x(n) for any non-negative integer h
	S = R(y(n+h,h)).resultant(x,R(n)).univariate_polynomial().roots()
	S = [item[0] for item in S if item[0] > 0]

	z = 1
	for j in range(0,len(S)):
		s = gcd(x,R(y(n+S[j],h)))
		x /= s
		y /= s(n-S[j],h)
		for i in range(1,S[j] + 1):
			z *= s(n-i,h)

	return (x,y,z)