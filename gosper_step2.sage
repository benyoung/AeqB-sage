def gosper_step2(r,n,h):
	# Takes as input a multivariate rational r defined in the polynomial variables n and
	# h. Returns polynomials that satisfy step 2 of Gosper's algorithm, as outlined in 
	# "A=B."    
	R = PolynomialRing(QQ,[n,h])
	x = R(r.numerator())    
	y = R(r.denominator())

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