def resultant_roots(x,y,n,h):
	return y(n=n+h).resultant(x,n).univariate_polynomial().roots()

def gosper_step2(r,n,h):
	# Takes as input a multivariate polynomial r defined in the polynomial variables n and
	# h. Returns polynomials that satisfy step 2 of Gosper's algorithm, as outlined in 
	# "A=B."
	x = r.numerator()
	y = r.denominator()

	S = resultant_roots(x,y,n,h)
	S = [item[0] for item in S if item[0] > 0]

	z = 1
	for j in range(0,len(S)):
		s = gcd(x,y(n=n+S[j]))
		x /= s
		y /= s(n=n-S[j])
		for i in range(1,S[j] + 1):
			z *= s(n=n-i)

	return (x,y,z)