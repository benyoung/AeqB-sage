def resultant_roots(x,y,n,h):
	return y(n=n+h).resultant(x,n).univariate_polynomial().roots()

def gosper_step2(r,n,h):
	x = r.numerator()
	y = r.denominator()
	
	S = resultant_roots(x,y,n,h)
	S = [item[0] for item in S if item[0] > 0]

	C = 1
	for j in range(0,len(S)):
		s = gcd(x,y(n=n+S[j]))
		x /= s
		y /= s(n=n-S[j])
		for i in range(1,S[j] + 1):
			C *= s(n=n-i)
	
	A = x
	B = y
	return (A,B,C)