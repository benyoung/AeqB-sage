def N(poly,i):
	
	# Shifts a polynomial in n to the left by i.
	 
	return poly.substitute({n:n+i})
	
def lc(poly):
	
	# Returns the leading coefficient of a polynomial in one variable.
	
	coefficients = poly.coefficients()
	return coefficients[len(coefficients) - 1][0]
	
def constant_coefficient(poly):
	
	# Returns the constant coefficient of a polynomial in any number of variables.
	
	for v in poly.variables():
		poly = poly.substitute({v:0})
	return int(poly)

def Poly(p,f, n = var('n')):

	# Input: A list p of polynomials p_i(n), for i = 0,1,...,r, and a polynomial f(n) with rational coefficients
	#
	# Output: The general polynomial solution of Ly = f, where L = sum_{i = 0}^r p_i(n)N^i

	r = len(p) - 1
	
	# Makes sure Sage knows these are polynomials in n
	f = f + 0*n
	for i in range(r+1): 
		p[i] = p[i] + 0*n

	# Step 1: Compute q_j's in array q

	q = [sum([binomial(i,j)*p[i] for i in range(j,r+1)]) for j in range(r+1)]

	# Step 2: Compute degree bound d for the polynomial solution y

	deg = [q[j].degree(n) for j in range(r+1)]
			
	b = max([deg[j] - j for j in range(r+1)])
	alpha = sum([lc(q[j])*n^j for j in range(r+1) if deg[j] - j == b])

	# d1 is the maximal integer root of the polynomial alpha
	
	if lc(alpha) == alpha: # checks if alpha is a constant polynomial
		d = max([f.degree(n) - b, -b - 1])
		
	else:
		c = constant_coefficient(alpha)
		
		if c == 0:
			d1 = 0
		else:
			m = prod([int(a[0].denominator()) for a in alpha.coefficients()])*c # all integral solutions must divide this
			m_divisors = [j for j in range(-abs(m), abs(m) + 1) if j != 0 and abs(m) % abs(j) == 0]
			d1 = max([j for j in m_divisors if alpha.substitute({n:j}) == 0])
			
		d = max([f.degree(n) - b, -b - 1, d1])
	
	if d < 0:
		
		return 'No polynomial solutions to this recurrence.'

	# Step 3: Using the method of undetermined coefficients, find all y(n) of the form y(n) = sum_{k = 0}^d c_k n^k satisfying Ly = f

	# create arbitrary polynomial y(n) = a0 + a1*n + ... + ad*n^d of degree d

	y = 0*n
	for i in range(0,d+1):
		v = var('a' + str(i))
		y  = y + v*n^i

	# Expand Ly = f, collect powers of n, and equate coefficients

	LHS = sum([p[j]*N(y,j) for j in range(r+1)])
	equations = [k[0] - constant_coefficient(k[0]) == constant_coefficient(k[0]) for k in (LHS - f).full_simplify().coefficients(n)]

	# Find solutions
	general_poly = y
	solutions = solve(equations, y.substitute({n:1}).variables())
	
	if len(solutions) > 0:
		coeff = solutions[0]
		
		# need a check that there is more than just the trivial solution...
		if type(solutions[0]) == list:
			
			for j in range(len(coeff)):
				general_poly = general_poly.substitute({coeff[j].lhs():coeff[j].rhs()})
				for k in coeff[j].rhs().variables():
					coeff[j] = coeff[j].substitute({k:1})
					y = y.substitute({coeff[j].lhs():coeff[j].rhs()})
			# return a general solution with coefficients given in terms of free variables and one with all free variables set to 1
			return general_poly.factor(), y
		
		else:
			return 'Only the trivial solution.'
		
	else:
		return 'No polynomial solutions to this recurrence.'
		
# Find polynomial solutions of the recurrence 3y(n+2) - ny(n+1) + (n-1)y(n) = 0
print(Poly([n-1,-n,3],0))
		
