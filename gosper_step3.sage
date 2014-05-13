##########################################################################################
# DESCRIPTION: Follows the explanation of step 3 of Gosper's algorithm described in "A=B" by 
# Marko Petkovsek, Herbert Wilf, and Doron Zeilberger.
#
# INPUT: Takes as input multivariate polynomials a, b, and c defined over variables including
# n, and input variable n.
#
# OUTPUT: Returns the polynomial x(n) such that a(n)*x(n+1) - b(n-1)*x(n) = c(n).
#
# AUTHOR: Kevin Wilson, kwilson8@uoregon.edu
##########################################################################################

def gosper_step3(a,b,c,n):
	deg_a = a.degree()
	deg_b = b.degree()
	deg_c = c.degree()
	lc_a = a.coefficients()[-1]
	lc_b = b.coefficients()[-1]
	
	# Determine degree of polynomial x
	# If degree of a doesn't equal degree of b or the the degrees are equal but the leading
	# coefficients aren't, then the degree of x is the degree of c minus the maximum of
	# the degrees of a and b since a(n)*x(n+1) - b(n-1)*x(n) = c(n).
	if (deg_a != deg_b) or lc_a != lc_b:
		d = (deg_c - max(deg_a,deg_b),0)
	# Otherwise, we examine the leading coefficients of the second power and use these
	# to determine the power of x.
	else:
		A = a.coefficients()[-2]
		B = b.coefficients()[-2]
		d = (deg_c - deg_a + 1, (B - A)/lc_a)
	
	if d[0] < 0 or int(d[1]) != d[1]:
		return "No polynomial solution exists"
	else:
		d = max(d[0],d[1])
		
	# Create the unknown coefficients of the powers of n in x.
	var_list = []
	for i in range(d+1):
		s = var("a" + str(i))
		var_list.append(s)
	
	# x is an expression at this point
	x = 0
	for i in range(len(var_list)):
		x += n**i * var_list[i]
	
	# Use method of undetermined coefficients to find explicit polynomial x that satisfies
	# the listed equation in the function header.
	R = PolynomialRing(QQ,[n])
	a = R(a)
	b = R(b)
	c = R(c)
	relations = []
	for elem in (a*x(n=n+1) - b(n-1)*x - c).full_simplify().coefficients(n):
		relations.append(elem[0] == 0)

	solved = solve(relations,var_list)
	if not solved:
		return "No polynomial solution exists"
	# x is now a polynomial over the variable n
	x = 0
	z = R(n)
 	for i in range(len(var_list)):
 		x += list(solved[i].iterator())[1]*(z^i)
	
	return x
