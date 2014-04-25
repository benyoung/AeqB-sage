##########################################################################################
# DESCRIPTION: Follows the explanation of creative telescoping described in "A=B" by 
# Marko Petkovsek, Herbert Wilf, and Doron Zeilberger.
#
# INPUT: Takes as input an expression identifier f, length of recurrence being checked J,       
# and variables n and k which identifier f is defined over. 
#
# OUTPUT: Returns the coefficients, ai, of the recurrence relation sum(ai*f(n+i,k),i,0,J) = g(n,k+1)-g(n,k). 
#
# Ex: F(n,k) = binomial(n,k)^2; creative_telescoping(F,1,n,k) ===> [a0 == -2*(2*n+1), a1 == n+1]
#
# AUTHOR: Kevin Wilson, kwilson8@uoregon.edu
#
# (Comment to editors: one could return the coefficients of the polynomial x that
# defines g(n,k) in terms of f(n,k) if that information is desired instead of the
# coefficients of the recurrence relation. Also, this method hasn't been extensively
# tested, but working on it.)
##########################################################################################

def creative_telescoping(f,J,n,k):    
    # Create to-be-determined coefficients of recurrence relation we find a0,...,aJ
    var_list = []
    for i in range(J+1):
        var_list.append("a" + str(i))
    vars = list(var(var_list))
    
    t(n,k) = 0
    for i in range(J+1):
        t += f(n+i,k)*vars[i]
    r1(n,k) = (f(n,k+1)/f(n,k)).full_simplify().numerator()
    r2(n,k) = (f(n,k+1)/f(n,k)).full_simplify().denominator()
    s1(n,k) = (f(n,k)/f(n-1,k)).full_simplify().numerator()
    s2(n,k) = (f(n,k)/f(n-1,k)).full_simplify().denominator()
    
    # Create p_0(n,k), r(n,k), and s(n,k) such that p_0(n,k+1)/p_0(n,k)*r(n,k)/s(n,k) = t(n,k+1)/t(n,k)
    p_0(n,k) = 0
    for j in range(J+1):
        prod1 = 1
        prod2 = 1
        for i in range(j):
            prod1 *= s1(n+j-i,k)
        for i in range(j+1,J+1):
            prod2 *= s2(n+i,k)
        p_0 += vars[j]*prod1*prod2
    
    r(n,k) = r1(n,k)
    for i in range(1,J+1):
        r *= s2(n+i,k)
    s(n,k) = r2(n,k)
    for i in range(1,J+1):
        s *= s2(n+i,k+1)
    q(n,k) = (r(n,k)/s(n,k)).full_simplify()
    
    # Factorize q(n,k) to get polynomials p1, p2, and p3 such that the following holds:
    # p1(n,k+1)/p1(n,k)*p2(n,k)/p3(n,k) = q(n,k)
    factorized = gosper_step2(q,n,k)
    p1(n,k) = (factorized[2])
    p2(n,k) = (factorized[0])
    p3(n,k) = (factorized[1])
    p(n,k) = p1(n,k)*p_0(n,k)
    x = recurse_solve(p2,p3,p,k,vars)
    if x == "No polynomial solution exists":
        return "Could not find recurrence for this value of J"
    else:
        return x

##########################################################################################
# DESCRIPTION: Find the polynomial x(n) such that a(n)x(n+1) - b(n-1)x(n) = c(n), and this
# solving simultaneously solves for the values of a0,...,aJ such that the following holds:
# sum(ai*f(n+i,k),i,0,J). This step of creative telescoping is described in "A=B" by
# Marko Petkovsek, Herbert Wilf, and Doron Zeilberger.
#
# INPUT: Takes as input expressions a, b, and c defined in variables n and k, variable n, 
# and variable list that contains a0,...,aJ.
#
# OUTPUT: Returns the coefficients a0,...,aJ of the recurrence relation sum(ai*f(n+i,k),i,0,J) = g(n,k+1)-g(n,k). 
#
# AUTHOR: Kevin Wilson, kwilson8@uoregon.edu
##########################################################################################

def recurse_solve(a,b,c,n,vars):
	# Find the degrees of expressions a, b, and c as well as the leading coefficients of
	# such expressions.
	a(n) = a
	b(n) = b
	c(n) = c
	deg_a = a.degree(n)
	deg_b = b.degree(n)
	deg_c = c.degree(n)
	lc_a = a.coefficient(a.variables()[0]^deg_a)
	lc_b = b.coefficient(b.variables()[0]^deg_b)
	
	# Determine degree of polynomial x(n).
	# If degree of a doesn't equal degree of b or the the degrees are equal but the leading
	# coefficients aren't, then the degree of x is the degree of c minus the maximum of
	# the degrees of a and b since a(n)*x(n+1) - b(n-1)*x(n) = c(n).
	if (deg_a != deg_b) or lc_a != lc_b:
		d = deg_c - max(deg_a,deg_b)
	# Otherwise, find the next highest power with unmatching coefficients, and set degree
	# of x to be the degree of c minus this power.
	else:
	    i = deg_a - 1
	    A = a.coefficient(a.variables()[0]^(i))
	    B = b.coefficient(b.variables()[0]^(i))
	    while (A == B):
	        i -= 1
	        A = a.coefficient(a.variables()[0]^(i))
	        B = b.coefficient(b.variables()[0]^(i))
	    d = deg_c - i
	if d < 0:
	    return "No polynomial solution exists"
	
	# Create coefficients of x.	
	var_list = []
	for i in range(d+1):
		s = var("b" + str(i))
		var_list.append(s)

	x(n) = 0
	for i in range(len(var_list)):
		x += n**i * var_list[i]
	
	# Using method of undetermined coefficients, solve for the polynomial x and the values
	# of the variables passed into the function in the vars parameter    
	relations = []
	var_list2 = vars + var_list
 	for elem in (a(n)*x(n+1) - b(n-1)*x(n) - c(n)).full_simplify().coefficients(n):
 		relations.append(elem[0] == 0)
 	intermed_solved = solve(relations,var_list2)[0]
 	if not intermed_solved:
 	    return "No polynomial solution exists"
 	
	# We only return the values of the variables passed into the function in the vars parameter
 	solved = []
 	for i in range(len(vars)):
 	    solved.append(intermed_solved[i])
 	return solved
 	

    