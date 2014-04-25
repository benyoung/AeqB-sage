def creative_telescoping(f,J,n,k):
    # Takes as input an expression identifier f, length of recurrence being checked J,
    # and variables n and k which identifier f is defined over. Returns the coefficients
    # ai of the recurrence relation sum(ai*f(n+i,k),i,0,J) = g(n,k+1)-g(n,k). Follows the
    # explanation of creative telescoping described in "A=B".
    #
    # Ex: F(n,k) = binomial(n,k)^2
    # creative_telescoping(F,1,n,k) ===> [a0 == -2*(2*n+1), a1 == n+1]
    #
    
    # (Comment to editors: one could return the coefficients of the polynomial x that
    # defines g(n,k) in terms of f(n,k) if that information is desired instead of the
    # coefficients of the recurrence relation. Also, this method hasn't been extensively
    # tested, but working on it.)
    
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
    
    # Factorize q(n,k) in the same fashion as that in Gosper's algorithm
    q(n,k) = (r(n,k)/s(n,k)).full_simplify()
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

def recurse_solve(a,b,c,n,vars):
	b(n) = b
	a(n) = a
	c(n) = c
	deg_a = a.degree(n)
	deg_b = b.degree(n)
	deg_c = c.degree(n)
	lc_a = a.coefficient(a.variables()[0]^deg_a)
	lc_b = b.coefficient(b.variables()[0]^deg_b)
	
	# Determine degree of polynomial x that will partially define g(n,k)
	if (deg_a != deg_b) or lc_a != lc_b:
		d = deg_c - max(deg_a,deg_b)
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
 	

    