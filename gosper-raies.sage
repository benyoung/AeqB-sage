###############################################################################
# If given a hypergeometric term, t, the goal of Goser's algorithm is to find a
# hypergeometric term, z, such that z(n+1)-z(n)=t(n) if one exists and to say that
# there is no such term if one does not exist.
# 
# If t is a hypergeometric term in n (created via t=...n...) then gosper(t,n) will
# return the term z if one exists and will return 0 otherwise.
#
# by Dan Raies (raies@uoregon.edu)
###############################################################################

###############################################################################
# REQUIRED FILES: gosper_step_two relies on a file called resultant.sage which allows
# for the computation of resultants of polynomials.
###############################################################################

###############################################################################
# INPUT: A hypergeometric term, t, and its variable, n.  NOTE: At this point, it must
# be the case that t(n+1)/t(n) is a rational function over the field of rational numbers.
# 
# OUTPUT: The hypergeometric term, z, satisfying z(n+1)-z(n)=t(n) if one exists and 0
# otherwise.  It should be noted that z will never be zero unless t was zero to begin
# with in which case you probably don't need this algorithm to figure that out.
###############################################################################
def gosper(t,n):
    # Step 1 of Gosper's algorithm as defined in AeqB.  See gosper_step_one below.
    r = gosper_step_one(t,n)
    
    # Step 2 of Gosper's algorithm as defined in AeqB.  See gosper_step_two below.
    poly_list = gosper_step_two(r,n)
    a = poly_list[0].expand()
    b = poly_list[1].expand()
    c = poly_list[2].expand()
    
    # Step 3 of Gosper's algorithm as defined in AeqB.  See gosper_step_three below.
    X = gosper_step_three(a,b,c,n)
    
    # If X==0 was returned then the algorithm has proved that no z exists and zero
    # is returned.
    if X==0:
        return 0
    # If X is a non-zero polynomial then the appropriate value of z is returned according
    # to Step 4 of Gosper's algorithm as defined in AeqB.
    else:
        return (b.subs(n==n-1)*X*t/c).full_simplify()

###############################################################################
# INPUT: A hypergeometric term, t, and its variable, n.
# 
# OUTPUT: The rational term, r, stipulated in step 1 of Gosper's algorithm.  This simply
# returns t(n+1)/t(n).
###############################################################################
def gosper_step_one(t,n):
    return (t.subs(n==n+1)/t).full_simplify()

###############################################################################
# INPUT: A rational function, r, and its variable, n.
# 
# OUTPUT: A list of polynomials, [a,b,c], such that r(n)=(a(n)/b(n))*(c(n+1)/c(n))
# and such that a, b, and c satisfy the conditions in step 2 of Gosper's algorithm as
# specified in AeqB.
###############################################################################
def gosper_step_two(r,n):
    # resultant.sage contains a method to calculate the resultant of two polynomials.
    load resultant-raies.sage
    
    # We will require that a(n) and b(n+h) are relatively prime for all values of h.
    var('h')
    
    # r_simp is used just in case.
    r_simp = r.full_simplify()
    
    # f, g, and Z are specified on page 80 of AeqB.  NOTE: coeffs_list_descending also
    # comes from gosper.sage
    f = r_simp.numerator()/coeffs_list_descending(r_simp.numerator(),n)[0]
    g = r_simp.denominator()/coeffs_list_descending(r_simp.denominator(),n)[0]
    Z = coeffs_list_descending(r_simp.numerator(),n)[0]/coeffs_list_descending(r_simp.denominator(),n)[0]
    
    # R is the resultant of f(n) and g(n+h).  This becomes a polynomial in h.
    R = resultant(f,g(n=n+h),n)
    
    # Here we grab all nonnegative integral solutions to R(h)=0 in ascending order
    # (multiplicities are ignored).
    list_of_roots = nonnegative_integer_roots(R,h)
    
    # The following loop finds the polynomials a, b, and c according to step 2 of
    # Gosper's algorithm as specified in AeqB (pseudo-code is provided on page 80).
    p = f
    q = g
    c = 1
    for j in range(len(list_of_roots)):
        h_val = list_of_roots[j]
        # I believe that this gcd function is the reason why t(n+1)/t(n) needs to be a
        # rational function over the field of rationals; gcd(_,_) expects two rational
        # polynomials.
        s = gcd(p,q(n=n+h_val))
        p = (p/s).full_simplify()
        q = (q/s(n=n-h_val)).full_simplify()
        for k in range(1,h_val+1):
            c = c*s(n=n-k)
    a = Z*p
    b = q
    
    # I was encountering a problem while troubleshooting.  If it turns out that any of
    # a, b, and c are constants then Sage will treat them as such and it won't know what
    # to do with something like a.degree(n) because it does not see a as a polynomial in
    # n.  This is my cheat to force Sage to treat a, b, and c as polynomials in n for
    # the rest of the steps.  Let me know if there is a better way than this.
    a = a+n-n
    b = b+n-n
    c = c+n-n
    
    # Finally, we ensure that b and c are monic.
    c = c/coeffs_list_descending(c,n)[0]
    a = a/coeffs_list_descending(b,n)[0]
    b = b/coeffs_list_descending(b,n)[0]

    return [a,b,c]
    
###############################################################################
# INPUT: The polynomials a, b, and c found by gosper_step_two and their variable, n.
# 
# OUTPUT: A polynomial, X, such that a(n)X(n+1)-b(n-1)X(n)=c(n) if one exists.  The
# number 0 is returned otherwise.
###############################################################################
def gosper_step_three(a,b,c,n):
    # gosper_find_degree_X finds the degree of X if such an X exists.  If not then it
    # returns -1
    d = gosper_find_degree_X(a,b,c,n)
    
    # If d==-1 then there is no solution, X, so we return 0.
    if d==-1:
        return 0
    # If d>-1 then gosper_undetermined_coefficients uses the method of undetermined
    # coefficients to find X.  It is possible that d>-1 but there is still no solution.
    # In this case, gosper_undetermined_coefficients returns a 0.
    else:
        return gosper_undetermined_coefficients(a,b,c,d,n)
    
###############################################################################
# INPUT: The polynomials a, b, and c found by gosper_step_two and their variable, n,
# as well as the degree of the polynomial, X, for which we are searching.  This degree
# can be found with gosper_find_degree_X.  NOTE: This method assumes that there actually
# is a solution of degree d.  If there isn't then the solve() command will not work
# properly.
# 
# OUTPUT: A polynomial, X, such that a(n)X(n+1)-b(n-1)X(n)=c(n).
###############################################################################
def gosper_undetermined_coefficients(a,b,c,d,n):
    aString = "any_string_will_do"
    
    # At this point, X is a generic polynomial of degree d with arbitrary coefficients.
    X = sum([var(aString+str(i))*n^i for i in range(d+1)])
    
    # the_variables simply holds the variables used for the coefficients of X.
    the_variables = [var(aString+str(i)) for i in range(d+1)]
    
    # poly is a polynomial that we would like to be zero.  The rest of the algorithm
    # looks for coefficients of X that will satisfy poly==0.
    poly = (a*X.subs(n==n+1)-b.subs(n==n-1)*X-c).full_simplify()
    
    # We proceed to solve the system of linear equations obtained by setting the 
    # coefficients of poly equal to zero
    the_coefficients = poly.coefficients(n)
    relations = []
    for [i,j] in the_coefficients:
        relations.append(i==0);
    # the_solutions is the solution to the system.  I'm not happy about making it a
    # dictionary but this is the only way I can get it to work.
    the_solutions_list = solve(relations, the_variables, solution_dict=True)
    
    # It is possible that there are no solutions.  If that is the case, zero is returned
    if (len(the_solutions_list) == 0):
        return 0
    
    # Otherwise, we take the first one
    the_solutions = the_solutions_list[0]
    
    # We now take all of the coefficients in X and substitute the solutions that were
    # found in the_solutions
    for avar in the_variables:
        X = X.subs(avar == the_solutions[avar])
        
    # Solving a system invokes some dummy variables and we want to set all of those to
    # zero.  However, we obviously don't want all of the variables in X to be zero.
    # all_variables is all of the variables in X whereas restricted_variables should be
    # all of the variables that show up in t (or, equivalently, a, b, and c).  The loop
    # sets any variable in all_variables to zero (in X) as long as that variable isn't in
    # restricted_variables.
    all_variables = X.variables()
    restricted_variables = get_variables_list([a,b,c,n])
    for avar in all_variables:
        if not(avar in restricted_variables):
            X=X.subs(avar == 0)
    
    # X is the desired polynomial and is returned.
    return X
    
###############################################################################
# INPUT: The polynomials a, b, and c found by gosper_step_two and their variable, n.
# 
# OUTPUT: If there exists a polynomial, X, such that a(n)X(n+1)-b(n-1)X(n)=c(n) then
# the degree of that polynomial is returned.  If there does not exist such a polynomial
# then -1 is returned.
###############################################################################
def gosper_find_degree_X(a,b,c,n):
    # resultant.sage contains a method to get a list of the coefficients of a polynomial.
    load resultant-raies.sage
    
    # a_coeffs holds the coefficients of a.
    a_coeffs = coeffs_list_descending(a,n)
    
    # lc_a and lc_b are the leading coefficients of a and b, respectively.
    lc_a = a_coeffs[0]
    lc_b = coeffs_list_descending(b,n)[0]
    
    # deg_a, deg_b, and deg_c are the degrees of a, b, and c, respectively.
    deg_a = a.degree(n)
    deg_b = b.degree(n)
    deg_c = c.degree(n)
    
    # The algorithm in AeqB for finding the degree of X seems to ignore the case when
    # both a and b are constant polynomials.  In this case they have the same degree so
    # one would think that they are handled by the else: case below.  However, in that
    # case, one of the possible degrees is (B-A)/lc_a and if a and b are constant then
    # B and A are undefined.  I believe it is straightforward to prove that X must have
    # degree at-most one more than the degree of c in that case.  NOTE: In case you
    # are curious, it does come up in certain test cases that a and b are constant.  In
    # fact, if you try t=n^4 then a=b=1.
    if ((deg_a == 0) and (deg_b == 0)):
        # If we find ourselves in this case then there is no reason to come up with other
        # possible degrees.  We simply return deg_c+1
        return deg_c+1
    
    # D will hold all of the possible degrees of X.  We append any integers that could
    # possibly be the degree of X to D and return the maximum of D.  In this manner, we
    # only consider positive integers and return -1 if there are no such integers.
    D = [-1]
    # The if case is what AeqB calls Case 1 (page 85).
    if ((lc_a != lc_b) or (deg_a != deg_b)):
        D.append(deg_c - max(deg_a,deg_b))
    # The else case is what AeqB calls Case 2 (page 85).
    else:
        D.append(deg_c-deg_a+1)
        A = a_coeffs[1]
        B = coeffs_list_descending(b.subs(n==n-1),n)[1]
        possible_degree = (B-A)/lc_a
        # possible_degree should only be considered an option if it is an integer.
        if possible_degree in ZZ:
            D.append(possible_degree)
    
    return max(D)
    
###############################################################################
# INPUT: A polynomial, p, and its variable, z.
# 
# OUTPUT: A list of all positive integers which are roots of p in ascending order.
###############################################################################
def nonnegative_integer_roots(p,z):
    # This algorithm uses the rational roots theorem.  That is, if $a_i$ is an integer
    # with $a_0 \neq 0$ and $a_n \neq 0$ then the only possible rational roots of
    #       $a_n x^n + a_{n-1} x^n + \ldots + a_1 x + a_0
    # are $\pm a/b$ where $a$ divides $a_0$ and $b$ divides $a_n$.
    #
    # It follows that the only possible positive integral roots of a rational polynomial
    # are the divisors of the constant term after clearing denominators.
    
    # The first step is to remove all factors of z from p as the rational roots theorem
    # only works on a polynomial whose constant term is not zero.  q will be the
    # polynomial p/z^k for the maximal value of k (among those for which p/z^k is
    # polynomial).
    q = p
    while q.subs(z==0)==0:
        q = (q/z).full_simplify()
    
    # As q is a rational polynomial, q.numerator will be the result of clearing the
    # denominators of q.  Then q.numerator().subs(z==0) will be the constant term of
    # that polynomial.  factors_raw is then a list containing all of the factors of
    # constant_term along with their multiplicities.
    constant_term = q.numerator().subs(z==0)
    factors_raw = list(factor(Integer(constant_term)))
    
    # factors and factors_raw hold the same information.  The difference is that 
    # factors_raw[i][0] is a factor and factors_raw[i][1] is the multiplicity of that
    # factor.  factors is a 1-dimensional list that simply has multiple entries to
    # represent multiplicities.
    factors = []
    for i in range(len(factors_raw)):
        for j in range(factors_raw[i][1]):
            factors.append(factors_raw[i][0])
    
    # subsets_of_factors contains all sub-multisets of the list of factors.  To obtain
    # all divisors of constant_term it then suffices to take all of the sets in
    # subsets_of_factors and multiply them together.
    subsets_of_factors = list(Subsets(factors,submultiset=True))
    
    # list_of_roots will hold all of the nonnegative integer roots of p.
    list_of_roots = []
    
    # 0 is added to list_of_roots if needed.
    if p.subs(z==0)==0:
        list_of_roots.append(0)
    
    # We loop through all of the subsets_of_factors.  At each stage, test_val is
    # constructed as the product of all of the elements in subsets_of_factors[i].  Then,
    # if p(test_val)=0 we add test_val to list_of_roots.
    for i in range(len(subsets_of_factors)):
        test_val = 1
        for j in range(len(subsets_of_factors[i])):
            test_val = test_val * subsets_of_factors[i][j]
        if p.subs(z==test_val)==0:
            list_of_roots.append(test_val)
    
    # Finally, a sorted version of list_of_roots is returned.
    return sorted(list_of_roots)
    
###############################################################################
# INPUT: A list of expressions.
# 
# OUTPUT: A set of all variables that appear in all of those expressions.
###############################################################################
def get_variables_list(list):
    newList = []
    for item in list:
        vars = item.variables()
        for v in vars:
            newList.append(v)
    return Set(newList)
    
###############################################################################
# INPUT: Expressions z and t along with their (shared) variable, n.
# 
# OUTPUT: The (simplified) value of z(n+1)-z(n)-t(n).  If z was obtained from gosper(t,n)
# and z is non-zero then it should be that gosper_test_sum(z,t,n) returns 0.
###############################################################################
def gosper_test_sum(z,t,n):
    return (z.subs(n==n+1)-z-t).full_simplify()

###############################################################################
# INPUT: A hypergeometric term, t, and its variable, n.  The input should be the same
# as that of gosper(t,n).
# 
# OUTPUT: A sentence describing the result of gosper(t,n).  I do not know what will
# happen if t is not hypergeometric.
###############################################################################
def gosper_test(t,n):
    z = gosper(t,n)
    if z==0:
        return "It was found that the term is not Gosper summable"
    else:
        k = gosper_test_sum(z,t,n)
        if k==0:
            return "A hypergeometric term was found satisfying t(n)=z(n+1)-z(n)"
        else:
            return "A hypergeometric term was found that does not satisfy t(n)=z(n+1)-z(n)"

###############################################################################
# INPUT: A list of expressions and a single variable, n.
# 
# OUTPUT: A univariate polynomial ring in n over the multivariate polynomial ring whose
# variables are those that show up in list (other than n).  If n is the only variable
# in list then real polynomials 
###############################################################################
def get_ring(list,n):
    # all of the variables that show up in the list
    all_variables = get_variables_list(list)
    
    expression_variables = []
    for element in all_variables:
        if element != n:
            expression_variables.append(element)
    if len(expression_variables) == 0:
        return PolynomialRing(RationalField(),n)
    else:
        return PolynomialRing(PolynomialRing(RationalField(),expression_variables),n)