###############################################################################
# The following is a collection of functions whose goal is to compute the resultant
# of two polynomials.  In short, if p and q are polynomials in x then resultant(p,q,x)
# will yield the resultant in Sage.
# 
# The current implementation relies on the computing the determinant of the Sylvester
# matrix.  This is probably not the most efficient way.
# 
# by Dan Raies (raies@uoregon.edu)
###############################################################################

###############################################################################
# INPUT: A univariate polynomial, p, and its variable, x.  It should be the case that
# p is some polynomial in R[x] for some ring, R.
# 
# OUTPUT: A list (object) containing the coefficients of p in order of descending powers
# of x.  For example, if p=3x^4+6x^2-2x+1 then coeffs_list_descending(p,x) should output
# [3, 0, 6, -2, 1]
###############################################################################
def coeffs_list_descending(p,x):
    # p_degree holds the degree of p so that the list should have a length of
    # (p_degree + 1).
    p_degree = p.degree(x)
    
    # p.coefficients(x) holds the x-coefficients of p along with some other stuff.  From
    # here on out, we pull out the coefficients and put them in p_coeffs (in the desired
    # order).
    p_coeffs_raw = p.coefficients(x)
    p_coeffs = []
    
    # Loop over the powers of p
    for i in range(p_degree+1):
        # At this point in the loop we will be looking for the coefficient of the
        # x^(this_power) term in p.  If we look for the coefficient of the x^i term here
        # then the list will be out of order at the end.
        this_power = p_degree - i
        
        # If the desired coefficient is actually 0 then it will have no entry in
        # p_coeffs_raw but we still want a zero to show up in p_coeffs.  If we actually
        # find a coefficient in p_coeffs_raw then we change found to true.  At the end,
        # if found is false then we put a 0 in p_coeffs
        found = false
        
        # We loop through all of the entries of p_coeffs_raw looking for a coefficient
        # in degree this_power.
        for j in range(len(p_coeffs_raw)):
            # This checks if the j-th entry of p_coeffs_raw is of degree this_power.
            # If it is, we add the coefficient to p_coeffs and change found to true.
            if this_power == p_coeffs_raw[j][1]:
                found = true
                p_coeffs.append(p_coeffs_raw[j][0])
        
        # If, after looping through all of the coefficients in p_coeffs_raw, we found no
        # coefficient in degree this_power then the coefficient must have been zero so
        # we add a zero to p_coeffs
        if not(found):
            p_coeffs.append(0)
    
    return p_coeffs

###############################################################################
# INPUT: Two univariate polynomials, p and q, and their variable.  It should be the case
# that p and q both belong to R[x] for some ring, R.
# 
# OUTPUT: The Sylvester matrix for those two polynomials.
###############################################################################
def sylvester_matrix(p,q,x):
    
    # We first fetch the coefficients of p and q.
    p_coeffs = coeffs_list_descending(p,x)
    q_coeffs = coeffs_list_descending(q,x)
    
    # The first N rows will contain shifts of the coefficients of p and
    # the next M rows will contain shifts of the coefficients of q.
    M = len(p_coeffs)-1
    N = len(q_coeffs)-1
    
    # list will contain a list of rows.  Matrix(list) will be returned
    list = []
    
    # The first loop adds the first N rows.
    for row_num in range(N):
        # this_row will be the row that is added.
        this_row = []
        
        # We first add the appropriate number of zeros to the row.
        for i in range(row_num):
            this_row.append(0)
        
        # Next, we add the coefficients of p to the row.
        for i in range(M+1):
            this_row.append(p_coeffs[i])
        
        # Finally, we fill in the row with zeros.
        for i in range(N-row_num-1):
            this_row.append(0)
        
        # this_row is added to list
        list.append(this_row)
    
    # The second loop adds the next M rows
    for row_num in range(M):
        # this_row will be the row that is added.
        this_row=[]
        
        # We first add the appropriate number of zeros to the row.
        for i in range(row_num):
            this_row.append(0)
        
        # Next, we add the coefficients of p to the row.
        for i in range(N+1):
            this_row.append(q_coeffs[i])
        
        # Finally, we fill in the row with zeros.
        for i in range(M-row_num-1):
            this_row.append(0)
        
        # this_row is added to list
        list.append(this_row)
    
    return Matrix(list)

###############################################################################
# INPUT: Two univariate polynomials, p and q, and their variable.  It should be the case
# that p and q both belong to R[x] for some ring, R.
# 
# OUTPUT: The resultant of these two polynomials (in the ring, R).
###############################################################################
def resultant(p,q,x):
    return sylvester_matrix(p,q,x).determinant()