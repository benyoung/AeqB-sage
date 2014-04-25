###############################################################################
# The goal here is to find a recurrence relation for a proper hypergeometric term using
# Sister Celine's algorithm.  It is assumed that F is created in the following manner:
# def F(n,k):
#   return ...
# 
# We look for a recurrence of the form
#       $\sum_{i=0}^{I} \sum_{j=0}^{J} a_{i,j} F(n+i,k+j) = 0$
# where $a_{i,j}$ is a rational expression in n.  Note: This notation is slightly 
# different than (4.3.1) in AeqB.
# 
# celine(F,n,k) will output a list of three things.  The second is I, the third is J, and
# the first is a dictionary.  If this dictionary is assigned to 'rec' in Sage - perhaps 
# by a command like rec = celine(F,n,k)[0] - then rec[(i,j)] will be the coefficient
# $a_{i,j}$ in the recurrence.
#
# by Dan Raies (raies@uoregon.edu)
###############################################################################

###############################################################################
# INPUT: A function F of two variables and its two variables.  It is assumed to be 
# defined as follows:
# def F(n,k):
#   return ...
# F should be a proper hypergeometric term.
#
# OUTPUT: A list of the form [rec, I, J].  I and J describe the depth of the recursion
# and rec[(i,j)] is the appropriate coefficient from Sister Celine's algorithm
###############################################################################

def celine(F,n,k):
    # Q is a function for convenience later.
    def Q(a,b):
        return (F(n + a, k + b) / F(n, k)).full_simplify()
    
    # There is a point later where we want to set all dummy variables equal to 1.
    # This is done by setting ALL variables equal to 1 that aren't in this list of
    # reserved_variables
    reserved_variables = Set(F(n,k).variables())
    
    # The first thing to do is find the required level of recursion.  That is the goal of
    # the while loop below.  I and J will be incremented at each step and when a suitable
    # depth is found, found_solution will be changed to true and the loop will stop.
    found_solution = false
    I = 0
    J = 0
    rational_numerator = 0
    while not(found_solution):
        I = I + 1
        J = J + 1
        Inds = CartesianProduct(range(I+1), range(J+1))
        
        # rational_recurrence is the rational expression that looks like
        #       a00*Q(0,0)+a01*Q(0,1)+a10*Q(1,0)+a11*Q(1,1)+...+aIJ*Q(I,J).
        rational_recurrence = sum([var("a"+str(i)+str(j))*Q(i,j) for i,j in Inds])
        
        # We then find a common denominator in rational_recurrence and assign the
        # numerator to rational_numerator.
        rational_numerator = rational_recurrence.factor().numerator()
        
        # rational_numerator is a polynomial in k.  Soon we will obtain a system of
        # equations by setting the k-coefficients of this polynomial equal to 0.  This
        # system will be solved for each of 'aij'.  Hence the number of equations will
        # be the number of non-zero k-coefficients and the number of variables will 
        # be the number of aij, or (I+1)*(J+1) 
        number_of_equations =  len(rational_numerator.coefficients(k))
        number_of_variables = (I+1)*(J+1)
        
        # If the number of variables is strictly larger than the number of the equations
        # then there is necessarily a non-trivial solution so we have found the proper
        # depth of recursion.  In that case we set found_solution to true and stop
        # looping.
        if (number_of_equations < number_of_variables):
            found_solution = true
    
    # the_variables holds the list of variables for which we are solving - that is, the
    # coefficients aij
    the_variables = [var("a"+str(i)+str(j)) for i,j in Inds]
    
    # We get equations by setting the k-coefficients of rational_numerator to zero.  Here,
    # the_coefficients holds these coefficients
    the_coefficients = rational_numerator.coefficients(k)
    
    # To solve the system, we set the expressions in the_coefficients equal to zero and
    # solve for the variables in the_variables.
    relations = []
    for [x,y] in the_coefficients:
        relations.append(x==0);
    the_solutions = solve(relations, the_variables)[0]
    
    ###########################################################################
    # NOTE: at this point, the system is solved and solutions are held in the_solutions.
    # The rest of the code puts these solutions into a dictionary
    ###########################################################################
    
    # recurrence will be a dictionary which holds the aij.  Here it is initialized so
    # that recurrence[(i,j)] is zero by default.
    from collections import defaultdict
    recurrence = defaultdict(int)
    
    # It is likely that the solutions contain dummy variables.  The first thing to do is
    # collect those dummy variables into a single set.  These are held in set_of_dummies.
    # A set is used instead of a list to avoid repetition
    dummies = []
    for soln in the_solutions:
        vars = soln.rhs().variables()
        # We want dummy variables to be all of the variables that show up in our solutions
        # except for those that were reserved in reserved_variables.
        for v in vars:
            if not(v in reserved_variables):
                dummies.append(v)
    set_of_dummies = Set(dummies)
    
    # The solutions in the_solutions are all of the form
    #       aij == ...some function of n and possibly other variables...
    # We wish to extract the right hand sides of these expressions, set all the dummy
    # variables equal to zero, and then store them in the appropriate place in the
    # dictionary.  Note that, by construction, the solutions in the_solutions are ordered
    # in whatever order Sage uses when it calls "for i,j in Inds".  I have no idea what
    # that order is, but I know how to call "for i,j in Inds".  Hence I turn the_solutions
    # into something over which I can iterate and then loop over Inds.  At each step the
    # .next() thing in the list will be the appropriate solution
    list_of_solutions = iter(the_solutions)
    for i,j in Inds:
        soln = list_of_solutions.next().rhs()
        for v in set_of_dummies:
            soln = soln.subs(v==1)
        recurrence[(i,j)] = soln
    
    return [recurrence,I,J]