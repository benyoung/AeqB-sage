##########################################################################################
# DESCRIPTION: Follows the explanation of Algorithm Poly as described in "A=B" by  
# Marko Petkovsek, Herbert Wilf, and Doron Zeilberger.
#
# INPUT: Takes as input a polynomial list that describes the recurrence operator
# L, the polynomial that the recurrence satisfies (if inhomogeneous), and the
# variable in which we view the polynomials as being functions of.
#
# OUTPUT: Returns the coefficients of the polynomial that satisfies the recurrence
# relation.
#
# AUTHOR: Kevin Wilson, kwilson8@uoregon.edu
##########################################################################################

def algPoly(polyList,f,n):
	# Construct polynomials qj that arise from treating the shift operator
    	# as the difference operator plus the identity operator.
    	qList = []
	for j in range(len(polyList)):
	   	qj = 0*n
	    	for i in range(j,len(polyList)):
			qj += binomial(i,j)*polyList[i]
	    	qList.append(qj)
	
	# Find all candidates for the degree bound on output polynomial.
	b = max([(qList[j].degree(n) - j) for j in range(len(qList))])
	first = f.degree(n) - b
	second = -1*b - 1
	lc = [qList[j].expand().coefficients()[-1][0] for j in range(len(qList))]
	
	alpha = 0*n
	for j in range(len(qList)):
		if (qList[j].degree(n) - j) == b:
			alpha += lc[j]*fallingFactorial(n,j)
	try:
    	    deg = alpha.roots()
	    third = max([deg[j][0] for j in range(len(deg))])
	    d = max(first,second,third,0)
        except:
	    d = max(first,second,0)
        
    	# Use method of undetermined coefficients to find the output polynomial.
	varList = [var('a' + str(j)) for j in range(d+1)]
    	pol = 0*n
	for i in range(d+1):
		pol += n**i * varList[i]
	solution = -1*f
	for i in range(len(polyList)):
	    	solution += pol.substitute(n==n+i)*polyList[i]
	coeff = solution.coefficients(n)
	solved = []
	for i in range(len(coeff)):
	    	solved.append(coeff[i][0] == 0)
	return solve(solved,varList)
	
##########################################################################################
# DESCRIPTION: Finds the falling factorial of the input value with j
# multiplications.
#
# INPUT: Takes as input expression or value x, and the number of terms of
# multiplication we desire in our falling factorial.
#
# OUTPUT: Returns the desired falling factorial.
##########################################################################################

def fallingFactorial(x,j):
	result = 1
    	for i in range(j):
		result *= (x-i)
	return result
