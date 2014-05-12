def algPoly(polyList,f,n):
	qList = []
	for j in range(len(polyList)):
	   	qj = 0*n
	    	for i in range(j,len(polyList)):
			qj += binomial(i,j)*polyList[i]
	    	qList.append(qj)
	
	b = max([(qList[j].degree(n) - j) for j in range(len(qList))])
	first = f.degree(n) - b
	second = -1*b - 1
	
	lc = [qList[j].coefficients()[-1][0] for j in range(len(qList))]
	alpha = 0*n
	for j in range(len(qList)):
		if (qList[j].degree(n) - j) == b:
			alpha += lc[j]*fallingFactorial(n,j)
	
	deg = alpha.roots()
	third = max([deg[j][0] for j in range(len(deg))])
	d = max(first,second,third)
    	
	varList = [var('a' + str(j)) for j in range(d+1)]
    	pol = 0*n
	for i in range(d+1):
		pol += n**i * varList[i]
	
	solution = -1*f
	for i in range(len(polyList)):
	    	solution += pol.substitute(n==n+i)*polyList[i]
	print solution.coefficients(n)
	coeff = solution.coefficients(n)
	solved = []
	for i in range(len(coeff)):
	    	solved.append(coeff[i][0] == 0)
	return solve(solved,varList)
	

def fallingFactorial(x,j):
	result = 1
    	for i in range(j):
		result *= (x-i)
	return result
