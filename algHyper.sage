def algHyper(polyList,n):
    if polyList[0].degree(n) == 1:
        aList = [polyList[0]/polyList[0].coefficient(n),1]
    else:
        aList = list(polyList[0].factor().iterator())
        for i in range(len(aList)):
	    aList[i] /= aList[i].expand().coefficients()[-1][0]
	if 1*n^0 not in aList:
	    aList.append(1*n^0)

    d = len(polyList) - 1
    if polyList[d].substitute(n==n-d+1).degree(n) == 1:
        bList = [polyList[d].substitute(n==n-d+1)/polyList[d].substitute(n==n-d+1).coefficient(n),1]
    else:
        bList = list(polyList[d].substitute(n==n-d+1).factor().iterator())
        for i in range(len(bList)):
	    bList[i] /= bList[i].expand().coefficients()[-1][0]
	if 1*n^0 not in aList:
	    bList.append(1*n^0)
    
    print aList, bList
    
    for aElem in aList:
        for bElem in bList:
	    PList = []
	    for i in range(d+1):
		p_i = polyList[i]
		for j in range(i):
		    p_i *= aElem.substitute(n==n+j)
		for j in range(i,d):
		    p_i *= bElem.substitute(n==n+j)
		PList.append(p_i)
	    
	    m = max([PList[i].degree(n) for i in range(len(PList))])
	    alpha = [PList[i].expand().coefficient(n**m) for i in range(len(PList))]
	    var('z')
	    zPol = 0*z
	    for i in range(len(alpha)):
		zPol += alpha[i]*z**i

	    vals = []
	    try:
	        vals = [elem[0] for elem in zPol.roots() if elem[0] != 0]
	    except:
	        pass
	    
	    for x in vals:
	        polyList2 = [x**i*PList[i] for i in range(len(PList))]
	        relations = algPoly(polyList2,0*n,n)[0]
	        coeff = [list(relations[i].iterator())[1] for i in range(len(relations))]
	 	if coeff != [0*n for i in range(len(coeff))]: 
	 	    c = 0*n
	 	    for i in range(len(coeff)):
		        c += n**i * coeff[i]
		    S = x*aElem/bElem*c.substitute(n==n+1)/c
		    return S
    return 0

