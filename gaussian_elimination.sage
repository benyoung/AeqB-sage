def gaussian_elimination(R,N):
    # Convert input matrix into a list with entries confined to input ring. If a conversion
    # error occurs, the program will raise an error and terminate.
    M = [[R(0) for i in range(N.ncols())] for j in range(N.nrows())]
    for i in range(N.nrows()):
    	for j in range(N.ncols()):
    		try:
    			M[i][j] = R(N[i][j])
    		except (TypeError):
    			raise AssertionError("Input matrix should have integer entries")
    
    m_rows = len(M)
    m_cols = len(M[0])
    
    # Go from top row to bottom row.
    for i in range(m_rows):
        pivot = -1
        for j in range(m_cols):
            if M[i][j] != 0:
                pivot = j
                break
        if pivot == -1:
            continue
        constant = M[i][pivot]
        for j in range(pivot, m_cols):
            M[i][j] /= (constant)
        for j in range(i + 1, m_rows):
            if M[j][pivot] != 0:
                constant2 = (M[j][pivot])
                for k in range(pivot,m_cols):
                    M[j][k] -= (M[i][k]*constant2)

    # Go from bottom row to top row.
    for i in range(m_rows - 1, -1, -1):
        pivot = -1
        for j in range(m_cols):
            if M[i][j] != 0:
                pivot = j
                break
        if pivot != -1:
            constant = M[i][pivot]
            for j in range(pivot,m_cols):
                M[i][j] /= (constant)
            for j in range(i - 1, -1, -1):
                if M[j][pivot] != 0:
                    constant2 = (M[j][pivot])
                    for k in range(pivot,m_cols):
                        M[j][k] -= (M[i][k]*constant2)
    
    return matrix(M)
