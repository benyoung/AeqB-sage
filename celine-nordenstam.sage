#Implementation of Sister Celine's algorithm
#If successful it returns coefficients aibj such that
#\sum_{i=0}^I \sum_{j=0}^J   aibj(v1) * f(v1+i, v2+i)  ==  0
#
#Example: 
#   n,k = var('n,k')
#   celine( k * binomial(n,k), n, k, 1, 1) 
#
def celine(f, v1, v2, I, J):
    acc = 0 
    variables = []
    for i in range(I+1):
        for j in range(J+1):
            x = var('a' + str(i) + 'b' + str(j))
            variables.append(x)
            acc += x * (f.substitute({v1:v1-i, v2:v2-j})/f).full_simplify()
    target = acc.factor().numerator()
    relations = []
    for [x,y] in target.coefficients(k):
        relations.append(x==0)
    return solve(relations, variables)

