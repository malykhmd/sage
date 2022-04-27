#Calculation of automorpisms of a hyperelliptic curve over CC.
#Package for sagemath, ver. 1.2, 13/6/2019, Malykh M.D.
def cross_ratio(a, b, c, d): 
    if a==b==c or a==c==d or b==c==d:
        ans='undef'
    elif a==c or b==d:
        ans=0
    elif b==c or a==d:
        ans=oo
    elif a==b or c==d:
        ans=1
    elif a==oo: 
        ans=(b - d)/(b - c)
    elif b==oo: 
        ans=(a - c)/(a - d)
    elif c==oo: 
        ans=(b - d)/(a - d)
    elif d==oo: 
        ans=(a - c)/(b - c)
    else:
        ans=((a-c)/(c-b))/((a-d)/(d-b))
    return ans
# This func. return the solution of the equation (a,b,c,xx)=A(x).
def cross_ratio_solution(a, b, c, A): 
    if a==b or a==c or b==c or (A in QQbar):
        ans='error: a==b or a==c or b==c'
    elif a==oo: 
        ans=-A*(b - c) + b
    elif b==oo: 
        ans=(A*a - a + c)/A
    elif c==oo: 
        ans=(A*a - b)/(A - 1)
    else:
        ans=((a*b - a*c)*A - a*b + b*c)/(A*(b - c) - a + c)
    return ans
def list_of_roots(p):
    KK=QQbar[x]
    L=KK(p).roots()
    if prod([m for [alpha,m] in L]) ==1:
        ans=[alpha for [alpha,m] in L]
        if is_odd(KK(p).degree()):
            ans.append(oo)
    else:
        ans='error: there are multiple roots'
    return ans
def list_of_permutations(p,L):
    ans=[m for m in Permutations(range(1,len(L)+1))[1:] if cross_ratio_test(L,m)==1]
    return ans
def cross_ratio_test(L,m):
# If 4 poits don't change then Moebius substitution is Id. 
    if sum(L[n]==L[m[n]-1] for n in range(len(L)))>3:
        ans=False
    else: 
        ans=True
        for n in range(3,len(L)):
            if cross_ratio(L[0],L[1],L[2],L[n])!=cross_ratio(L[m[0]-1],L[m[1]-1],L[m[2]-1],L[m[n]-1]):
                ans=False
                break 
    return ans
def list_of_substitutions(p):
    L=list_of_roots(p)
    if L=='error: there are multiple roots':
        ans='error: p(x)=0 has multiple roots'
    else:
        n=QQbar[x](p).degree()
        M=list_of_permutations(p,L)
        ans=[[xx==QQbar[x](x), yy==QQbar[x,y](y)], [xx==QQbar[x](x), yy==-QQbar[x,y](y)]]
        F=(1/QQbar[x](x)).parent()
        for m in M:
            X=F(cross_ratio_solution(L[m[0]-1],L[m[1]-1],L[m[2]-1],cross_ratio(L[0],L[1],L[2],QQbar[x](x))))
            pp=p.subs(x==X)
            Y=sqrt(F(pp/p))*y
            ans+=[[xx==X, yy==Y], [xx==X, yy==-Y]]
    return ans
def pi_group_of_automorphisms(p):
    L=list_of_roots(p)
    if L=='error: there are multiple roots':
        ans='error: p(x)=0 has multiple roots'
    else:
        S=list_of_permutations(p,L)
        ans=PermutationGroup([Permutation(s).to_cycles() for s in S])
    return ans
