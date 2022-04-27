#Lagutinski package, v. 1.6
#Добавлен mu
def lagutinski_matrix(R,D,B,N):
    phi=B[:N]
    L=[phi]
    for n in range(1,N): 
        phi = [D(phi[i]) for i in range(N)]	
        L.append(phi)
    return matrix(L)
def lagutinski_det_random(R,D,B,N,p=0):
    L={R.gens()[i]:floor(random()*100) for i in range(len(R.gens()))}
    if p==0:
        L2=det(lagutinski_matrix(R,D,B,N).subs(L))
    else:
        A=lagutinski_matrix(R,D,B,N).subs(L)
        L2=matrix(ZZ.quo(p), N, lambda n,m: A[n][m]).det()
    return L2
def lagutinski_gcd_random(R,D,B,N):
    L={R.gens()[i]:floor(random()*100) for i in range(len(R.gens()))}
    L2=det(lagutinski_matrix(R,D,B,N).subs(L))
    L3=det(lagutinski_matrix(R,D,B,N+1).subs(L))
    L4=det(lagutinski_matrix(R,D,B,N+2).subs(L))
    return gcd([int(L2),int(L3),int(L3)])
def lagutinski_det(R,D,B,N):
    phi=B[:N]
    L=[phi]
    for n in range(1,N): 
        phi = [D(phi[i]) for i in range(N)]	
        L.append(phi)
    return det(matrix(L))
def lagutinski_integral(R,D,B,N):
    phi=B[:N]
    L=[phi]
    for n in range(1,N-1): 
       phi = [D(phi[i]) for i in range(N)]	
       L.append(phi)
    A = matrix(L)
    psi2 = det(A.delete_columns([N-1]))
    if psi2==0:
        r='Reduce N'
    else:        
        n=0
        psi1 = det(A.delete_columns([n]))
# For the case when psi1/psi2 is a constant. 
        while psi1/psi2 in CC:
             n = n+1
             psi1 = det(A.delete_columns([n])) 
        r = psi1/psi2
    return r
def lagutinski_integrals(R,D,B,N):
    phi=B[:N]
    L=[phi]
    for n in range(1,N-1): 
       phi = [D(phi[i]) for i in range(N)]	
       L.append(phi)
    A = matrix(L)
    psi2 = det(A.delete_columns([N-1]))
    if psi2==0:
        ans='Reduce N'
    else:        
        ans=[det(A.delete_columns([n]))/psi2 for n in range(N-2)]
    return ans
def lagutinski_ab(p,q):
	s=solve([p==0,q==0],[x,y])[0]
	x0=x.subs(s[0])
	y0=y.subs(s[1])
	a11=diff(p,x).subs(x=x0).subs(y=y0)
	a12=diff(p,y).subs(x=x0).subs(y=y0)
	a21=diff(q,x).subs(x=x0).subs(y=y0)
	a22=diff(q,y).subs(x=x0).subs(y=y0)
	P=lambda x,y: a11*(x-x0)+a12*(y-y0)
	Q=lambda x,y: a21*(x-x0)+a22*(y-y0)
	var('xi, eta, dxi, deta, alpha')
	d= lambda f: diff(f,xi)*dxi+diff(f,eta)*deta
	xx=x0+xi+alpha*eta
	yy=y0+eta
	omega=P(xx,yy)*d(xx)+Q(xx,yy)*d(yy)
	a=solve(diff(diff(omega,deta),eta)==0,alpha)[0]
	if 	diff(omega.subs(a),xi,deta)==0: 
		t=diff(omega.subs(a),eta,dxi) in QQ
	else:
		t=diff(omega.subs(a),eta,dxi)/diff(omega.subs(a),xi,deta) in QQ
	return t
def lagutinski_micronomial(R,D,B,N,M):
    if M==3:
        B3=[[n,m,k] for n in B[:N] for m in B[:N] for k in B[:N] if (n<m<k)]
    elif M==4:
        B3=[[n,m,k,l] for n in B[:N] for m in B[:N] for k in B[:N] for l in B[:N] if (n<m<k<l)]
    elif M==5:
        B3=[[n,m,k,l,ll] for n in B[:N] for m in B[:N] for k in B[:N] for l in B[:N] for ll in B[:N] if (n<m<k<l<ll)]
# Распараллеливание входа к lagutinski_det.
    @parallel
    def lagutinski_micronomial_N(n):
        return (lagutinski_det(R,D,B3[n],len(B3[n]))==0)
# Рандомизация не приводит к заметному ускорению 
    N3=[kk for kk in list(lagutinski_micronomial_N(range(len(B3)))) if kk[1]]
    if N3==[]:
        ans=[]
    else:
        n=N3[0][0][0][0]
        ans=lagutinski_integral(R,D,B3[n],M)
    return ans
def lagutinski_uv(R,p,q,N):
    B = sorted (((1 + x + y )^(2*N) + v *(1 + x + y )^(2*N)). monomials(), reverse =0)
    r=v*(diff(p,y)*q-diff(q,y)*p) - (diff(p,y) - diff(q,x))*diff(q,y)+q*diff((diff(p,y) - diff(q,x)),y)
#    r= v*q^2*diff(p/q,y) + q^2 * diff(1/q *(diff(p,y) - diff(q,x)),y)
    D = lambda F: q^2*diff(F,x)-p*q*diff(F,y) + r*diff(F,v)
    L=lagutinski_det(R,D,B,N)
    if L==0:
        L0=list(lagutinski_integral(R,D,B,N).factor())+list((lagutinski_integral(R,D,B,N)-1).factor())
# Что если $v$ будет в некю степени?
    else: 
        L0=list(L.factor())
    L1=[a for (a,b) in L0 if a.degree(v)==1] 
    L2=[f for f in L1 if D(f)/f in R]   
    V=[(-f.subs(v=0))/(f.coefficient(v)) for f in L2] 
    L4=[[(p*s +diff(p,y) - diff(q,x))/q,s] for s in V]    
    return L4
def moses(R,p,q):
    r1 = (diff((diff(p,y)-diff(q,x))/q,y)==0)
    r2 = (diff((diff(q,x)-diff(p,y))/p,x)==0)
    return [r1,r2]
def lagutinski_uv_det(R,p,q,N):
    B = sorted (((1 + x + y )^(2*N) + v *(1 + x + y )^(2*N)). monomials(), reverse =0)
    r=v*(diff(p,y)*q-diff(q,y)*p) - (diff(p,y) - diff(q,x))*diff(q,y)+q*diff((diff(p,y) - diff(q,x)),y)
    D = lambda F: q^2*diff(F,x)-p*q*diff(F,y) + r*diff(F,v)
    L=lagutinski_det_random(R,D,B,N)
    return L
def lagutinski_uv_random(R,p,q,N):
    B = sorted (((1 + x + y )^(2*N) + v *(1 + x + y )^(2*N)). monomials(), reverse =0)
    r=v*(diff(p,y)*q-diff(q,y)*p) - (diff(p,y) - diff(q,x))*diff(q,y)+q*diff((diff(p,y) - diff(q,x)),y)
    D = lambda F: q^2*diff(F,x)-p*q*diff(F,y) + r*diff(F,v)
    x0=floor(random()*100)
    y0=floor(random()*100)
    v0=floor(random()*100)
    v1=floor(random()*100)
    L1=gcd([int(det(lagutinski_matrix(R,D,B,N).subs({x:x0,y:y0,v:v0}))), int(det(lagutinski_matrix(R,D,B,N+1).subs({x:x0,y:y0,v:v0}))), int(det(lagutinski_matrix(R,D,B,N+2).subs({x:x0,y:y0,v:v0})))])
    L2=gcd([int(det(lagutinski_matrix(R,D,B,N).subs({x:x0,y:y0,v:v1}))), int(det(lagutinski_matrix(R,D,B,N+1).subs({x:x0,y:y0,v:v1}))), int(det(lagutinski_matrix(R,D,B,N+2).subs({x:x0,y:y0,v:v1})))])
    return [L1.factor(),L2.factor()]
