var("x,y,dx")
f=x^3-y^3+2*x*y+x-2*y+1
x1=0
x2=1
y1=QQbar[y](f.subs(x=x1)).roots(multiplicities=False)[0]
y2=QQbar[y](f.subs(x=x2)).roots(multiplicities=False)[0]

def iii_eqs(f,p1,p2):
    [x1,y1]=p1
    [x2,y2]=p2
    r=QQ[x,y](f).degree()
    m=[x^i*y^j for i in range(r) for j in range(r) if i+j<r]
    c=var(['c'+str(i) for i in range(len(m))])
    E=sum([c[i]*m[i] for i in range(len(m))])
    R1=[z for z in QQbar[y](f.subs(x=x1)).roots(multiplicities=False)]
    R2=[z for z in QQbar[y](f.subs(x=x2)).roots(multiplicities=False)]
    eqs1=[]
    for z in R1:
        if z==y1:
            eqs1.append((E-(x2-x1)*diff(f,y)).subs([x==x1,y==y1]))
        else:
            eqs1.append(E.subs([x==x1,y==z]))
    eqs2=[]
    for z in R2:
        if z==y2:
            eqs2.append((E-(x2-x1)*diff(f,y)).subs([x==x2,y==y2]))
        else:
            eqs2.append(E.subs([x==x2,y==z]))    
    return eqs1+eqs2

def iii_eqs_sym(f,p1,p2):
    [x1,y1]=p1
    [x2,y2]=p2
    r=QQ[x,y](f).degree()
    m=[x^i*y^j for i in range(r) for j in range(r) if i+j<r]
    c=var(['c'+str(i) for i in range(len(m))])
    E=sum([c[i]*m[i] for i in range(len(m))])
    R1=[z for z in QQbar[y](f.subs(x=x1)).roots(multiplicities=False)]
    R2=[z for z in QQbar[y](f.subs(x=x2)).roots(multiplicities=False)]
    eqs1=[]
    for z in R1:
        if z==y1:
            eqs1.append((E-(x2-x1)*diff(f,y)).subs([x==x1,y==y1]))
        else:
            eqs1.append(E.subs([x==x1,y==z]))
    teqs1=[sum([R1[s]^m*eqs1[s] for s in range(len(R1))]) for m in range(r)]
    eqs2=[]
    for z in R2:
        if z==y2:
            eqs2.append((E-(x2-x1)*diff(f,y)).subs([x==x2,y==y2]))
        else:
            eqs2.append(E.subs([x==x2,y==z]))    
    teqs2=[sum([R2[s]^m*eqs2[s] for s in range(len(R2))]) for m in range(r)]
    return teqs1+teqs2
    
def iii(f,p1,p2):
    [x1,y1]=p1
    [x2,y2]=p2
    r=QQ[x,y](f).degree()
    m=[x^i*y^j for i in range(r) for j in range(r) if i+j<r]
    c=var(['c'+str(i) for i in range(len(m))])
    E=sum([c[i]*m[i] for i in range(len(m))])
    R1=[z for z in QQbar[y](f.subs(x=x1)).roots(multiplicities=False)]
    R2=[z for z in QQbar[y](f.subs(x=x2)).roots(multiplicities=False)]
    eqs1=[]
    for z in R1:
        if z==y1:
            eqs1.append((E-(x2-x1)*diff(f,y)).subs([x==x1,y==y1]))
        else:
            eqs1.append(E.subs([x==x1,y==z]))
    teqs1=[sum([R1[s]^m*eqs1[s] for s in range(len(R1))]) for m in range(r)]
    eqs2=[]
    for z in R2:
        if z==y2:
            eqs2.append((E-(x2-x1)*diff(f,y)).subs([x==x2,y==y2]))
        else:
            eqs2.append(E.subs([x==x2,y==z]))    
    teqs2=[sum([R2[s]^m*eqs2[s] for s in range(len(R2))]) for m in range(r)]
    eqs=teqs1+teqs2
    S=lsolve(eqs,c)
    ans=E.subs(S)/(x-x1)/(x2-x)/diff(f,y)*dx
    return ans
       
def haupt_fuction_eval(f,p1,p2,pp,A):
    u=iii(f,p1,p2)
    r=QQ[x,y](f).degree()
    m=[x^i*y^j for i in range(r) for j in range(r) if i+j<r-2]
    if len(m)==1:
        c=[var("c0")]
    else:
        c=var(['c'+str(i) for i in range(len(m))])
    E=sum([c[i]*m[i] for i in range(len(m))])
    eqs=[u.subs([x==a,dx==1]).subs(y=b) + (E*dx/diff(f,y)).subs([x==a,dx==1]).subs(y=b) for [a,b] in A]
    S=lsolve(eqs,c)
    E = QQbar(sum([(c[i]).subs(S)*(m[i]).subs([x==pp[0], y==pp[1]]) for i in range(len(m))]))
    d = QQbar((1/diff(f,y)).subs(x=pp[0]).subs(y=pp[1]))
    u = QQbar(u.subs(dx=1).subs(x=pp[0]).subs(y=pp[1]))
    return u+E*d

def lsolve(eqs,x,K=QQbar,T=False):
    n=len(eqs)
    m=len(x)
    A=matrix(K, n, m, lambda i, j: eqs[i].coefficient(x[j]))
    b=matrix(K, n, 1, lambda i, j: -eqs[i].subs([xx==0 for xx in x]))
    x_eval=A \ b
    if T==True:
        ans=x_eval
    else:
        ans=[x[i]==K(x_eval[i][0]) for i in range(m)]
    return ans
