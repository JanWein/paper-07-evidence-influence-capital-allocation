"""Exact rational polynomials and Sturm isolation. No optimizer or sampling grid."""
from fractions import Fraction as F
from math import gcd,lcm
from functools import reduce

def p(x):
    if isinstance(x,(int,float,F)):return [F(str(x))]
    return trim([F(str(z)) for z in x])
def trim(x):
    x=list(x)
    while len(x)>1 and not x[-1]:x.pop()
    return x or [F(0)]
def zero(x):return not any(x)
def add(a,b):
    a=p(a);b=p(b);return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(max(len(a),len(b)))])
def neg(a):return [-x for x in p(a)]
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    a=p(a);b=p(b);z=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]+=x*y
    return trim(z)
def scale(a,k):return mul(a,[k])
def powp(a,n):
    z=p(1)
    for _ in range(n):z=mul(z,a)
    return z
def der(a):return trim([i*a[i] for i in range(1,len(a))])
def val(a,x):
    out=0
    for y in reversed(a):out=out*x+y
    return out
def primitive(a):
    a=p(a)
    if zero(a):return a
    den=reduce(lcm,(x.denominator for x in a),1)
    nums=[int(x*den) for x in a];g=reduce(gcd,(abs(x) for x in nums),0)
    return [F(x//g) for x in nums] # positive scaling, preserving all signs

def divrem(a,b):
    a=p(a);b=p(b);assert not zero(b)
    q=[F(0)]*max(1,len(a)-len(b)+1)
    while not zero(a) and len(a)>=len(b):
        k=len(a)-len(b);c=a[-1]/b[-1];q[k]+=c
        a=sub(a,[F(0)]*k+scale(b,c))
    return trim(q),trim(a)
def roots01(a,bits=48):
    """All distinct real roots strictly between 0 and 1, exact Sturm counts."""
    a=primitive(a)
    if zero(a):return [] # identity handled as constant-branch case by caller
    while len(a)>1 and a[0]==0:a=a[1:]
    while len(a)>1 and val(a,F(1))==0:a,_=divrem(a,[-1,1])
    if len(a)<2:return []
    seq=[primitive(a),primitive(der(a))]
    while not zero(seq[-1]):
        _,rem=divrem(seq[-2],seq[-1])
        if zero(rem):break
        seq.append(primitive(neg(rem)))
    def variations(x):
        signs=[1 if v>0 else -1 for y in seq if (v:=val(y,x))!=0]
        return sum(a!=b for a,b in zip(signs,signs[1:]))
    initial=(F(0),F(1));stack=[(*initial,variations(initial[0])-variations(initial[1]))];out=[]
    while stack:
        lo,hi,count=stack.pop()
        if not count:continue
        if hi-lo<F(1,2**bits) and count==1:
            out.append(float((lo+hi)/2));continue
        mid=(lo+hi)/2;left=variations(lo)-variations(mid)
        stack.append((mid,hi,count-left));stack.append((lo,mid,left))
    return sorted(out)

# Polynomials in d with univariate-polynomial coefficients in tau.
def bp(x):return [p(x)]
def bt(a):
    a=list(a)
    while len(a)>1 and zero(a[-1]):a.pop()
    return a
def ba(a,b):
    return bt([add(a[i] if i<len(a) else p(0),b[i] if i<len(b) else p(0)) for i in range(max(len(a),len(b)))])
def bn(a):return [neg(x) for x in a]
def bs(a,b):return ba(a,bn(b))
def bm(a,b):
    z=[p(0) for _ in range(len(a)+len(b)-1)]
    for i,x in enumerate(a):
        for j,y in enumerate(b):z[i+j]=add(z[i+j],mul(x,y))
    return bt(z)
def bscale(a,k):return [mul(x,k) for x in a]
def quad_resultant(E,R):
    """Pseudo-reduce R modulo quadratic E, then exact quadratic-linear resultant.
    tau=0 can be an extraneous root; roots01 excludes it. Nonzero roots are filtered.
    """
    rem=bt(R)
    while len(rem)>2:
        k=len(rem)-1;lead=rem[-1]
        rem=bs(bscale(rem,E[2]),[p(0)]*(k-2)+bscale(E,lead))
    M=rem[0];L=rem[1] if len(rem)>1 else p(0)
    return primitive(add(sub(mul(E[2],mul(M,M)),mul(E[1],mul(L,M))),mul(E[0],mul(L,L))))
