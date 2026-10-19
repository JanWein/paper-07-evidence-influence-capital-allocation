#!/usr/bin/env python3
"""Paper 07 R2. Synthetic equilibria only. Python standard library, no dependencies."""
from dataclasses import dataclass, replace
from math import sqrt, isclose
import json, csv, random
from pathlib import Path
@dataclass(frozen=True)
class Params:
    p:float=.4; v:float=2.; ell:float=2.; r:float=1.; f:float=.2; C:float=1.2
    k:float=.02; e:float=.03; a:float=.04; K:float=.06
    def validate(self):
        assert 0<self.p<1 and min(self.v,self.ell,self.r,self.f,self.C)>0
        assert self.f/self.C < self.p*self.v/((1-self.p)*self.ell) < self.r/self.C < 1
        assert min(self.k,self.e,self.a,self.K)>=0

def h(x,d):
    assert d>0
    return max(0.,min(1.,(x+d)/(2*d)))
def B(x,d):
    if d==0:return max(0.,x)
    if x<=-d:return 0.
    if x>=d:return x
    return (x+d)**2/(4*d)
def nu(s,z):return (z.p*z.v-(1-z.p)*z.ell*s)/(z.p+(1-z.p)*s)
def bisect(fun,lo,hi):
    # Independent monotone root check, not used for the closed-form solution.
    for _ in range(100):
        mid=(lo+hi)/2
        if fun(mid)>0:hi=mid
        else:lo=mid
    return (lo+hi)/2

def audit_s(d,z):
    if d==0:return z.p*z.v/((1-z.p)*z.ell)
    rho=z.r/z.C
    aa=2*d*(1-z.p)*rho
    bb=2*d*z.p+(max(z.ell-d,0))*(1-z.p)*rho
    cc=z.p*(z.v+min(d,z.ell))
    return rho*(2*cc/(bb+sqrt(bb*bb+4*aa*cc)))

def mixing_interval(d,z,lower):
    n=h(-z.ell,d); b=z.f/z.r
    if abs(n+b-1)<2e-12:
        upper=min(1.,z.p*(z.v-d)/((1-z.p)*(z.ell+d)))
        return (lower,max(lower,upper))
    target=d*(2*(n+b)-1)
    s=z.p*(z.v-target)/((1-z.p)*(z.ell+target))
    return (max(lower,min(1.,s)),)*2

def solve(d,A,L,z=Params(),endpoint=0):
    z.validate(); assert d>=0 and A in (0,1) and L in (0,1)
    rho=z.r/z.C; alpha=z.f/z.C; b=z.f/z.r
    if d==0:
        T=z.p*z.v/((1-z.p)*z.ell)
        s=0 if A and L else T
        M=0 if L else (T if A else alpha)
        I=0 if A else s-M
        n=0.; u=1. if A and L else (T/rho if A else b)
        interval=(s,s); regime='baseline'
    else:
        n=h(-z.ell,d)
        D=lambda s:h(nu(s,z),d)-n
        if A:
            s=0. if L else audit_s(d,z)
            M=s; I=0.; interval=(s,s);regime='revealing' if L else 'metric-only'
        else:
            lower=0. if L else alpha
            if D(1)>b+1e-12:
                s=1.; M=lower; I=s-M; interval=(s,s);regime='pooling'
            elif D(lower)<b-1e-12:
                s=0. if L else audit_s(d,z)
                M=0. if L else s;I=0.;interval=(s,s);regime='revealing' if L else 'metric-only'
            else:
                interval=mixing_interval(d,z,lower)
                s=interval[endpoint];M=lower;I=max(0.,s-M);regime='mixed' if interval[0]==interval[1] else 'plateau'
        u=h(nu(s,z),d)
    P=z.p+(1-z.p)*s
    Q=z.p+(1-z.p)*(M+I)
    prep=(1-z.p)*(z.C*M*M/2+z.f*I)
    gross=P*B(nu(s,z),d)+(1-z.p)*(1-s)*B(-z.ell,d)
    W=gross-z.k-(z.e+A*z.a)*Q-prep-L*z.K
    return dict(A=A,L=L,delta=d,s=s,M=M,I=I,N=1-M-I,u=u,n=n,
                posterior=z.p/P,packet=Q,gross=gross,prep=prep,W=W,
                no_evidence=B(nu(1,z),d),regime=regime,s_interval=interval)

def counterfactual_accounting(d,A,L,M,I,z):
    """All information sets and FAILED packets, also outside equilibrium."""
    s=M+(1-A)*I; passed=z.p+(1-z.p)*s
    failed=(1-z.p)*A*I; negative=(1-z.p)*(1-M-I)
    Q=passed+failed
    gross=passed*B(nu(s,z),d)+(failed+negative)*B(-z.ell,d)
    return gross-z.k-(z.e+A*z.a)*Q-(1-z.p)*(z.C*M*M/2+z.f*I)-L*z.K

def check(d,A,L,z):
    for endpoint in (0,1):
        q=solve(d,A,L,z,endpoint);s=q['s'];n=q['n'];u=q['u'];M=q['M'];I=q['I']
        assert min(M,I,q['N'])>-1e-9 and abs(M+I+q['N']-1)<1e-9
        if d:
            assert abs(u-h(nu(s,z),d))<1e-10
            if A and not L:
                alt=bisect(lambda x:x-z.r/z.C*(h(nu(x,z),d)-n),0,z.r/z.C)
                assert abs(s-alt)<1e-9
        # Check individual best responses over private types, with point cutoffs excluded.
        for j in range(301):
            c=z.C*(j+.371)/301
            UN=z.r*n; UI=z.r*(n if A else u)-z.f
            UM=z.r*u-c if not L else -1e100
            best=max(UN,UI,UM)
            if c<z.C*M-1e-9: chosen=UM
            elif I>1e-9 and (q['N']<1e-9 or UI>UN+1e-8):chosen=UI
            else: chosen=UN if q['N']>1e-9 else UI
            assert best-chosen<2e-8,(d,A,L,q,c,best,chosen)
            if I>1e-8:assert UI>=max(UN,UM if c>=z.f else -1e100)-1e-8 or c<z.f
        assert abs(q['W']-counterfactual_accounting(d,A,L,M,I,z))<1e-9
        if d:
            # Independent midpoint integration of receiver payoff on each information set.
            # Integrate conditional expected value, not realized quality-based decisions.
            def quadrature(x):
                count=4000
                return sum(max(x-(-d+2*d*(j+.5)/count),0) for j in range(count))/count
            P=z.p+(1-z.p)*s
            numeric=P*quadrature(nu(s,z))+(1-z.p)*(1-s)*quadrature(-z.ell)
            assert abs(q['gross']-numeric)<max(1e-7,d/2000000)
    return q

def main():
    root=Path(__file__).resolve().parents[1];out=root/'results';out.mkdir(exist_ok=True)
    z=Params();rows=[]
    grid=[0,.0001,.01,.1,.3,.5,2/3,1,2,3,4,5,6,7,8,10,20,100]
    for d in grid:
        for A,L in [(0,0),(1,0),(0,1),(1,1)]:rows.append(check(d,A,L,z))
    # Genuine aggregate multiplicity at upper clipping with locked metrics.
    plateau=Params(p=.14,v=10.,ell=2.,r=1.,f=.75,C=1.1)
    q=check(4.,0,1,plateau)
    assert q['regime']=='plateau'
    assert abs(q['s_interval'][0])<1e-10
    assert abs(q['s_interval'][1]-.14*6/(.86*6))<1e-10
    (out/'clipping_multiplicity.json').write_text(json.dumps({'parameters':plateau.__dict__,'delta':4.,'endpoints':[solve(4.,0,1,plateau,j) for j in (0,1)]},indent=2)+'\n')
    rng=random.Random(7017);random_tests=0;valuable=0
    for _ in range(150):
        p=rng.uniform(.1,.65);ell=rng.uniform(.5,4);r=rng.uniform(.1,3);rho=rng.uniform(.2,.95)
        T=rng.uniform(.04,rho*.98);alpha=rng.uniform(.001,T*.98)
        zz=Params(p=p,v=T*(1-p)*ell/p,ell=ell,r=r,f=alpha*r/rho,C=r/rho,k=0,e=0,a=0,K=0)
        d=10**rng.uniform(-3,2)*ell
        old=check(d,0,0,zz);new=check(d,1,0,zz)
        check(d,0,1,zz);check(d,1,1,zz);random_tests+=4
        if old['W']>=old['no_evidence']-1e-10:
            valuable+=1;assert new['W']>=old['W']-1e-9
    # Exact and adjacent endpoints of all regimes in the displayed reference family.
    for d in [2/3,4.,8.,10.]:
        for dd in [d*(1-1e-7),d,d*(1+1e-7)]:
            for A,L in [(0,0),(1,0),(0,1),(1,1)]:check(dd,A,L,z)
    # Positive fees and positive resource costs in both policy examples.
    cases={'audited_imperfect_optimum':(1.,replace(z,K=1.)),
           'unaudited_imperfect_optimum_with_fee_loss':(.1,replace(z,r=.02,f=.004,C=.024,k=.0001,e=.0001,a=.1,K=1.))}
    examples={}
    for name,(d,zz) in cases.items():
        vals=[solve(d,A,L,zz) for A,L in [(0,0),(1,0),(0,1),(1,1)]]
        q0,qA=vals[:2];P=qA['packet']
        threshold=(qA['gross']-q0['gross']+q0['prep']-qA['prep']+zz.e*(q0['packet']-P))/P
        examples[name]={'parameters':zz.__dict__,'delta':d,'no_evidence':q0['no_evidence'],
          'policies':vals,'audit_fee_threshold':threshold,'audit_costless_change':threshold*P}
        target=1 if name.startswith('audited') else 0
        assert vals[target]['W']>max([q0['no_evidence']]+[x['W'] for i,x in enumerate(vals) if i!=target])
    # Failed packets still paid for and may be FUNDED at bad posterior when d>ell.
    zz=replace(z,k=0,e=.03,a=.04,K=0);d=5
    assert h(-zz.ell,d)>0
    value=counterfactual_accounting(d,1,0,.2,.3,zz)
    assert abs(value-(counterfactual_accounting(d,1,0,.2,0,zz)-(1-zz.p)*.3*(zz.f+zz.e+zz.a)))<1e-10
    (out/'examples.json').write_text(json.dumps(examples,indent=2)+'\n')
    with (out/'equilibria.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=rows[0]);writer.writeheader();writer.writerows(rows)
    report={'status':'PASS','deterministic_policy_rows':len(rows),'random_policy_checks':random_tests,
            'valuable_random_baselines':valuable,'checks':'closed root versus independent bisection; individual deviations; all information sets; quadrature; regime endpoints; costly failed packets; strict policy examples; boundary implication'}
    (out/'checks.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    print(json.dumps(examples,indent=2))
if __name__=='__main__':main()
