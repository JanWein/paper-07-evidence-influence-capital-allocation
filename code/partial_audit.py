#!/usr/bin/env python3
"""R3 open-metric random audits. All results are synthetic.
Global selection uses exact rational polynomial candidates/Sturm isolation, not a grid.
"""
from dataclasses import replace
from fractions import Fraction as F
from pathlib import Path
import json,csv,random,sys
from polynomials import *
from reproduce_R2 import Params,h,B,nu,audit_s,solve as legacy_solve


def solve(tau,delta,z):
    z.validate();assert delta>=0 and 0<=tau<=1
    if delta==0:return solve_zero(tau,z)
    rho=z.r/z.C;alpha=z.f/z.C;b=z.f/z.r;n=h(-z.ell,delta)
    D=lambda s:h(nu(s,z),delta)-n
    y=audit_s(delta,z);dy=y/rho;star=max(0.,1-b/dy)
    if tau>=star-1e-13:
        s=y;M=y;I=0.;reg='metric-only'
    else:
        d0=b/(1-tau);upper=1-tau+tau*rho*d0
        if D(upper)<=d0+1e-13:
            x=delta*(2*(n+d0)-1)
            s=z.p*(z.v-x)/((1-z.p)*(z.ell+x));M=rho*d0;I=max(0.,(s-M)/(1-tau));reg='mix'
        else:
            lo=0.;hi=1.
            for _ in range(90):
                mid=(lo+hi)/2
                if mid-(1-tau+tau*alpha+tau*tau*rho*D(mid))>0:hi=mid
                else:lo=mid
            s=(lo+hi)/2;di=D(s);M=alpha+tau*rho*di;I=1-M
            reg='strict-upper' if nu(s,z)>=delta-1e-12 else 'strict-interior'
    di=D(s);u=n+di;Q=z.p+(1-z.p)*(M+I);P=z.p+(1-z.p)*s
    G=P*B(nu(s,z),delta)+(1-z.p)*(1-s)*B(-z.ell,delta)
    prep=(1-z.p)*(z.C*M*M/2+z.f*I)
    W=G-z.k-z.e*Q-z.a*tau*Q-prep
    return dict(tau=tau,delta=delta,regime=reg,s=s,M=M,I=I,N=1-M-I,d=di,u=u,n=n,
        Q=Q,passed=P,failed=(1-z.p)*tau*I,posterior=z.p/P,gross=G,prep=prep,
        evaluation=z.e*Q,audit=z.a*tau*Q,production=z.k,W=W,y=y,tau_star=star,
        no_evidence=B(nu(1,z),delta))

def solve_zero(t,z):
    rho=z.r/z.C;al=z.f/z.C;b=z.f/z.r;T=z.p*z.v/((1-z.p)*z.ell)
    star=1-al/T
    if t>=star-1e-13:s=T;u=T/rho;M=T;I=0.;reg='metric-only'
    else:
        d0=b/(1-t);upper=1-t+t*rho*d0
        if T<=upper+1e-13:s=T;u=d0;M=rho*d0;I=(s-M)/(1-t);reg='mix'
        else:
            ut=(T-1+t*(1-al))/(rho*t*t)
            if ut<=1:s=T;u=ut;reg='strict-indifferent'
            else:u=1.;s=1-t+t*al+rho*t*t;reg='strict-upper'
            M=al+t*rho*u;I=1-M
    Q=z.p+(1-z.p)*(M+I);P=z.p+(1-z.p)*s
    G=u*(z.p*z.v-(1-z.p)*z.ell*s);prep=(1-z.p)*(z.C*M*M/2+z.f*I)
    W=G-z.k-z.e*Q-z.a*t*Q-prep
    return dict(tau=t,delta=0.,regime=reg,s=s,M=M,I=I,N=1-M-I,d=u,u=u,n=0.,Q=Q,
       passed=P,failed=(1-z.p)*t*I,posterior=z.p/P,gross=G,prep=prep,evaluation=z.e*Q,
       audit=z.a*t*Q,production=z.k,W=W,y=T,tau_star=star,no_evidence=0.)

def global_zero(z):
    zz=rational(z);pp=zz['p'];r=zz['r'];C=zz['C'];rho=r/C;al=zz['f']/C;T=pp*zz['v']/((1-pp)*zz['ell'])
    star=float(1-al/T);cand=[(0.,'zero'),(star,'deterrence'),(1.,'full')]
    for t in roots01([T-1,2-T-al,-1]):
        if t<=star+1e-9 and abs(solve_zero(t,z)['N'])<1e-7:cand.append((t,'mix/strict'))
    for t in roots01([1-T,-(1-al),rho]):
        if t<=star+1e-9 and solve_zero(t,z)['regime'].startswith('strict'):cand.append((t,'clipping'))
    peak=((1-pp)*zz['ell']*(1-al)-zz['a'])/((1-pp)*(2*zz['ell']*rho+C*rho*rho))
    if 0<float(peak)<star and solve_zero(float(peak),z)['regime']=='strict-upper':cand.append((float(peak),'upper-stationary'))
    rows=[]
    for t,kind in sorted(cand):
        if any(abs(t-x['tau'])<1e-10 for x in rows):continue
        row=solve_zero(t,z);row['candidate']=kind;rows.append(row)
    best=max(x['W'] for x in rows)
    return dict(candidates=rows,best=max(rows,key=lambda x:x['W']),all_optimal_candidates=[x['tau'] for x in rows if abs(x['W']-best)<1e-10],stationary_resultant_degree=0,stationary_resultant_coefficients=[],constant_interior_branch=False,optimal_intervals=([[star,1.]] if z.a==0 and abs(solve_zero(star,z)['W']-best)<1e-10 else []))

def derivative(t,delta,z):
    rr=solve(t,delta,z);d=rr['d'];rho=z.r/z.C;al=z.f/z.C;n=rr['n'];q=z.p/(1-z.p)
    lam=z.p*(z.v+z.ell)/(1-z.p);L=max(z.ell-delta,0);H=L+2*delta*d
    if delta==0:
        S=z.p*z.v/((1-z.p)*z.ell);Sp=0.
    else:S=lam/H-q;Sp=-2*delta*lam/H**2
    Gd=2*(1-z.p)*lam*delta*(L*(d+n)+delta*d*d)/H**2
    if rr['regime']=='metric-only':return -z.a*rr['Q']
    if rr['regime']=='mix':
        b=z.f/z.r;Qp=(1-z.p)*(rho+(S+d*Sp-2*rho*d)/b)
        wd=Gd-(z.e+z.a*t)*Qp-z.a*b*rr['Q']/d**2-(1-z.p)*z.r*(S+d*Sp-rho*d)
        return wd*d*d/b
    if rr['regime']=='strict-indifferent':
        return -z.a-(1-z.p)*z.C*((1-al)-(1-S)/t)*(1-S)/t**2
    if rr['regime']=='strict-upper':
        sp=-(1-al)+2*rho*t*d
        return -(1-z.p)*(z.ell+B(-z.ell,delta))*sp-z.a-(1-z.p)*z.C*rho*rho*t*d*d
    dp=(1-al-2*rho*t*d)/(rho*t*t-Sp)
    J=(1-z.p)*z.C*rho*rho
    return (Gd-J*t*t*d)*dp-z.a-J*t*d*d

def rational(z):return {k:F(str(v)) for k,v in z.__dict__.items()}
def global_candidates(delta,z):
    """Construct ALL candidate families in the theorem, exact rational coefficients.
    Every algebraic root is isolated by Sturm; spurious algebraic branches are filtered
    against actual equilibrium and stationary/regime conditions.
    """
    if delta==0:return global_zero(z)
    zz=rational(z);pp=zz['p'];r=zz['r'];C=zz['C'];f=zz['f'];ell=zz['ell'];a=zz['a'];de=F(str(delta))
    rho=r/C;al=f/C;b=f/r;q=pp/(1-pp);lam=pp*(zz['v']+ell)/(1-pp)
    n=max(F(0),(de-ell)/(2*de));L=max(ell-de,F(0));dm=1-n
    H=p([L,2*de]);dc=p([0,1]);tc=p([0,1])
    y=audit_s(delta,z);star=max(0.,1-float(b)/(y/float(rho)))
    cand=[(0.,'zero'),(star,'deterrence'),(1.,'full')]
    # Mixing/strict transition cubic T(d), transformed d=b/(1-tau).
    Td=sub(scale(dc,lam),mul(add(add(scale(dc,q-b*rho),scale(powp(dc,2),rho)),b),H))
    # Exact substitution clears (1-tau)^degree.
    def substit_d(poly):
        deg=len(poly)-1;res=p(0)
        for j,c in enumerate(poly):res=add(res,scale(powp([1,-1],deg-j),c*b**j))
        return primitive(res)
    for t in roots01(substit_d(Td)):
        if t<=star+1e-9:
            rr=solve(t,delta,z)
            if abs(rr['N'])<1e-7 and abs((1-t)*rr['d']-float(b))<1e-7:cand.append((t,'mix/strict'))
    # Strict/interior equation E(tau,d)=(q+1-(1-alpha)t+rho*t^2*d)H-lambda.
    V=p([q+1,-(1-al)]);tt=powp(tc,2)
    E=[sub(scale(V,L),lam),add(scale(V,2*de),scale(tt,L*rho)),scale(tt,2*de*rho)]
    # Upper-clipping transitions: E(tau,dmax)=0.
    clip=p(0)
    for j,co in enumerate(E):clip=add(clip,scale(co,dm**j))
    for t in roots01(clip):
        rr=solve(t,delta,z)
        if t<star+1e-9 and rr['regime'].startswith('strict') and abs(nu(rr['s'],z)-delta)<1e-7:cand.append((t,'clipping'))
    # Strict upper-clipped stationary point, an explicitly concave quadratic.
    Bbad=de*n*n
    peak=((1-pp)*(ell+Bbad)*(1-al)-a)/((1-pp)*(2*(ell+Bbad)*rho*dm+C*rho*rho*dm*dm))
    if 0<float(peak)<star:
        rr=solve(float(peak),delta,z)
        if rr['regime']=='strict-upper':cand.append((float(peak),'upper-stationary'))
    # Stationary polynomial R for the strict interior branch, in d with tau coeffs.
    HH=bm([p(L),p(2*de)],[p(L),p(2*de)])
    J=(1-pp)*C*rho*rho
    GN=[p(2*(1-pp)*lam*de*L*n),p(2*(1-pp)*lam*de*L),p(2*(1-pp)*lam*de*de)]
    part=bs(GN,[p(0)]+bscale(HH,scale(tt,J)))
    NN=[p(1-al),scale(tc,-2*rho)]
    denom=ba(bscale(HH,scale(tt,rho)),bp(2*de*lam))
    RR=bs(bm(part,NN),bm([p(a),p(0),scale(tc,J)],denom))
    resultant=quad_resultant(E,RR)
    # If resultant is identically zero, irreducibility of E makes W constant on
    # its interior branch. Endpoints already cover that case (proof in manuscript).
    root_audit=[]
    for t in roots01(resultant):
        checkrow=solve(t,delta,z)
        root_audit.append({'tau':t,'actual_regime':checkrow['regime'],'actual_derivative':derivative(t,delta,z),'feasible_interior_stationary':bool(t<star-1e-10 and checkrow['regime']=='strict-interior' and abs(derivative(t,delta,z))<2e-6)})
        if t<star-1e-10:
            rr=solve(t,delta,z)
            if rr['regime']=='strict-interior' and abs(derivative(t,delta,z))<2e-6:cand.append((t,'interior-stationary'))
    rows=[]
    for t,kind in sorted(cand):
        if any(abs(t-x['tau'])<1e-10 for x in rows):continue
        row=solve(t,delta,z);row['candidate']=kind;rows.append(row)
    best=max(x['W'] for x in rows)
    intervals=[]
    if z.a==0 and abs(solve(star,delta,z)['W']-best)<1e-10:intervals.append([star,1.])
    if zero(resultant):
        edges=sorted(set([0.,star]+[t for t,kind in cand if kind in ['mix/strict','clipping']]))
        for lo,hi in zip(edges,edges[1:]):
            mid=solve((lo+hi)/2,delta,z)
            if mid['regime']=='strict-interior' and abs(mid['W']-best)<1e-10:intervals.append([lo,hi])
    return dict(candidates=rows,best=max(rows,key=lambda x:x['W']),optimal_intervals=intervals,stationary_root_audit=root_audit,
        all_optimal_candidates=[x['tau'] for x in rows if abs(x['W']-best)<1e-10],
        stationary_resultant_degree=len(resultant)-1,
        stationary_resultant_coefficients=[str(x) for x in resultant],
        constant_interior_branch=zero(resultant))


def raw_check(t,delta,z):
    rr=solve(t,delta,z);M=rr['M'];I=rr['I'];N=rr['N']
    assert min(M,I,N)>-1e-8 and abs(M+I+N-1)<1e-8
    # Enumerate quality, action, and Bernoulli audit outcomes BEFORE observing messages.
    # Receiver sees only +, failed, negative. Checked pass and unchecked pass are pooled.
    events=[]
    for quality,action,mass in [('G','G',z.p),('B','M',(1-z.p)*M),('B','I',(1-z.p)*I),('B','N',(1-z.p)*N)]:
        if action=='N':events.append((quality,action,False,'negative',mass));continue
        for checked,prob in [(False,1-t),(True,t)]:
            hist='failed' if checked and action=='I' else '+'
            events.append((quality,action,checked,hist,mass*prob))
    assert abs(sum(x[4] for x in events)-1)<1e-10
    gross=0.;fee=0.;obs={}
    for hist in ['+','failed','negative']:
        subevents=[x for x in events if x[3]==hist];mass=sum(x[4] for x in subevents)
        pg=sum(x[4] for x in subevents if x[0]=='G')/mass if mass>1e-14 else (0 if hist!='+' else 1)
        mean=pg*z.v-(1-pg)*z.ell
        if mass>1e-14:
            assert abs(pg-(rr['posterior'] if hist=='+' else 0))<1e-9
        # Raw quality payoffs under the SAME history-based decision, not omniscience.
        nt=600
        value=0.
        for j in range(nt):
            shock=-delta+2*delta*(j+.5)/nt
            funded=(rr['u'] if hist=='+' and abs(mean)<1e-10 else float(mean>0)) if delta==0 else float(mean>shock)
            value+=sum(prob*((z.v if quality=='G' else -z.ell)-shock)*funded for quality,action,checked,hist0,prob in subevents)/nt
        gross+=value;obs[hist]={'mass':mass,'good_posterior':pg,'fund_probability':((rr['u'] if hist=='+' else 0.) if delta==0 else h(mean,delta))}
    # Exact cost enumeration over physical packets and completed checks.
    fee=sum(prob*(z.e*(action!='N')+z.a*checked) for quality,action,checked,hist,prob in events)
    prep=(1-z.p)*(z.C*M*M/2+z.f*I)
    assert abs(fee-rr['evaluation']-rr['audit'])<1e-10
    assert abs(gross-rr['gross'])<max(1e-6,delta/100000)
    assert abs(gross-fee-prep-z.k-rr['W'])<max(1e-6,delta/100000)
    # Each type's action is checked against both deviations from physical audit outcomes.
    u=rr['u'];n=rr['n'];cut=z.C*M
    for j in range(151):
        c=z.C*(j+.271)/151
        vals={'N':z.r*n,'I':z.r*((1-t)*u+t*n)-z.f,'M':z.r*u-c}
        best=max(vals.values())
        if c<cut-1e-8:chosen='M'
        elif rr['regime'].startswith('strict'):chosen='I'
        else:chosen='N' if N>1e-8 else 'I'
        assert best-vals[chosen]<2e-8,(rr,c,vals,chosen)
        if I>1e-8 and c>cut+1e-8:assert best-vals['I']<2e-8
    if rr['regime'].startswith('strict'):
        assert abs(cut-(z.f+t*z.r*rr['d']))<1e-8
    if rr['regime']=='mix':assert abs((1-t)*z.r*rr['d']-z.f)<1e-8
    return obs


def main():
    out=Path(__file__).resolve().parents[1]/'results';out.mkdir(exist_ok=True)
    z=Params()
    known=mul(mul([-1,3],[-2,3]),mul([-1,3],[-1,1]))
    isolated=roots01(known)
    assert len(isolated)==2 and abs(isolated[0]-1/3)<1e-12 and abs(isolated[1]-2/3)<1e-12
    cases={
      'R2_A':(1.,replace(z,K=1.)),
      'R2_B':(.1,replace(z,r=.02,f=.004,C=.024,k=.0001,e=.0001,a=.1,K=1.)),
      'interior_partial_inflation':(.02,Params(p=.5,v=1.6,ell=2.,r=.1,f=.00525,C=.105,k=.001,e=.001,a=.005,K=1.)),
      'zero_shock_partial_inflation':(0.,Params(p=.5,v=1.6,ell=2.,r=.1,f=.00525,C=.105,k=.001,e=.001,a=.005,K=1.)),
      'zero_audit_fee_partial_inflation':(.02,Params(p=.5,v=1.6,ell=2.,r=.1,f=.00525,C=.105,k=.001,e=.001,a=0.,K=1.)),
      'zero_audit_fee_tie':(1.,replace(z,K=1.,a=0.)),
      'wide_shock_negative_funding':(5.,replace(z,K=1.)),
      'inflation_absent_at_zero':(12.,z)}
    output={};checked=0;derivative_checks=0
    for name,(delta,zz) in cases.items():
        glob=global_candidates(delta,zz);star=solve(0,delta,zz)['tau_star']
        # Endpoint consistency with the preserved R2 solution.
        for t,A in [(0.,0),(1.,1)]:
            old=legacy_solve(delta,A,0,zz);new=solve(t,delta,zz)
            for key in ['s','M','I','W','gross','prep']:assert abs(old[key]-new[key])<1e-9
        rows=glob['candidates']
        boundaries=sorted(set([0.,1.,star]+[x['tau'] for x in rows if x['candidate'] in ['mix/strict','clipping']]))
        probe=sorted(set([x['tau'] for x in rows]+[(a+b)/2 for a,b in zip(boundaries,boundaries[1:])]))
        physical=[]
        for t in probe:
            obs=raw_check(t,delta,zz);checked+=1;physical.append({'tau':t,'histories':obs})
            if t>1e-5 and t<1-1e-5:
                eps=1e-6
                if solve(t-eps,delta,zz)['regime']==solve(t+eps,delta,zz)['regime']:
                    numeric=(solve(t+eps,delta,zz)['W']-solve(t-eps,delta,zz)['W'])/(2*eps)
                    assert abs(numeric-derivative(t,delta,zz))<2e-5;derivative_checks+=1
        # Main menu expands open audits only; locked choices remain the R2 0/1 choices.
        menu={'no_evidence':solve(0,delta,zz)['no_evidence'],
              'optimized_open':glob['best']['W'],
              'locked_no_audit':legacy_solve(delta,0,1,zz)['W'],
              'locked_full_audit':legacy_solve(delta,1,1,zz)['W']}
        glob['parameters']=zz.__dict__;glob['delta']=delta;glob['deterrence']=solve(star,delta,zz)
        glob['physical_checks']=physical;glob['menu']=menu;glob['optimal_menu']=max(menu,key=menu.get)
        output[name]=glob
    assert output['zero_audit_fee_tie']['optimal_intervals']==[[output['zero_audit_fee_tie']['deterrence']['tau_star'],1.]]
    assert abs(output['R2_A']['deterrence']['W']-.11901537693455363)<1e-10
    assert abs(output['R2_B']['deterrence']['W']+.015522826262085593)<1e-10
    assert abs(output['zero_audit_fee_partial_inflation']['best']['tau']-399/820)<1e-12
    assert output['zero_audit_fee_partial_inflation']['best']['W']>0
    best=output['interior_partial_inflation']['best']
    assert best['candidate']=='upper-stationary' and 0<best['tau']<best['tau_star'] and best['I']>0
    # Every candidate and all actual branch midpoints of the finite partition checked.
    (out/'partial_audit.json').write_text(json.dumps(output,indent=2)+'\n')
    with (out/'partial_candidates.csv').open('w',newline='') as f:
        rows=[dict(case=name,**x) for name,g in output.items() for x in g['candidates']]
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    report={'status':'PASS','exact_global_cases':len(cases),'raw_state_checks':checked,'derivative_checks':derivative_checks,
      'global_method':'exact rational polynomial candidate construction, Sturm root isolation, equilibrium filtering; no optimizer or grid',
      'audit_status_information':'receiver observes pooled positive, failed, or unchanged negative, not checked-pass flag',
      'attribution':'R2 A/B deterrence checks supplied by central review, reproduced here; not author-origin discoveries'}
    (out/'checks_R3.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
    for name,g in output.items():print(name,'tau',g['best']['tau'],'W',g['best']['W'],'candidate',g['best']['candidate'],'menu',g['optimal_menu'],'degree',g['stationary_resultant_degree'])
if __name__=='__main__':main()
