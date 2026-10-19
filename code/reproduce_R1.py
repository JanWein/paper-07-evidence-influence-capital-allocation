"""R1 Bayesian equilibrium checks. Run: python code/reproduce.py. No real data."""
from pathlib import Path
import csv,json,platform
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
import scipy
O=Path(__file__).resolve().parents[1]/'results'/'r1';O.mkdir(parents=True,exist_ok=True)
p=.4; v=2.; ell=2.; r=1.; f=.2; C=1.2; kappa=.02; e=.03; a=.04; K=.06
class Dist:
 def __init__(self,name):self.name=name
 def F(self,x):
  t=np.clip(x/C,0,1)
  return {'uniform':t,'rising_triangle':t*t,'falling_triangle':2*t-t*t}[self.name]
 def pdf(self,x):
  if x<0 or x>C:return 0.
  return {'uniform':1/C,'rising_triangle':2*x/C**2,'falling_triangle':2*(1-x/C)/C}[self.name]
 def quantile(self,t):return brentq(lambda x:self.F(x)-t,0,C,xtol=1e-14)
 def I(self,x):return quad(lambda c:c*self.pdf(c),0,x,epsabs=1e-12)[0]
def solve(dist,p=p):
 T=p*v/((1-p)*ell);alpha=dist.F(f);ct=dist.quantile(T)
 assert 0<T<1 and 0<f<r and alpha<T<dist.F(r)
 P=p+(1-p)*T;If=dist.I(f);It=dist.I(ct);J=If+f*(T-alpha)
 rows=[]
 for audit,locked in [(0,0),(1,0),(0,1),(1,1)]:
  if audit and locked:u=1.;m=0.;b=0.;n=1.;cost=0.;s=0.;post=1.;gross=p*v;pack=p
  elif locked:u=f/r;m=0.;b=T;n=1-T;cost=f*T;s=T;post=p/(p+(1-p)*T);gross=0.;pack=P
  elif audit:u=ct/r;m=T;b=0.;n=1-T;cost=It;s=T;post=p/(p+(1-p)*T);gross=0.;pack=P
  else:u=f/r;m=alpha;b=T-alpha;n=1-T;cost=J;s=T;post=p/(p+(1-p)*T);gross=0.;pack=P
  W=gross-kappa-(e+audit*a)*pack-(1-p)*cost-locked*K
  rows.append(dict(distribution=dist.name,audit=audit,locked=locked,u=u,metric_bad=m,
   inflation_bad=b,negative_bad=n,posterior_positive=post,positive_probability=pack,
   good_funded=p*u,bad_funded=(1-p)*s*u,gross_value=gross,
   preparation_cost=(1-p)*cost,production_cost=kappa,evaluation_cost=e*pack,
   audit_cost=audit*a*pack,governance_cost=locked*K,net_value=W,
   informed_allocation_regret=p*v-gross))
  for c in np.linspace(0,C,2401):
   payN=0.;payF=(1-audit)*r*u-f;payM=r*u-c if not locked else -np.inf
   best=max(payN,payF,payM)
   if audit and locked:chosen=payN
   elif locked:chosen=payF
   elif audit:chosen=payM if c<ct else payN
   else:chosen=payM if c<f else payF
   assert best-chosen<1e-10,(dist.name,audit,locked,c,best,chosen)
  value=post*v-(1-post)*ell
  assert abs(value)<1e-10 or (u==1 and value>=0)
  assert abs(p*u*v-(1-p)*s*u*ell-gross)<1e-10
 delta=f*alpha-If
 pred_audit=-(1-p)*quad(lambda c:(c-f)*dist.pdf(c),f,ct)[0]-a*P
 assert abs((rows[1]['net_value']-rows[0]['net_value'])-pred_audit)<1e-12
 interaction=rows[3]['net_value']-rows[2]['net_value']-rows[1]['net_value']+rows[0]['net_value']
 pred_inter=p*v+(e+a)*(1-p)*T+(1-p)*(It+delta)
 assert abs(interaction-pred_inter)<1e-12
 assert np.argmax([0]+[x['net_value'] for x in rows]) in (0,4)
 return rows,dict(distribution=dist.name,T=T,alpha=alpha,ct=ct,J=J,It=It,
  audit_effect_open=pred_audit,metric_only_effect=-(1-p)*delta-K,
  interaction=interaction,K_active_threshold=p*v+e*(1-p)*T-a*p+(1-p)*J,
  K_optout_threshold=p*(v-e-a)-kappa)
allrows=[];summ=[]
for name in ['uniform','rising_triangle','falling_triangle']:
 rows,s=solve(Dist(name));allrows+=rows;summ.append(s)
for filename,rows in [('policies.csv',allrows),('identities.csv',summ)]:
 with (O/filename).open('w') as h:
  wr=csv.DictWriter(h,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
# Genuine sponsor signal with P(S=+)=P(S=-)=.5 and conditional priors .45 and .35.
signal=[]
for ps in (.35,.45):signal+=solve(Dist('uniform'),p=ps)[0]
blind=solve(Dist('uniform'))[0];comparison=[]
for j in range(4):
 vis=.5*signal[j]['net_value']+.5*signal[j+4]['net_value']
 comparison.append(dict(audit=blind[j]['audit'],locked=blind[j]['locked'],
  blind_value=blind[j]['net_value'],visible_value=vis,visible_minus_blind=vis-blind[j]['net_value']))
with (O/'sponsor_signal.csv').open('w') as h:
 wr=csv.DictWriter(h,fieldnames=comparison[0]);wr.writeheader();wr.writerows(comparison)
assert comparison[1]['visible_minus_blind']<0
max_root_error=0.
for name in ['uniform','rising_triangle','falling_triangle']:
 d=Dist(name);T=p*v/((1-p)*ell)
 root=brentq(lambda u:p*v-(1-p)*ell*d.F(r*u),1e-10,1.)
 max_root_error=max(max_root_error,abs(root-d.quantile(T)/r))
assert max_root_error<1e-10
# Retained synthetic off-policy fixed-scale bounds example.
y=np.array([3.,0.,0.,-1.]);obs=np.array([1,0,0,1],bool);dd=np.array([-1.,1.,1.,-1.])
lo=float(np.sum(dd[obs]*y[obs])+np.minimum(dd[~obs]*-2,dd[~obs]*4).sum())
hi=float(np.sum(dd[obs]*y[obs])+np.maximum(dd[~obs]*-2,dd[~obs]*4).sum())
assert (lo,hi)==(-6.,6.)
meta=dict(status='synthetic Bayesian equilibrium calculations; no firm data',
 python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
 tested_cost_points_per_policy=2401,policies_per_distribution=4,distributions=3,
 max_equilibrium_root_error=max_root_error,off_policy_bounds=[lo,hi],
 parameters=dict(p=p,v=v,ell=ell,r=r,f=f,cost_support=C,kappa=kappa,e=e,a=a,K=K),
 policy_optimum='full audit plus locked metric for supplied synthetic parameters')
(O/'validation.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(meta,indent=2));print('Uniform values:',[x['net_value'] for x in allrows[:4]])
print('Identities:',summ);print('Sponsor signal:',comparison)
# Strict evaluator decisions: private continuous outside-option value t~U[-d,d].
# This changes equilibria; no old actions are imposed on the full-audit equilibrium.
regular=[]
for name in ['uniform','rising_triangle','falling_triangle']:
 D=Dist(name);T=p*v/((1-p)*ell);al=D.F(f);If=D.I(f);u0=f/r
 for d in (.001,.01,.05):
  def capacity(u):
   t=d*(2*u-1)
   return p*(v-t)/((1-p)*(ell+t))
  s0=capacity(u0)
  assert al<s0<1
  ua=brentq(lambda u:D.F(r*u)-capacity(u),u0,1)
  sa=D.F(r*ua)
  Ws=[]
  for audit,locked in [(0,0),(1,0),(0,1),(1,1)]:
   if audit and locked:u=1.;s=0.;prep=0.;P=p;gross=p*v
   else:
    u=ua if audit else u0;s=sa if audit else s0;P=p+(1-p)*s
    prep=(1-p)*(D.I(r*u) if audit else (f*s if locked else If+f*(s-al)))
    gross=P*d*u*u
    mu=(p*v-(1-p)*ell*s)/P
    assert abs(mu-d*(2*u-1))<1e-9
    # Integrate strict ex post receiver payoff, splitting at optimal cutoff.
    integrated=quad(lambda t:max(mu-t,0)/(2*d),-d,d,points=[mu])[0]*P
    assert abs(integrated-gross)<1e-10
   W=gross-kappa-(e+audit*a)*P-prep-locked*K;Ws.append(W)
   regular.append(dict(distribution=name,delta=d,audit=audit,locked=locked,u=u,
     bad_positive=s,net_value=W,gross_value=gross,preparation_cost=prep))
  assert Ws[1]<Ws[0] and Ws[2]<Ws[0]
  assert Ws[3]-Ws[2]-Ws[1]+Ws[0]>0
  assert max(0,*Ws)==max(0,Ws[3])
with (O/'strict_receiver.csv').open('w') as h:
 wr=csv.DictWriter(h,fieldnames=regular[0]);wr.writeheader();wr.writerows(regular)
meta['strict_receiver_equilibria']=len(regular)
(O/'validation.json').write_text(json.dumps(meta,indent=2)+'\n')
print('Strict evaluator extension:',len(regular),'policy equilibria checked')
