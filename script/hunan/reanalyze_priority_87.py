"""Reproducible 2015–2025 analysis of 87 observed sentinel site-years.

Stages: describe, models, finalize, all. No imputed outcome enters any likelihood.
Exact trajectory permutations preserve site observation-pattern groups.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import sysconfig
from fractions import Fraction

os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
if not (Path(sysconfig.get_paths()['include'])/'Python.h').exists():
    os.environ.setdefault('PYTENSOR_FLAGS','cxx=')
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'data'
OUT=ROOT/'results/priority_87_2026-10-22'
OUT.mkdir(parents=True,exist_ok=True)
STREAMS={'human':('positive','tested',-3.4),'culture':('leptospira_culture_positive','effective_cultures',-3.0)}
RULES={'Culture':('leptospira_culture_positive','effective_cultures'),'Capture':('rodents_caught','trap_setting_count'),'Compound':('leptospira_culture_positive','trap_setting_count')}
SEED=20261022

def dump(obj,path):
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')

def save(d,name):d.to_csv(OUT/name,index=False,encoding='utf-8-sig')

def load_panel():
    p=pd.read_csv(DATA/'analysis_panel.csv').sort_values(['year','site']).reset_index(drop=True)
    assert len(p)==87 and p.site.nunique()==8 and p.year.nunique()==11
    assert not ((p.site=='Yuanjiang')&(p.year==2018)).any()
    assert tuple(p[['positive','tested','leptospira_culture_positive','effective_cultures']].sum())==(291,8834,317,6824)
    assert not p.duplicated(['site','year']).any()
    return p

def aggregate(p):
    a=p.groupby('site').agg(human_pos=('positive','sum'),human_n=('tested','sum'),culture_pos=('leptospira_culture_positive','sum'),culture_n=('effective_cultures','sum'),caught=('rodents_caught','sum'),trap_effort=('trap_setting_count','sum'),observed_years=('year','nunique')).reset_index()
    for prefix,k,n in [('human','human_pos','human_n'),('culture','culture_pos','culture_n'),('capture','caught','trap_effort')]:
        a[prefix+'_fraction']=a[k]/a[n]
        z=1.959963984540054;q=a[prefix+'_fraction'];den=1+z*z/a[n]
        mid=(q+z*z/(2*a[n]))/den;half=z*np.sqrt(q*(1-q)/a[n]+z*z/(4*a[n]**2))/den
        a[prefix+'_low']=mid-half;a[prefix+'_high']=mid+half
    return a

def rho_p(x,y):
    x=rankdata(x);y=rankdata(y);xc=x-x.mean();yc=y-y.mean();scale=np.linalg.norm(xc)*np.linalg.norm(yc)
    vals=np.asarray(list(itertools.permutations(yc)))@xc/scale
    obs=float(yc@xc/scale)
    return {'rho':obs,'p_two_sided':float(np.mean(np.abs(vals)>=abs(obs)-1e-12)),'permutations':len(vals)}

def selections(a,num,den,k):
    scores={s:Fraction(int(a.loc[s,num]),int(a.loc[s,den])) for s in a.index}
    threshold=sorted(scores.values(),reverse=True)[k-1]
    required=[s for s in a.index if scores[s]>threshold];tied=[s for s in a.index if scores[s]==threshold]
    sets=[tuple(sorted(required+list(extra))) for extra in itertools.combinations(tied,k-len(required))]
    weights=pd.Series({s:sum(s in v for v in sets)/len(sets) for s in a.index})
    assert np.isclose(weights.sum(),k)
    return sets,weights

def descriptive(p,tag):
    sites=sorted(p.site.unique());years=sorted(p.year.unique());a=aggregate(p)
    save(a,f'{tag}_site_summary.csv')
    ann=p.groupby('year')[['positive','tested','leptospira_culture_positive','effective_cultures','rodents_caught','trap_setting_count']].sum().reset_index()
    ann['sites']=p.groupby('year').site.nunique().to_numpy();save(ann,f'{tag}_annual_counts.csv')
    corr={s:rho_p(a.human_fraction,a[s+'_fraction']) for s in ['culture','capture']}
    fold=[];weight_rows=[];allweights={rule:np.zeros((len(years),len(sites))) for rule in RULES}
    for yi,yr in enumerate(years):
        test=p[p.year==yr].set_index('site');train=p[p.year!=yr].groupby('site').sum(numeric_only=True).loc[sorted(test.index)]
        for rule,(num,den) in RULES.items():
            sets,w=selections(train,num,den,3);chosen=w[w>0].index
            ps=float((test.positive*w).sum());ns=float((test.tested*w).sum());pa=int(test.positive.sum());na=int(test.tested.sum())
            fold.append({'held_out_year':int(yr),'rule':rule,'eligible_sites':len(test),'selected_sites':';'.join(chosen),'admissible_sets':len(sets),'all_tied_sets':'|'.join(';'.join(v) for v in sets),'positive_selected':ps,'positive_all':pa,'tested_selected':ns,'tested_all':na,'coverage':ps/pa,'testing_share':ns/na,'enrichment':(ps/pa)/(ns/na)})
            for s in sites:
                v=float(w.get(s,0));allweights[rule][yi,sites.index(s)]=v
                weight_rows.append({'year':int(yr),'site':s,'rule':rule,'selection_probability':v,'available':s in test.index})
    fold=pd.DataFrame(fold)
    if tag=='primary87':
        # Selected sets, not the full culture ranking, are fold-invariant.
        # Capture's eight-site order is invariant; 2018 only substitutes the fourth site.
        assert set(fold.loc[fold.rule=='Culture','selected_sites'])=={'Lianyuan;Shuangfeng;Xiangtan'}
        assert set(fold.loc[fold.rule=='Compound','selected_sites'])=={'Lianyuan;Shuangfeng;Xiangtan'}
        cap=fold.loc[fold.rule=='Capture',['held_out_year','selected_sites']]
        assert set(cap.loc[cap.held_out_year!=2018,'selected_sites'])=={'Liuyang;Shuangfeng;Yuanjiang'}
        assert set(cap.loc[cap.held_out_year==2018,'selected_sites'])=={'Lianyuan;Liuyang;Shuangfeng'}
        orders={}
        for yr in years:
            test=p[p.year==yr].set_index('site')
            train=p[p.year!=yr].groupby('site').sum(numeric_only=True).loc[sorted(test.index)]
            for rule,(num,den) in RULES.items():
                scores={s:Fraction(int(train.loc[s,num]),int(train.loc[s,den])) for s in train.index}
                orders.setdefault(rule,[]).append((int(yr),tuple(sorted(scores,key=lambda s:(-scores[s],s)))))
        cap_full=[order for yr,order in orders['Capture'] if yr!=2018]
        assert len(set(cap_full))==1 and cap_full[0][3]=='Lianyuan'
        cult_2020=dict(orders['Culture'])[2020]
        assert cult_2020[:3]==('Xiangtan','Shuangfeng','Lianyuan')
    save(fold,f'{tag}_loyo_by_year.csv');save(pd.DataFrame(weight_rows),f'{tag}_loyo_selection.csv')
    pool=fold.groupby('rule')[['positive_selected','positive_all','tested_selected','tested_all']].sum().reset_index()
    pool['coverage']=pool.positive_selected/pool.positive_all;pool['testing_share']=pool.tested_selected/pool.tested_all;pool['enrichment']=pool.coverage/pool.testing_share
    save(pool,f'{tag}_loyo_pooled.csv')
    # A permutation may only move complete human trajectories between sites with
    # the same observed-year pattern. Yuanjiang stays fixed in the 87-row panel.
    y=p.pivot(index='year',columns='site',values='positive').reindex(index=years,columns=sites)
    n=p.pivot(index='year',columns='site',values='tested').reindex(index=years,columns=sites)
    patterns={}
    for j,s in enumerate(sites): patterns.setdefault(tuple(y[s].notna()),[]).append(j)
    groups=list(patterns.values());perms=[]
    for group_perms in itertools.product(*[list(itertools.permutations(g)) for g in groups]):
        mapping=np.arange(len(sites))
        for g,pg in zip(groups,group_perms):mapping[g]=pg
        perms.append(mapping)
    perms=np.asarray(perms);yv=y.fillna(0).to_numpy();nv=n.fillna(0).to_numpy()
    pc={};pn={}
    for rule,w in allweights.items():
        pc[rule]=np.einsum('yps,ys->p',yv[:,perms],w)
        pn[rule]=np.einsum('yps,ys->p',nv[:,perms],w)
    delta_e=(pc['Culture']/pn['Culture']-pc['Capture']/pn['Capture'])*nv.sum()/yv.sum()
    delta_c=(pc['Culture']-pc['Capture'])/yv.sum()
    pooled=pool.set_index('rule');obs_e=float(pooled.loc['Culture','enrichment']-pooled.loc['Capture','enrichment']);obs_c=float(pooled.loc['Culture','coverage']-pooled.loc['Capture','coverage'])
    perm={'n_permutations':len(perms),'exchangeable_groups':[[sites[i] for i in g] for g in groups]}
    for name,values,obs in [('enrichment',delta_e,obs_e),('coverage',delta_c,obs_c)]:
        center=float(np.median(values));perm[name]={'observed_difference':obs,'null_median':center,'null_interval':np.quantile(values,[.025,.975]).tolist(),'p_two_sided':float(np.mean(abs(values-center)>=abs(obs-center)-1e-12))}
    save(pd.DataFrame({'delta_enrichment':delta_e,'delta_coverage':delta_c}),f'{tag}_permutations.csv')
    contribution=[]
    for yi,yr in enumerate(years):
        for j,s in enumerate(sites):
            if pd.notna(y.loc[yr,s]):contribution.append({'year':int(yr),'site':s,'human_positive':int(y.loc[yr,s]),'contribution':float((allweights['Culture'][yi,j]-allweights['Capture'][yi,j])*y.loc[yr,s])})
    save(pd.DataFrame(contribution),f'{tag}_contributions.csv')
    result={'n':len(p),'years':[int(v) for v in years],'totals':{c:int(p[c].sum()) for c in ['positive','tested','leptospira_culture_positive','effective_cultures','rodents_caught','trap_setting_count']},'correlations':corr,'permutation':perm,'loyo':pool.to_dict('records')}
    dump(result,OUT/f'{tag}_summary.json')
    return result

def describe():
    p=load_panel();save(p,'analysis_panel.csv')
    sets={'primary87':p,'balanced80':p[p.year!=2018],'early63':p[p.year<=2022],'recent24':p[p.year>=2023],'exclude2015':p[p.year!=2015]}
    alt=p.copy();alt.effective_cultures=alt.cultures_total-alt.contaminated
    sets['arithmetic_effective']=alt
    summaries={tag:descriptive(frame,tag) for tag,frame in sets.items()}
    dump(summaries,OUT/'descriptive_sensitivity_summary.json')
    a=aggregate(p).set_index('site');sites=sorted(a.index);influence=[]
    for excluded in ['None']+sites:
        d=a if excluded=='None' else a.drop(excluded)
        influence.append({'excluded_site':excluded,**{f'{stream}_{key}':v for stream in ['culture','capture'] for key,v in rho_p(d.human_fraction,d[stream+'_fraction']).items()}})
    save(pd.DataFrame(influence),'leave_one_site_out.csv')
    ports=[];rules=[]
    for k in range(1,8):
        for combo in itertools.combinations(sites,k):
            d=a.loc[list(combo)];c=d.human_pos.sum()/a.human_pos.sum();t=d.human_n.sum()/a.human_n.sum()
            ports.append({'k':k,'sites':';'.join(combo),'coverage':c,'testing_share':t,'enrichment':c/t})
        for rule,num,den in [('Culture','culture_pos','culture_n'),('Capture','caught','trap_effort')]:
            combos,w=selections(a,num,den,k);vals=[a.loc[list(c),'human_pos'].sum()/a.human_pos.sum() for c in combos]
            rules.append({'k':k,'rule':rule,'n_tied_sets':len(combos),'coverage':np.mean(vals),'minimum':min(vals),'maximum':max(vals),'selected_sets':'|'.join(';'.join(v) for v in combos)})
    save(pd.DataFrame(ports),'all_portfolios.csv');save(pd.DataFrame(rules),'pooled_rule_frontier.csv')
    print(json.dumps(summaries['primary87'],indent=2),flush=True)

def basis(n):return np.linalg.qr(np.vstack([np.eye(n-1),-np.ones((1,n-1))]))[0][:,:n-1]

def fit_model(p,stream,family,label,site_tau=.7,year_tau=.5,phi_rate=.1,draws=1500,tune=1500):
    import arviz as az
    import pymc as pm
    import pytensor.tensor as pt
    sites=sorted(p.site.unique());years=sorted(p.year.unique());p=p.sort_values(['year','site']).reset_index(drop=True)
    si=pd.Categorical(p.site,categories=sites).codes;yi=pd.Categorical(p.year,categories=years).codes
    sb,yb=basis(len(sites)),basis(len(years));kcol,ncol,prior=STREAMS[stream]
    k=p[kcol].to_numpy(int);n=p[ncol].to_numpy(int)
    signature={'counts_sha256':hashlib.sha256(p[['site','year',kcol,ncol]].to_csv(index=False).encode()).hexdigest(),'family':family,'site_tau':site_tau,'year_tau':year_tau,'phi_rate':phi_rate,'draws':draws,'tune':tune,'chains':4,'sampler':'nutpie','sites':sites,'years':[int(y) for y in years]}
    nc=OUT/f'{label}.nc';meta=OUT/f'{label}_fit.json'
    if nc.exists() and meta.exists() and json.loads(meta.read_text())==signature:
        print('USE completed posterior',label,flush=True);return
    seed=SEED+int(hashlib.sha256(label.encode()).hexdigest()[:6],16)%100000
    signature['seed']=seed
    # Metadata seed is deterministic; use it when resuming.
    if nc.exists() and meta.exists() and json.loads(meta.read_text())==signature:
        print('USE completed posterior',label,flush=True);return
    print('FIT',label,'observations',len(p),'seed',seed,flush=True)
    with pm.Model(coords={'obs':np.arange(len(p)),'site_contrast':np.arange(len(sites)-1),'year_contrast':np.arange(len(years)-1)}) as model:
        intercept=pm.Normal('intercept',prior,2)
        st=pm.Normal('site_theta',0,site_tau,dims='site_contrast');yt=pm.Normal('year_theta',0,year_tau,dims='year_contrast')
        mu=pm.Deterministic('mu',pm.math.sigmoid(intercept+pt.dot(sb[si],st)+pt.dot(yb[yi],yt)),dims='obs')
        if family=='betabinomial':
            phi=pm.Gamma('concentration',alpha=2,beta=phi_rate)
            pm.BetaBinomial('y',n=n,alpha=mu*phi,beta=(1-mu)*phi,observed=k,dims='obs')
        else:pm.Binomial('y',n=n,p=mu,observed=k,dims='obs')
        idata=pm.sample(draws=draws,tune=tune,chains=4,cores=4,target_accept=.95,nuts_sampler="nutpie",random_seed=seed,progressbar=False,return_inferencedata=True,compute_convergence_checks=True)
    idata.attrs.update({'label':label,'site_order':json.dumps(sites),'year_order':json.dumps([int(y) for y in years]),'counts_sha256':signature['counts_sha256']})
    # nutpie records nested inference metadata; NetCDF attributes must be scalar/array types.
    for attrs in [idata.attrs]+[getattr(idata,g).attrs for g in idata.groups()]:
        for key,value in list(attrs.items()):
            if isinstance(value,dict): attrs[key]=json.dumps(value,default=str)
    idata.to_netcdf(nc);dump(signature,meta)
    diagnostics=az.summary(idata,var_names=['intercept','site_theta','year_theta']+(['concentration'] if family=='betabinomial' else []),kind='diagnostics',round_to='none')
    save(diagnostics.reset_index(names='parameter'),f'{label}_diagnostics.csv')
    print('DONE',label,'max Rhat',diagnostics.r_hat.max(),'min bulk ESS',diagnostics.ess_bulk.min(),'divergences',int(idata.sample_stats.diverging.sum()),flush=True)

def models():
    p=load_panel()
    for family in ['betabinomial','binomial']:
        for stream in STREAMS:fit_model(p,stream,family,f'{stream}_{family}')
    for label,kw in [('site035',{'site_tau':.35}),('site140',{'site_tau':1.4}),('year100',{'year_tau':1.0}),('phi_mean10',{'phi_rate':.2}),('phi_mean40',{'phi_rate':.05})]:
        fit_model(p,'human','betabinomial',f'human_sensitivity_{label}',draws=1000,tune=1000,**kw)
    arithmetic=p.copy();arithmetic.effective_cultures=arithmetic.cultures_total-arithmetic.contaminated
    fit_model(arithmetic,'culture','betabinomial','culture_arithmetic_effective',draws=1000,tune=1000)
    for stream in STREAMS:
        fit_model(p[p.year!=2018],stream,'betabinomial',f'{stream}_balanced80',draws=1000,tune=1000)
        fit_model(p[p.year<=2022],stream,'betabinomial',f'{stream}_early63',draws=1000,tune=1000)

def posterior(label):
    import arviz as az
    z=az.from_netcdf(OUT/f'{label}.nc');meta=json.loads((OUT/f'{label}_fit.json').read_text());post=z.posterior
    sb=basis(len(meta['sites']));yb=basis(len(meta['years']))
    st=post.site_theta.values.reshape(-1,sb.shape[1])@sb.T;yt=post.year_theta.values.reshape(-1,yb.shape[1])@yb.T
    rates=expit(post.intercept.values.reshape(-1,1,1)+st[:,:,None]+yt[:,None,:]).mean(axis=2)
    return z,meta,rates

def distribution(x):return {'median':float(np.median(x)),'mean':float(np.mean(x)),'low':float(np.quantile(x,.025)),'high':float(np.quantile(x,.975))}

def finalize():
    import arviz as az
    p=load_panel();sites=sorted(p.site.unique());agg=aggregate(p).set_index('site').loc[sites];site_records=[];rank_records=[];model_records=[];joint_records=[];ppc_records=[];cell_records=[];rho_draws=[];diagnostic_records=[]
    cache={};rng=np.random.default_rng(SEED+500)
    for label in [q.stem.replace('_fit','') for q in sorted(OUT.glob('*_fit.json'))]:
        z,meta,rates=posterior(label);cache[label]=rates
        diag=pd.read_csv(OUT/f'{label}_diagnostics.csv');diagnostic_records.append({'fit':label,'max_rhat':float(diag.r_hat.max()),'min_bulk_ess':float(diag.ess_bulk.min()),'min_tail_ess':float(diag.ess_tail.min()),'divergences':int(z.sample_stats.diverging.sum()),'draws':rates.shape[0],'max_tree_depth':int(z.sample_stats['tree_depth' if 'tree_depth' in z.sample_stats else 'depth'].max())})
        stream='human' if label.startswith('human') else 'culture';family=meta['family'];kcol,ncol,_=STREAMS[stream]
        frame=p[p.year.isin(meta['years'])].copy();frame=frame.sort_values(['year','site'])
        if label=='culture_arithmetic_effective': frame.effective_cultures=frame.cultures_total-frame.contaminated
        observed=frame[kcol].to_numpy(int);n=frame[ncol].to_numpy(int)
        assert np.array_equal(z.observed_data.y.values,observed)
        mu=z.posterior.mu.values.reshape(-1,len(frame))
        phi=z.posterior.concentration.values.reshape(-1,1) if family=='betabinomial' else None
        prob=rng.beta(mu*phi,(1-mu)*phi) if phi is not None else mu
        rep=rng.binomial(n,prob);low,high=np.quantile(rep,[.025,.975],axis=0);pbar=mu.mean(axis=0)
        po=np.sum((observed-n*pbar)**2/(n*pbar*(1-pbar)+1e-12));pr=np.sum((rep-n*pbar)**2/(n*pbar*(1-pbar)+1e-12),axis=1)
        rec={'fit':label,'stream':stream,'family':family,'cells':len(frame),'observed_total':int(observed.sum()),'observed_zeros':int((observed==0).sum()),'total_low':float(np.quantile(rep.sum(axis=1),.025)),'total_high':float(np.quantile(rep.sum(axis=1),.975)),'zeros_low':float(np.quantile((rep==0).sum(axis=1),.025)),'zeros_high':float(np.quantile((rep==0).sum(axis=1),.975)),'ppp_zeros':float(np.mean((rep==0).sum(axis=1)>=(observed==0).sum())),'ppp_pearson':float(np.mean(pr>=po)),'cells_covered':int(((observed>=low)&(observed<=high)).sum())}
        if phi is not None:rec.update({f'concentration_{k}':v for k,v in distribution(phi.ravel()).items()})
        ppc_records.append(rec)
        if label not in [f'{s}_{f}' for s in STREAMS for f in ['binomial','betabinomial']]:continue
        save(pd.DataFrame({'total':rep.sum(axis=1),'zeros':(rep==0).sum(axis=1)}),f'{label}_predictive_totals.csv')
        order=np.argsort(-rates,axis=1);ranks=np.argsort(order,axis=1)+1
        for j,s in enumerate(sites):
            site_records.append({'model':family,'stream':stream,'site':s,**distribution(rates[:,j]),'top3':float(np.mean(ranks[:,j]<=3))})
            for k in range(1,9):rank_records.append({'model':family,'stream':stream,'site':s,'rank':k,'probability':float(np.mean(ranks[:,j]==k))})
        for i,row in enumerate(frame.itertuples()):cell_records.append({'model':family,'stream':stream,'site':row.site,'year':row.year,'index':i+1,'observed':int(observed[i]),'low':float(low[i]),'high':float(high[i])})
    for family in ['binomial','betabinomial']:
        h,c=cache[f'human_{family}'],cache[f'culture_{family}'];rh=rankdata(h,axis=1);rc=rankdata(c,axis=1);corr=((rh-4.5)*(rc-4.5)).sum(axis=1)/42
        oh=np.argsort(-h,axis=1);oc=np.argsort(-c,axis=1);same=np.all(np.sort(oh[:,:3],axis=1)==np.sort(oc[:,:3],axis=1),axis=1)
        model_records.append({'model':family,**{f'rho_{k}':v for k,v in distribution(corr).items()},'probability_positive':float(np.mean(corr>0)),'same_top3':float(same.mean())})
        for j,s in enumerate(sites):joint_records.append({'model':family,'site':s,'joint_top3':float(np.mean(np.any(oh[:,:3]==j,axis=1)&np.any(oc[:,:3]==j,axis=1)))})
        rho_draws.extend({'model':family,'rho':float(v)} for v in corr)
    window_records=[]
    for scenario,hn,cn in [('balanced80','human_balanced80','culture_balanced80'),('early63','human_early63','culture_early63'),('arithmetic_effective','human_betabinomial','culture_arithmetic_effective')]:
        h,c=cache[hn],cache[cn];nd=min(len(h),len(c));h,c=h[:nd],c[:nd]
        rh,rc=rankdata(h,axis=1),rankdata(c,axis=1);rv=((rh-4.5)*(rc-4.5)).sum(axis=1)/42
        window_records.append({'scenario':scenario,**distribution(rv),'probability_positive':float(np.mean(rv>0))})
    save(pd.DataFrame(window_records),'posterior_correlation_sensitivity.csv')
    for name,records in [('posterior_site_rates.csv',site_records),('posterior_ranks.csv',rank_records),('posterior_joint.csv',joint_records),('posterior_correlations.csv',model_records),('posterior_rho_draws.csv',rho_draws),('posterior_predictive_checks.csv',ppc_records),('posterior_predictive_cells.csv',cell_records),('sampling_diagnostics.csv',diagnostic_records)]:save(pd.DataFrame(records),name)
    # Decision utilities use fixed observed site testing volumes and all 56 triples.
    combos=list(itertools.combinations(range(8),3));weights=agg.human_n.to_numpy();candidate={}
    for name,num,den in [('Culture','culture_pos','culture_n'),('Capture','caught','trap_effort')]:
        chosen,_=selections(agg,num,den,3);assert len(chosen)==1
        candidate[name]=combos.index(tuple(sorted(sites.index(s) for s in chosen[0])))
    decisions=[];contrasts=[];portfolio_rows=[]
    for label,rates in cache.items():
        if not label.startswith('human'):continue
        w=weights.copy()
        if label.endswith('balanced80'): w=aggregate(p[p.year!=2018]).set_index('site').loc[sites,'human_n'].to_numpy()
        if label.endswith('early63'): w=aggregate(p[p.year<=2022]).set_index('site').loc[sites,'human_n'].to_numpy()
        weighted=rates*w;u=np.stack([weighted[:,list(c)].sum(axis=1)/weighted.sum(axis=1) for c in combos],axis=1);best=u.max(axis=1)
        selected={**candidate,'Posterior optimum':int(u.mean(axis=0).argmax())}
        for name,j in selected.items():
            vals=u[:,j];regret=best-vals
            decisions.append({'fit':label,'strategy':name,'sites':';'.join(sites[i] for i in combos[j]),**distribution(vals),'probability_best':float(np.mean(vals==best)),'expected_regret':float(regret.mean()),'regret_95':float(np.quantile(regret,.95))})
        vals=u.mean(axis=1);decisions.append({'fit':label,'strategy':'Uniform random','sites':'all 56 equally weighted',**distribution(vals),'probability_best':1/56,'expected_regret':float(np.mean(best-vals)),'regret_95':float(np.quantile(best-vals,.95))})
        delta=u[:,candidate['Culture']]-u[:,candidate['Capture']]
        contrasts.append({'fit':label,**distribution(delta),'probability_culture_better':float(np.mean(delta>0))})
        if label=='human_betabinomial':
            save(pd.DataFrame({'delta_utility':delta}),'paired_utility_draws.csv')
            for j,c in enumerate(combos):portfolio_rows.append({'sites':';'.join(sites[i] for i in c),'expected_utility':float(u[:,j].mean()),'expected_regret':float(np.mean(best-u[:,j])),'probability_optimal':float(np.mean(u[:,j]==best))})
    save(pd.DataFrame(decisions),'posterior_decisions.csv');save(pd.DataFrame(contrasts),'posterior_decision_sensitivity.csv');save(pd.DataFrame(portfolio_rows),'posterior_all_triples.csv')
    dump({'primary_panel':87,'years':11,'sampling_diagnostics':diagnostic_records,'models':model_records,'decision_contrasts':contrasts,'primary_count_hash':hashlib.sha256((OUT/'analysis_panel.csv').read_bytes()).hexdigest()},OUT/'model_summary.json')
    print(pd.DataFrame(model_records).to_string(index=False),flush=True)
    print(pd.DataFrame(contrasts).to_string(index=False),flush=True)
    bad=[r for r in diagnostic_records if r['max_rhat']>1.01 or r['min_bulk_ess']<400 or r['divergences']>0]
    if bad:raise RuntimeError('Sampling diagnostics require attention: '+json.dumps(bad))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['describe','models','finalize','all'],default='all');arg=ap.parse_args()
    if arg.stage in ['describe','all']:describe()
    if arg.stage in ['models','all']:models()
    if arg.stage in ['finalize','all']:finalize()
