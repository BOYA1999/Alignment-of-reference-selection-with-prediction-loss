from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.linalg import solve
from scipy.special import expit
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT/'artifacts/revision_002'
OLD = ROOT/'artifacts/confirmatory/confirm_001/full_recovery_001'
CACHE = ROOT/'artifacts/revision_001/forensic/local_feature_cache.npz'
MODELS = ['ecfp_rf', 'ecfp_desc_xgb']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    checks, solver_rows, summaries = [], [], []
    pred = pd.read_parquet(OLD/'predictions.parquet')
    replay = pd.read_parquet(OLD/'row_replay_seed41001.parquet')
    gen = pd.read_parquet(OLD/'generator_vectors.parquet')
    cache = np.load(CACHE)
    names = json.loads((OLD/'structure_audit.json').read_text())['descriptor_names']
    ipc = 2048+names.index('Ipc')
    for phase in ['original','bridge']:
        folder = BASE/phase
        estimates = pd.read_csv(folder/'estimates.csv',float_precision='round_trip')
        aggregate = pd.read_csv(folder/'summary.csv',float_precision='round_trip')
        manifest = json.loads((folder/'manifest.json').read_text())
        assert all(sha(OLD/n)==h for n,h in manifest['frozen_sha256'].items())
        assert all(sha(folder/n)==h for n,h in manifest['output_hashes'].items())
        assert len(estimates)==(30400 if phase=='original' else 14400)
        seen_states = set()
        max_error = max_expanded = 0.
        vector_hashes = {}
        for (ep,anchor,mechanism,seed), rows in estimates.groupby(['endpoint','anchor','mechanism','seed'],sort=True):
            vv = folder/'row_vectors'/ep/f'{anchor}_{mechanism}'/f'{seed}.npz'
            vector_hashes[vv.relative_to(folder).as_posix()] = sha(vv)
            v = np.load(vv)
            first = replay[replay.endpoint.eq(ep)].sort_values('source_row_index').drop_duplicates('source_row_index')
            ids, folds = first.source_row_index.to_numpy(), first.fold.to_numpy()
            y = np.column_stack([pred[pred.endpoint.eq(ep)&pred.model.eq(m)].sort_values('source_row_index').loss for m in MODELS])
            n = len(y)
            case = int(np.flatnonzero(ids==2945)[0])
            if phase=='original':
                p = gen[gen.endpoint.eq(ep)&gen.anchor.eq(anchor)&gen.mechanism.eq(mechanism)].sort_values('source_row_index').p_oracle.to_numpy()
            else:
                mw = first.audit_MolWt.to_numpy()
                score = (1 if mechanism.endswith('positive') else -1)*2*(mw-mw.mean())/mw.std()
                coverage = gen[gen.endpoint.eq(ep)&gen.anchor.eq(anchor)&gen.mechanism.eq('mar')].p_oracle.mean()
                left,right = -50.,50.
                for _ in range(200):
                    mid = (left+right)/2
                    if np.clip(expit(score+mid),.1,.9).mean()<coverage:
                        left = mid
                    else:
                        right = mid
                p = np.clip(expit(score+(left+right)/2),.1,.9)
            assert np.allclose(p,v['p_oracle'],atol=2e-15,rtol=0)
            r = (np.random.default_rng(seed).random(n)<p).astype(int)
            assert np.array_equal(r,v['R'])
            for row in rows.itertuples():
                k = MODELS.index(row.model)
                loss = y[:,k]
                if row.method=='complete_case':
                    value = np.sum(r*loss)/r.sum()
                else:
                    w = r/v['p_'+row.propensity]
                    if row.method=='ht':
                        value = np.dot(w,loss)/n
                    elif row.method=='hajek':
                        value = np.dot(w,loss)/w.sum()
                    else:
                        m = v['m_'+row.loss_variant][:,k]
                        expanded = (1-w)*m+w*loss
                        value = np.sum(expanded)/n
                        other = m.mean()+np.dot(w,loss-m)/n
                        max_expanded = max(max_expanded,abs(value-other)/max(1.,abs(value)))
                        assert np.isclose(row.case_pseudo_contribution,expanded[case]/n,rtol=1e-10,atol=1e-10)
                        assert np.isclose(row.case_abs_share,abs(expanded[case])/np.abs(expanded).sum(),rtol=1e-10,atol=1e-12)
                        assert np.isclose(row.case_m_coefficient,1-w[case])
                assert np.isclose(row.truth,loss.mean(),atol=1e-14,rtol=0)
                delta = abs(value-row.estimate)/max(1.,abs(row.estimate))
                max_error = max(max_error,delta)
                assert delta<1e-10
            statekey = (ep,anchor,mechanism,int(r[case]))
            if phase=='original' and statekey not in seen_states:
                seen_states.add(statekey)
                full = np.column_stack([cache['features'][np.searchsorted(cache['source_rows'],ids)],first.audit_max_train_similarity])
                train, valid = folds!=folds[case],folds==folds[case]
                for variant in ['raw','drop']:
                    x = full if variant=='raw' else np.delete(full,ipc,axis=1)
                    mu,sd = x[train].mean(0),x[train].std(0)
                    sd[sd==0] = 1
                    xx = ((x-mu)/sd)[train&r.astype(bool)]
                    yy = y[train&r.astype(bool)]
                    xm,ym = xx.mean(0),yy.mean(0)
                    xc,yc = xx-xm,yy-ym
                    dual = solve(xc@xc.T+10*np.eye(len(xc)),yc,assume_a='pos')
                    coef = xc.T@dual
                    dense = (((x[case]-mu)/sd)-xm)@coef+ym
                    normal = xc.T@(xc@coef-yc)+10*coef
                    residual = np.linalg.norm(normal)/max(1.,np.linalg.norm(xc.T@yc))
                    for k,model in enumerate(MODELS):
                        lsqr = v['m_'+variant][case,k]
                        solver_rows.append(dict(endpoint=ep,anchor=anchor,mechanism=mechanism,seed=seed,case_R=int(r[case]),
                            model=model,variant=variant,dense_case_prediction=dense[k],lsqr_case_prediction=lsqr,
                            relative_difference=abs(dense[k]-lsqr)/max(1.,abs(dense[k])),normal_equation_relative_residual=residual))
                    assert residual<1e-10
        (folder/'row_vector_hashes.json').write_text(json.dumps(vector_hashes,indent=2)+'\n')
        for row in aggregate.itertuples():
            mask = np.ones(len(estimates),dtype=bool)
            for key in ['endpoint','anchor','mechanism','model','method','propensity','loss_variant']:
                mask &= estimates[key].eq(getattr(row,key))
            part = estimates[mask]
            error = part.estimate.to_numpy()-part.truth.to_numpy()
            expected = dict(signed_bias=error.mean(),mc_sd=np.std(error,ddof=1),rmse=np.sqrt(np.mean(error**2)),
                            median_abs_error=np.median(abs(error)),mean_abs_error=np.mean(abs(error)))
            for key,value in expected.items():
                assert np.isclose(getattr(row,key),value,rtol=1e-10,atol=1e-12)
            summaries.append({'phase':phase,'rows':len(part),**expected})
        checks.append(dict(phase=phase,estimate_rows=len(estimates),vector_files=len(vector_hashes),summary_rows=len(aggregate),
                           max_relative_arithmetic_difference=max_error,max_relative_expanded_identity_difference=max_expanded))
    aa = BASE/'alignment'
    rep = pd.read_parquet(aa/'replicates.parquet')
    ag = pd.read_csv(aa/'summary.csv',float_precision='round_trip')
    des = pd.read_csv(aa/'design.csv')
    assert des.groupby(['endpoint','anchor','model']).sorted_probability_sha256.nunique().eq(1).all()
    for row in ag.itertuples():
        values = rep[(rep.endpoint==row.endpoint)&(rep.anchor==row.anchor)&(rep.model==row.model)&(rep.pairing==row.pairing)][row.method].to_numpy()
        e = values-row.truth
        assert len(values)==10000
        assert np.isclose(np.median(abs(e)),row.median_abs_error,rtol=1e-11,atol=1e-12)
        assert np.isclose(np.sqrt(np.mean(e*e)),row.rmse,rtol=1e-11,atol=1e-12)
    pd.DataFrame(solver_rows).to_csv(BASE/'independent_solver.csv',index=False,float_format='%.17g')
    report = dict(status='PASS_SPECIFIED_COMPUTATIONAL_CHECKS',not_external_human_validation=True,checks=checks,
        alignment_summary_rows=len(ag),solver_rows=len(solver_rows),solver_max_relative_residual=max(r['normal_equation_relative_residual'] for r in solver_rows),
        solver_max_relative_case_difference=max(r['relative_difference'] for r in solver_rows),
        frozen_files_unchanged=True,script_sha256=sha(Path(__file__)))
    (BASE/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    with threadpool_limits(limits=1):
        main()
