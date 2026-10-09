from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from revision_risk_diagnostics import metrics, probability_loss

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'artifacts/confirmatory/confirm_001/full_recovery_001'
OUT = ROOT/'artifacts/revision_002/alignment'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(exist_ok=False)
    started = time.time()
    original = {p.name:sha(p) for p in SOURCE.iterdir()}
    pred = pd.read_parquet(SOURCE/'predictions.parquet')
    gen = pd.read_parquet(SOURCE/'generator_vectors.parquet')
    replay = pd.read_parquet(SOURCE/'row_replay_seed41001.parquet')
    replay = replay[replay.propensity_representation.eq('lowdim_11') & replay.loss_representation.eq('lowdim_11') & replay.mechanism.eq('mar')]
    designs, summaries, replicates = [], [], []
    identity_max = 0.
    for endpoint_index, (endpoint, ep) in enumerate(pred.groupby('endpoint',sort=True)):
        n = ep.source_row_index.nunique()
        u = np.random.default_rng(73001+endpoint_index).random((10000,n))
        permutation = np.random.default_rng(72001+endpoint_index).permutation(n)
        for model, pp in ep.groupby('model',sort=True):
            pp = pp.sort_values('source_row_index')
            loss = pp.loss.to_numpy()
            truth = loss.mean()
            order = np.argsort(loss,kind='stable')
            for anchor, gg in gen[gen.endpoint.eq(endpoint) & gen.mechanism.eq('mar')].groupby('anchor',sort=True):
                source_p = gg.sort_values('source_row_index').p_oracle.to_numpy()
                fixed_m = replay[replay.endpoint.eq(endpoint) & replay.model.eq(model) & replay.anchor.eq(anchor)].sort_values('source_row_index').mhat.to_numpy()
                aligned, reversed_p = np.empty(n), np.empty(n)
                aligned[order], reversed_p[order] = np.sort(source_p), np.sort(source_p)[::-1]
                for pairing, p in [('original',source_p),('random',source_p[permutation]),('aligned',aligned),('reversed',reversed_p)]:
                    assert np.array_equal(np.sort(p),np.sort(source_p))
                    ident = {'endpoint':endpoint,'anchor':anchor,'model':model,'pairing':pairing}
                    designs.append({**ident,**probability_loss(loss,p),'probability_sd':p.std(),
                                    'sorted_probability_sha256':hashlib.sha256(np.sort(p).astype('<f8').tobytes()).hexdigest(),
                                    'permutation_seed':72001+endpoint_index,'mask_seed':73001+endpoint_index,'replicates':10000})
                    parts = []
                    for start in range(0,10000,500):
                        r = u[start:start+500]<p
                        k = r.sum(axis=1)
                        w = r/p
                        cc = r@loss/k
                        values = {'complete_case':cc,'oracle_ht':w@loss/n,'oracle_hajek':(w@loss)/w.sum(axis=1),
                                  'oracle_aipw_fixed_mean':truth+(w@(loss-truth))/n,
                                  'oracle_aipw_fixed_replay':(fixed_m+w*(loss-fixed_m)).mean(axis=1),
                                  'constant_q_true_mean':truth+(r@(loss-truth))/(n*p.mean())}
                        diff = values['constant_q_true_mean']-truth-k/(n*p.mean())*(cc-truth)
                        identity_max = max(identity_max,float(np.abs(diff).max()))
                        assert np.abs(diff).max()<1e-12
                        parts.append(pd.DataFrame({'replicate':np.arange(start,start+len(r)),**values}))
                    frame = pd.concat(parts,ignore_index=True)
                    cc_med = np.median(np.abs(frame.complete_case-truth))
                    for method in frame.columns[1:]:
                        stat = metrics(frame[method],truth)
                        summaries.append({**ident,'method':method,**stat,'median_abs_error_ratio_vs_cc':stat['median_abs_error']/cc_med})
                    replicates.append(frame.assign(**ident))
    pd.DataFrame(designs).to_csv(OUT/'design.csv',index=False,float_format='%.17g')
    pd.DataFrame(summaries).to_csv(OUT/'summary.csv',index=False,float_format='%.17g')
    pd.concat(replicates,ignore_index=True).to_parquet(OUT/'replicates.parquet',index=False)
    assert original=={p.name:sha(p) for p in SOURCE.iterdir()}
    record = {'status':'COMPLETED_POST_RESULT_FIXED_MARGINAL_DIAGNOSTIC','settings':len(designs),
              'replicate_rows':sum(len(x) for x in replicates),'summary_rows':len(summaries),
              'constant_q_identity_max_abs_difference':identity_max,'all_probability_marginals_exact':True,
              'frozen_files_unchanged':True,'duration_seconds':time.time()-started,
              'finished_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(Path(__file__)),
              'protocol_sha256':sha(OUT.parent/'protocol.md'),'input_sha256':original,
              'output_sha256':{p.name:sha(p) for p in OUT.iterdir() if p.is_file()}}
    (OUT/'manifest.json').write_text(json.dumps(record,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in record.items() if not k.endswith('sha256')},indent=2))


if __name__=='__main__':
    with threadpool_limits(limits=1):
        main()
