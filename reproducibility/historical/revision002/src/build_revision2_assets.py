from pathlib import Path
import json
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import build_paper_assets as old

ROOT = Path(__file__).resolve().parents[1]
REV = ROOT/'artifacts/revision_002'
PAPER = ROOT/'paper'
KEYS = ['endpoint','anchor','mechanism','model']
LABEL = {'LOG_HLM':'HLM','LOG_RLM':'RLM','LOG_SOL':'SOL','LOG_MDR1_ER':'MDR1','ecfp_rf':'RF','ecfp_desc_xgb':'XGB'}
COLORS = ['#557F83','#A16D68','#8B8972','#7686A0','#777777']


def save_csv(frame,name):
    frame.to_csv(PAPER/'tables'/name,index=False,float_format='%.17g')


def rr(s,d=3):
    return f'{s.min():.{d}f} to {s.max():.{d}f}'


def export(fig,name):
    fig.savefig(PAPER/'figures'/f'{name}.png',dpi=260,bbox_inches='tight')
    fig.savefig(PAPER/'figures'/f'{name}.pdf',bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None})
    fig.savefig(PAPER/'figures'/f'{name}.svg',bbox_inches='tight',metadata={'Date':None})
    plt.close(fig)


def main():
    old.style()
    plt.rcParams.update({'svg.hashsalt':'jpa-revision002','axes.spines.top':False,'axes.spines.right':False})
    original = pd.read_csv(REV/'original/estimates.csv',float_precision='round_trip')
    state = original[original.method.eq('aipw')].groupby(['propensity','loss_variant','case_R']).agg(
        estimates=('estimate','size'),median_abs_error=('abs_error','median'),median_case_share=('case_abs_share','median'),
        min_case_share=('case_abs_share','min'),max_case_share=('case_abs_share','max'),
        top_case_count=('largest_source_row',lambda x:int((x==2945).sum())),median_coefficient=('case_m_coefficient','median')).reset_index()
    save_csv(state,'Table_S13_state_summary.csv')
    display = state[['propensity','loss_variant','case_R','estimates','median_abs_error','median_case_share','top_case_count']].copy()
    display.columns=['P','Loss input','R','Estimates','Median error','Case share','Case top count']
    display['Median error']=display['Median error'].map(lambda x:f'{x:.6g}')
    display['Case share']=display['Case share'].map(lambda x:f'{100*x:.3f}%')
    save_csv(display,'Table_S13_display.csv')
    save_csv(original.groupby(KEYS+['propensity','loss_variant','method','case_R'],dropna=False).agg(n=('estimate','size'),
        median_abs_error=('abs_error','median'),case_share_median=('case_abs_share','median')).reset_index(),'Table_S13_setting_states.csv')
    b = pd.read_csv(REV/'bridge/summary.csv',float_precision='round_trip')
    cc = b[b.method.eq('complete_case')].set_index(KEYS).median_abs_error
    b['median_error_ratio_vs_cc']=[r.median_abs_error/cc.loc[tuple(getattr(r,k) for k in KEYS)] for r in b.itertuples()]
    save_csv(b,'Table_S14_all_fitted_comparisons.csv')
    bd = pd.read_csv(REV/'bridge/design.csv',float_precision='round_trip')
    save_csv(bd,'Table_S14_selection_design.csv')
    rows=[]
    for key,g in b.groupby(KEYS,sort=True):
        ep,anchor,mech,model=key
        row={'Setting':f'{LABEL[ep]} {LABEL[anchor]} {LABEL[model]}','MW sign':'+' if mech.endswith('positive') else '-'}
        for label,method,p,m in [('CC','complete_case','none','none'),('Low M','aipw','low','low'),('Drop Ipc','aipw','low','drop'),
                                ('Oracle low','aipw','oracle','low'),('Oracle drop','aipw','oracle','drop'),('Oracle HT','ht','oracle','none'),('Oracle Hajek','hajek','oracle','none')]:
            value=g[g.method.eq(method)&g.propensity.eq(p)&g.loss_variant.eq(m)].iloc[0]
            row[label]=f'{value.median_abs_error:.5f}'
        rows.append(row)
    save_csv(pd.DataFrame(rows),'Table_S14_display.csv')
    al = pd.read_csv(REV/'alignment/summary.csv',float_precision='round_trip')
    design = pd.read_csv(REV/'alignment/design.csv',float_precision='round_trip')
    save_csv(al,'Table_S15_alignment_all.csv')
    save_csv(design,'Table_S15_alignment_design.csv')
    rows=[]
    for (pairing,method),g in al.groupby(['pairing','method'],sort=True):
        rows.append({'Pairing':pairing,'Reference':method.replace('oracle_','').replace('aipw_fixed_','AIPW ').replace('_',' '),
            'Settings':len(g),'Bias range':rr(g.signed_bias,4),'RMSE range':rr(g.rmse,4),'Error ratio range':rr(g.median_abs_error_ratio_vs_cc)})
    save_csv(pd.DataFrame(rows),'Table_S15_display.csv')
    solver=pd.read_csv(REV/'independent_solver.csv')
    save_csv(solver,'Table_S16_independent_solver.csv')
    low=b[b.method.eq('aipw')&b.propensity.eq('low')&b.loss_variant.eq('low')]
    drop=b[b.method.eq('aipw')&b.propensity.eq('low')&b.loss_variant.eq('drop')]
    evidence={'bridge_low_ratio_range':[low.median_error_ratio_vs_cc.min(),low.median_error_ratio_vs_cc.max()],
        'bridge_drop_ratio_range':[drop.median_error_ratio_vs_cc.min(),drop.median_error_ratio_vs_cc.max()],
        'bridge_low_better_cc':int((low.median_error_ratio_vs_cc<1).sum()),'bridge_drop_better_cc':int((drop.median_error_ratio_vs_cc<1).sum()),
        'bridge_low_25percent':int((low.median_error_ratio_vs_cc<=.75).sum()),'bridge_drop_25percent':int((drop.median_error_ratio_vs_cc<=.75).sum()),
        'drop_better_low':int((drop.set_index(KEYS).median_abs_error<low.set_index(KEYS).median_abs_error).sum()),
        'original_lowp_pooled_medians':original[original.method.eq('aipw')&original.propensity.eq('low')].groupby('loss_variant').abs_error.median().to_dict(),
        'bridge_probability_loss_correlation_range':[bd.p_loss_correlation.min(),bd.p_loss_correlation.max()],
        'bridge_probability_weighted_shift_range':[bd.probability_weighted_shift.min(),bd.probability_weighted_shift.max()]}
    (REV/'paper_numbers.json').write_text(json.dumps(evidence,indent=2)+'\n')
    save_csv(pd.DataFrame([
        ['Original MAR','8 cells × 100 masks','CC signed bias','-0.002879 to 0.003924','Fixed measured cohorts'],
        ['Observable MW selection','16 cells × 100 refits','Low M / CC median-error ratio',rr(low.median_error_ratio_vs_cc),'Actual nuisance refits'],
        ['Observable MW selection','16 cells × 100 refits','Drop Ipc / CC median-error ratio',rr(drop.median_error_ratio_vs_cc),'Same low-dimensional P'],
        ['Raw full loss nuisance','1,600 estimates per P','Known case ranks first','1,600 / 1,600 for all three P','100 original masks per cell'],
        ['Fixed marginal alignment','8 cells × 10,000 masks per pairing','Oracle Hajek / CC ratio',rr(al[al.method.eq('oracle_hajek')&al.pairing.isin(['aligned','reversed'])].median_abs_error_ratio_vs_cc),'Idealized loss-informed pairing'],
    ],columns=['Analysis','Denominator','Quantity','Result','Scope']),'Table_2_quantitative_evidence.csv')
    original.to_csv(PAPER/'source_data/figure3_all_seed_estimates.csv',index=False,float_format='%.17g')
    b.to_csv(PAPER/'source_data/figure4_fitted_bridge.csv',index=False,float_format='%.17g')
    al.to_csv(PAPER/'source_data/figure2_fixed_marginal.csv',index=False,float_format='%.17g')
    fig,ax=plt.subplots(figsize=(7.2,3.1)); ax.axis('off')
    blocks=[(.01,.62,.29,.31,'Known finite-cohort loss','HLM 929 and RLM 922\nFixed RF and XGB predictors\nPreviously measured labels only'),
            (.36,.62,.28,.31,'Selective references','Original MCAR / MAR\nMolecular-weight selection\nProbability alignment controls'),
            (.70,.62,.29,.31,'Risk estimation','Complete case, HT and Hajek\nAIPW with fitted or oracle p\nLow, raw and processed loss')]
    from matplotlib.patches import Rectangle
    for x,y,w,h,title,body in blocks:
        ax.add_patch(Rectangle((x,y),w,h,transform=ax.transAxes,facecolor='#F4F5F4',edgecolor='#B5BEBC',lw=.7))
        ax.text(x+.015,y+h-.035,title,weight='bold',va='top',fontsize=8.5)
        ax.text(x+.015,y+h-.12,body,va='top',fontsize=7.2,linespacing=1.6)
    for a,z in [(.305,.355),(.645,.695)]:
        ax.annotate('',xy=(z,.775),xytext=(a,.775),arrowprops={'arrowstyle':'->','color':'#777777'})
    ax.text(.5,.47,r'$\psi_i=(1-R_i/\hat p_i)\hat m_i+R_iL_i/\hat p_i$',ha='center',fontsize=14)
    ax.text(.23,.24,'Reference absent  R = 0\nPredicted loss enters with coefficient 1',ha='center',fontsize=9,linespacing=1.6)
    ax.text(.75,.24,'Reference observed  R = 1\nCoefficient 1 - 1/p need not be zero',ha='center',fontsize=9,linespacing=1.6)
    ax.text(.5,.04,'Separate absolute accuracy, benefit over complete cases, and the uncertainty target',ha='center',fontsize=8,color='#5E6368')
    export(fig,'Figure_1_risk_mechanism')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.3)); fig.subplots_adjust(bottom=.26,wspace=.42,top=.86)
    pairings=['original','random','aligned','reversed']
    for i,pairing in enumerate(pairings):
        g=design[design.pairing.eq(pairing)]
        axes[0].scatter(i+np.linspace(-.13,.13,8),g.probability_weighted_shift,s=19,color=COLORS[0])
    axes[0].axhline(0,color='#999999',lw=.7)
    axes[0].set(xticks=range(4),xticklabels=['Original','Random','Aligned','Reversed'],ylabel='Probability-weighted risk shift')
    axes[0].set_title('a  Identical probability marginals',loc='left',weight='bold')
    methods=['oracle_hajek','oracle_ht','oracle_aipw_fixed_replay']
    for j,(method,label) in enumerate(zip(methods,['Hajek','HT','AIPW fixed replay'])):
        for i,pairing in enumerate(pairings):
            values=al[al.method.eq(method)&al.pairing.eq(pairing)].median_abs_error_ratio_vs_cc
            axes[1].scatter(i+(j-1)*.2+np.linspace(-.045,.045,8),values,s=14,color=COLORS[j],marker=['o','^','s'][j],label=label if i==0 else None)
    axes[1].axhline(1,color='#777777',lw=.8)
    axes[1].axhline(.75,color='#999999',lw=.7,ls='--')
    axes[1].set(xticks=range(4),xticklabels=['Original','Random','Aligned','Reversed'],ylabel='Median absolute-error ratio vs CC',ylim=(0,2.3))
    axes[1].set_title('b  Oracle references and comparators',loc='left',weight='bold')
    axes[1].legend(loc='upper center',bbox_to_anchor=(.45,-.20),ncol=3,fontsize=6.8,handletextpad=.1,columnspacing=.8)
    export(fig,'Figure_2_selection_alignment')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.3)); fig.subplots_adjust(bottom=.23,wspace=.44,top=.85)
    raw=original[original.method.eq('aipw')&original.propensity.eq('low')&original.loss_variant.eq('raw')]
    axes[0].scatter(raw.case_m_coefficient,raw.case_pseudo_contribution,s=6,alpha=.32,color=COLORS[1],rasterized=True)
    axes[0].set_yscale('symlog',linthresh=.1)
    axes[0].set(xlabel='Case coefficient 1 - R/p',ylabel='Case contribution to estimated risk')
    axes[0].set_title('a  Both observation states propagate',loc='left',weight='bold')
    axes[0].text(.06,.49,'Low P; 1,600 estimates\nR = 0: 820; R = 1: 780',transform=axes[0].transAxes,fontsize=6.6)
    for j,state_value in enumerate([0,1]):
        for i,variant in enumerate(['raw','heldout_clip','drop','low']):
            g=original[original.method.eq('aipw')&original.propensity.eq('low')&original.loss_variant.eq(variant)&original.case_R.eq(state_value)].abs_error
            lo,md,hi=np.quantile(g,[.1,.5,.9])
            axes[1].errorbar(i+(j-.5)*.17,md,yerr=[[md-lo],[hi-md]],color=COLORS[j],marker=['o','s'][j],capsize=2,ms=4,label=f'R = {state_value}' if i==0 else None)
    axes[1].set_yscale('log')
    axes[1].set(xticks=range(4),xticklabels=['Raw','Held-out\nclip','Drop Ipc\nrefit','Low M'],ylabel='Absolute risk error',ylim=(.0003,10000))
    axes[1].set_title('b  Processing versus sensitivity',loc='left',weight='bold')
    axes[1].legend(loc='upper right',fontsize=7)
    export(fig,'Figure_4_ipc_propagation')
    fig,axes=plt.subplots(1,2,figsize=(7.2,4.3)); fig.subplots_adjust(left=.15,right=.98,bottom=.18,top=.88,wspace=.35)
    chosen=[('low','low','Fitted low M','o'),('low','drop','Fitted drop Ipc','s'),('oracle','low','Oracle low M','^'),('oracle','drop','Oracle drop Ipc','d')]
    order=list(b.groupby(KEYS,sort=True).groups)
    labels=[f'{LABEL[e]} {LABEL[a]} {LABEL[m]} '+('+' if mech.endswith('positive') else '-') for e,a,mech,m in order]
    for j,(p,m,label,marker) in enumerate(chosen):
        g=b[b.method.eq('aipw')&b.propensity.eq(p)&b.loss_variant.eq(m)].set_index(KEYS).loc[order]
        y=np.arange(16)+(j-1.5)*.14
        axes[0].scatter(g.signed_bias,y,s=17,color=COLORS[j],marker=marker,label=label)
        axes[1].scatter(g.median_error_ratio_vs_cc,y,s=17,color=COLORS[j],marker=marker)
    axes[0].axvline(0,color='#999999',lw=.7)
    axes[1].axvline(1,color='#777777',lw=.8)
    axes[1].axvline(.75,color='#999999',lw=.7,ls='--')
    for ax in axes:
        ax.set(yticks=range(16),ylim=(15.7,-.7)); ax.grid(axis='y',color='#ECECEC',lw=.4)
    axes[0].set(yticklabels=labels,xlabel='Monte Carlo signed bias')
    axes[1].set(yticklabels=[],xlabel='Median absolute-error ratio vs CC')
    axes[0].set_title('a  Actual refitted auxiliaries',loc='left',weight='bold')
    axes[1].set_title('b  Benefit varies across settings',loc='left',weight='bold')
    fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.52,.01),ncol=4,fontsize=7)
    export(fig,'Figure_3_fitted_bridge')
    fig = old.graphical_abstract()
    updates = {
        '1  Measure selection bias':'1  Hold probability marginals',
        'Original MAR: weak bias':'Same probabilities; changed pairing',
        '|Corr(probability, loss)| < 0.107':'Hajek / CC: 0.883–1.104 original',
        'Oracle also misses 25% benefit':'0.113–0.168 aligned / reversed',
        'Benefit is distinct from accuracy':'Loss-informed oracle diagnostic',
        '2  Trace extrapolation':'2  Trace all original seeds',
        'Single-seed Ipc diagnostic':'100 seeds; both observation states',
        'Limit only held-out Ipc range':'Raw Ipc dropped before full refit',
        'Low-P median error: 702.65 → 0.00595':'Low-P median error: 323.07 → 0.00552',
        'Configuration-specific finding':'Fixed drug-property models',
        '3  Audit PXR records':'3  Refit under MW selection',
        'Dose-table membership':'Observable molecular variables only',
        'Hits: 869/873':'Low M beats CC: 11/16 settings',
        'Non-hits: 1,853/9,874':'Drop Ipc beats CC: 10/16 settings',
        'Full-library risk not estimated':'Drop beats low M in only 1/16',
        'Original frozen results retained; bias, interval and extrapolation diagnostics are post-result':'Original frozen results retained; new refits and alignment controls are post-result',
    }
    for artist in fig.axes[0].texts:
        if artist.get_text() in updates:
            artist.set_text(updates[artist.get_text()])
            if artist.get_text().startswith(('1  ','2  ','3  ')):
                artist.set_fontsize(9)
    old.FIGURES = PAPER/'figures'
    old.save_figure(fig,'Graphical_Abstract')
    for source,target in [('Figure_2_confirmatory_signatures','Figure_S2_original_signatures'),('Figure_3_admission_translation_pxr','Figure_S3_original_selection'),('Figure_4_diagnostic_mechanisms','Figure_S4_conditional_diagnostics')]:
        for ext in ['png','pdf','svg']:
            a, b = PAPER/'figures'/f'{source}.{ext}', PAPER/'figures'/f'{target}.{ext}'
            if a.exists():
                shutil.copy2(a,b)
            else:
                assert b.is_file(), 'Historical supplementary figure must be supplied for aggregate-only replay'
    print(json.dumps(evidence,indent=2))


if __name__=='__main__':
    main()
