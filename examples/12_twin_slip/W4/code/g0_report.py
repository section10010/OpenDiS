"""Chinese G0 review and diagnostic figures from actual native saved outputs."""
import argparse,html,json,csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from g0_analyze import read
from g0_campaign import FAMILIES,MASKS


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--analysis',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.root;o=a.out;o.mkdir(parents=True,exist_ok=True)
    data=read(a.analysis/'case_metrics.json');by={x['case_id']:x for x in data};review=read(a.analysis/'gate_review.json');refine=read(a.analysis/'primary_refinement.json');resolution=read(a.analysis/'contrast_resolution.json');contrasts=read(a.analysis/'four_wall_contrasts.json')
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    colors=dict(none='#647686',lower='#2675a9',upper='#ba8738',both='#b33c43')
    def save(fig,name):
        fig.savefig(o/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(o/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10,7))
    for col,(kind,mode) in enumerate(FAMILIES):
        for mask in MASKS:
            for label,ls in [('cap4000','-'),('uncapped_diagnostic','--')]:
                folder=r/'cases'/f'G0_{kind}_{mode}_{mask}_{label}'
                with (folder/'history.csv').open() as f:h=list(csv.DictReader(f))
                axes[0,col].plot([float(x['time_ns']) for x in h],[float(x['swept_area_nm2'])/1e6 for x in h],ls,color=colors[mask],lw=1.,label=mask+' / '+('cap' if ls=='-' else 'uncapped'))
        cap=[by[f'G0_{kind}_{mode}_{mask}_cap4000'] for mask in MASKS];x=np.arange(4)
        axes[1,col].bar(x-.18,[z['cap_line_time_fraction']*100 for z in cap],.36,color='#2675a9',label='Line x accepted-time fraction')
        axes[1,col].bar(x+.18,[z['cap_any_record_fraction']*100 for z in cap],.36,color='#b7bfc7',label='Any capped node / record count')
        axes[1,col].set(xticks=x,xticklabels=MASKS,ylim=(0,105),ylabel='Fraction (%)');axes[0,col].set(title=kind+' / '+mode,xlabel='Actual time (ns)',ylabel='Integration sweep (um2)');axes[0,col].grid(alpha=.16)
    axes[0,1].legend(frameon=False,fontsize=6,ncol=2);axes[1,0].legend(frameon=False,fontsize=7);fig.suptitle('G0 real native diagnostics: capped vs uncapped; uncapped is not material calibration');fig.tight_layout(rect=(0,0,1,.96));save(fig,'G0_cap_diagnostics')
    labels=['space_20_to_10','time_half_at_10','space_10_to_5','time_half_at_5'];short=['20 -> 10 nm','dt/2 at 10','10 -> 5 nm','dt/2 at 5']
    fig,axes=plt.subplots(2,3,figsize=(12,6.7))
    for row,(kind,mode) in enumerate(FAMILIES):
        for col,metric in enumerate(('area_nm2','line_nm','tail_nm')):
            ax=axes[row,col]
            for j,mask in enumerate(MASKS):
                xx=[];yy=[]
                for x in refine:
                    if x['kind']==kind and x['mode']==mode and x['wall']==mask and x['relative_errors'][metric] is not None:
                        xx.append(labels.index(x['comparison'])+(j-1.5)*.035);yy.append(100*x['relative_errors'][metric])
                ax.plot(xx,yy,'o-',color=colors[mask],ms=3,lw=.7,label=mask)
            ax.axhline(2,color='#aa3c41',ls='--',lw=.9);ax.set(xticks=range(4),xticklabels=short,ylabel='Relative difference (%)',title=kind+' / '+mode+' / '+metric);ax.tick_params(axis='x',labelrotation=20);ax.grid(alpha=.14)
    axes[0,2].legend(frameon=False,fontsize=7);fig.suptitle('All measured refinement errors retained; missing tests are not plotted as zero');fig.tight_layout(rect=(0,0,1,.96));save(fig,'G0_refinement')
    fig,axes=plt.subplots(2,2,figsize=(10,6.5))
    for row,(kind,mode) in enumerate(FAMILIES):
        res=next((x for x in resolution if x['kind']==kind and x['mode']==mode),None)
        labels2=['cap4000','uncapped_diagnostic','tau250_uncapped']+([res['finest_label']] if res else [])
        for col,key in enumerate(('second_upper_nm2','nonadditivity_nm2')):
            vals=[]
            for label in labels2:
                x=next(z for z in contrasts if z['kind']==kind and z['mode']==mode and z['label']==label);vals.append(x[key]/1e6)
            axes[row,col].plot(range(len(vals)),vals,'o',color='#2675a9')
            if res:
                bound=res['second_refinement_error_sum_nm2'] if col==0 else res['observed_refinement_error_sum_nm2']
                axes[row,col].errorbar([len(vals)-1],[vals[-1]],yerr=[bound/1e6],fmt='none',color='#ad3f43',capsize=4)
            axes[row,col].axhline(0,color='#8a959e',lw=.6);axes[row,col].set(xticks=range(len(vals)),xticklabels=labels2,ylabel='Area contrast (um2)',title=kind+' / '+mode+' / '+('second wall' if col==0 else 'nonadditivity'));axes[row,col].tick_params(axis='x',labelrotation=20);axes[row,col].grid(alpha=.15)
    fig.suptitle('Native contrasts; finest red bars = observed numerical-refinement estimate, not seed CI');fig.tight_layout(rect=(0,0,1,.95));save(fig,'G0_contrasts')
    md=[];web=[]
    def h(s,l=2):md.append('#'*l+' '+s+'\n');web.append(f'<h{l}>{html.escape(s)}</h{l}>')
    def text(s):md.append(s+'\n');web.append('<p>'+html.escape(s)+'</p>')
    def table(headers,rows):
        md.append('| '+' | '.join(headers)+' |\n| '+' | '.join(['---']*len(headers))+' |\n'+'\n'.join('| '+' | '.join(map(str,x))+' |' for x in rows)+'\n')
        web.append('<table><tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in headers)+'</tr>'+''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in x)+'</tr>' for x in rows)+'</table>')
    def pic(n,s):md.append(f'![{s}]({n}.png)\n');web.append(f'<figure><img src="{n}.png"><figcaption>{html.escape(s)}</figcaption></figure>')
    h('W4-G0 中文实施与独立复核报告',1)
    text('2026-10-30；有效FCC、非周期单通道I0基准。只审查真实ExaDiS轨迹，没有生成新论文F1–F10图，没有开展四厚度集体计算。W1–W3及规划原输入未改写。')
    table(['审查项','实际裁决'],[[k,str(review[k])] for k in ('status','screen_completed','case_counts','native_numerical_invariants','primary_refinement','dynamic_release_regression','second_wall_mechanism','actual_time_cap_insensitivity','VCoNi','G1_admission')])
    h('1. 实际入口与源码变更')
    text('已现场核实VM仓库 /home/sakura/opendis/src、构建 /home/sakura/opendis/build、Python /opt/miniforge3/envs/anvil/bin/python。W3历史快照19个Python文件与生产源码逐字节相同；原生W1 pin补丁仍为14增3删，没有更改上游C++。新代码只在W4副本，新增诊断、批次、细化、审查模块。')
    text('W3旧副本和W4修改后副本分别重新完成A0：5项配置、10项原生响应（含3项pin契约）、12项上游原生单元和2项FR流程回归通过。W4配置/谱系/诊断8项通过，参考测试46项在VM通过；Windows缺pytest的首次记录保留，不能写成Windows参考测试PASS。新增原生两模式反力重建、三速度区分及6项接触验收通过。真实正例、非法初始接触FAIL和资源CENSORED负控也通过。')
    text('逐节点保存原生GLIDE未约束速度、接触后速度和最终限速速度；反力来自同一个原生迁移率MN，测试以M(f+反力)独立重建接触后速度。每个接受步记录其真实dt、tag、坐标和半连接线长权重；不使用稀疏图帧插值生成动力学。')
    h('2. 24例原始筛查与限速结论')
    text('24例来自原规划两构型×四壁×三个限速/载荷变体，逐例原样保留输入。另做空间、时间、容差、碰撞半径、拓扑、重网格和迁移率时钟对照。每个四壁族初态数组逐项相同，其余物理与数值参数一致；只允许case名称和壁面开关不同。')
    caprows=[x for x in data if x['group']=='screen' and x['case_id'].endswith('cap4000')]
    table(['构型/壁面','真实线长-时间限速占比','至少一节点限速的记录比例','原始最大速度m/s','接触后最大速度m/s'],[[x['kind']+'/'+x['wall'],f'{100*x["cap_line_time_fraction"]:.3f}%',f'{100*x["cap_any_record_fraction"]:.3f}%',f'{x["max_raw_velocity_m_s"]:.1f}',f'{x["max_constrained_velocity_m_s"]:.1f}'] for x in caprows])
    text('两个比例的分母不同，不能混称。loop/HI双壁原始速度超过4000 m/s，但接触后未超过上限，限速实际占比为零。反之threader/HII双壁的记录比例远大于真实线长-时间占比。原始速度高可能主要是被法向反力抵消的方向，不能直接作为实际传播速度。')
    cap=read(a.analysis/'cap_sensitivity.json');table(['构型/壁面','取消4000上限后的面积相对变化','固定真实时间的限速不敏感资格'],[[x['kind']+'/'+x['wall'],f'{x["relative_area_change"]*100:.3f}%',x['timing_cap_insensitive']] for x in cap])
    text('取消限速是数学诊断，不意味着更真实的VCoNi。高速度线性过阻尼/准静态弹性适用性未被标定；不能把本次时间轴映射到实际准静态材料。M/2并将时间和最大dt加倍的8个对照检验几何时钟关系，不冒充新的材料机制。')
    pic('G0_cap_diagnostics','真实接受步曲线与两种限速比例。')
    h('3. 数值细化及失败保留')
    failed=[x for x in refine if x['status']!='PASS']
    table(['构型/壁面','细化','面积差','线长差','尾段差','状态'],[[x['kind']+'/'+x['wall'],x['comparison']]+[f'{100*x["relative_errors"][k]:.3f}%' if x['relative_errors'][k] is not None else 'UNRESOLVED' for k in ('area_nm2','line_nm','tail_nm')]+[x['status']] for x in failed])
    text('原准则主要量最后两级≤2%、事件时刻≤5%没有放宽。所有不满足的粗级结果仍保留，不把后续细化改写为历史PASS。5/2 nm及独立时间再减半只针对真实未收敛构型，扩展在运行前单独冻结，预算1800 s；若仍不足则保留UNRESOLVED/CENSORED，不无限运行。')
    pic('G0_refinement','全部实际细化误差；缺失的试验不填零，2%线是事前门槛。')
    h('4. 第二界面和非加性')
    table(['构型','变体','第二壁增量µm²','非加性µm²','相对单下壁面积减少'],[[x['kind'],x['label'],f'{x["second_upper_nm2"]/1e6:.6f}',f'{x["nonadditivity_nm2"]/1e6:.6f}',f'{100*x["reduction_vs_lower"]:.2f}%'] for x in contrasts if x['label'] in ('cap4000','uncapped_diagnostic','tau250_uncapped','mesh10_timehalf','mesh5_timequarter')])
    table(['构型','最细非加性µm²','四项细化差之和µm²','非加性数值符号','第二壁可辨认'],[[x['kind'],f'{x["finest_nonadditivity_nm2"]/1e6:.6f}',f'{x["observed_refinement_error_sum_nm2"]/1e6:.6f}',x['nonadditivity_numerical_sign'],str(x['second_wall_detectable'])] for x in resolution])
    text('4000 m/s限速与取消限速之间，小非加性发生符号变化；不能据W3负的微小差值声称稳定的协同强化。误差条来自已观察的空间/时间细化差，既不是严格数学上界，也不是独立种子置信区间。第二壁的大增量与小非加性必须分别判读；没有分辨非加性不代表不存在双壁约束。')
    pic('G0_contrasts','只显示真实四壁差；最细级附经验数值误差，不制造预期符号。')
    h('5. 守恒、能量及未完成能力')
    good=[x for x in data if x['status']=='PASS']
    text(f'已完成案例最大Burgers残差为{max(x["max_burgers_residual"] for x in good):.3e} b，最大归一化接触功残差{max(x["max_contact_power_residual"] for x in good):.3e}，最大全局塑性闭合残差{max(x["plastic_closure"] for x in good):.3e}。每次积分、碰撞、拓扑、Remesh输出均检查无非法穿透，事件账本区分几何接触成员与Remesh换tag。')
    text('所有接受步的完整三速度档案重新审查限速重建、实际dt、线长权重与正乘子；额外保留各原生操作的全局面积张量增量、接触功、非负耗散和由原生塑性计算的外功。尚无独立弹性势能API，未将由力积分的功改名为独立能量守恒；绝对弹性势能为NOT_RUN。局域拓扑塑性张量场、I1反应、多层接触与VCoNi资格也没有由这些记录自动通过。')
    if (r/'dynamic_release_review_r2.json').exists():
        rel=read(r/'dynamic_release_review_r2.json')
        table(['释放构型','实际终态','湮灭前无接触状态数','首次证实释放/ns','末态塑性闭合'],[[x['kind']+'/'+x['mode'],x['terminal_state'],x['released_nonempty_states'],f'{x["first_released_time_ns"]:.6f}',f'{x["final_native_plastic_closure"]:.3e}'] for x in rel['checks']])
        text('取消限速的反向加载出现合法湮灭终态，首次审查脚本对空数组求最小值的错误已保留。修复仅涉及审查与空面积处理：不重跑已存在HI轨迹，不改原result。通过湮灭前的真实非空释放状态和末态原生塑性闭合验收；空网络间隙保留null，不填写无穷大或伪零。')
    h('6. 放行边界与恢复')
    text('本次G0判定为 '+review['status']+'，G1准入为 '+review['G1_admission']+'。即便第二壁因果信号已可辨认，整体数值门禁未过时仍不得进入G1。若G0通过，下一步也只允许共享CTB小基准，不允许直接集体RVE。')
    text('资源截断保存checkpoint.npz及累计塑性/扫掠/时间，状态CENSORED；精确恢复完整累计账本尚未验收，不把分段重启直接拼接成无间断长轨迹。历史输出不可覆盖；失败保存first-invalid与last-legal，旧日志、首次Windows文件句柄错误及修复后测试均保留。')
    text('原始数据位于VM /home/sakura/opendis/w4_g0_20261030_001/production。每例有配置、精确命令、代码快照、线程/种子、退出码、初末构型、逐步history、接触事件、三速度压缩块；复跑使用新case ID。分析脚本为W4/code/g0_analyze.py，图/报告脚本为g0_report.py。长期状态入口为DDD_Long_Term_Plan/PROJECT_STATUS.md和NEXT_ACTION.md。')
    (o/'W4_G0_REPORT_ZH.md').write_text('\n'.join(md),encoding='utf-8')
    css='body{max-width:1180px;margin:36px auto;padding:0 24px;color:#26394a;font:16px/1.75 "Microsoft YaHei",sans-serif}h1,h2{color:#1d405b}h2{border-bottom:1px solid #ccd5df;padding-top:18px}table{border-collapse:collapse;width:100%;font-size:13px;margin:20px 0}td,th{border:1px solid #d0dae2;padding:8px;text-align:left}th{background:#edf2f6}img{max-width:100%}figure{margin:26px 0}figcaption{font-size:14px;color:#526474}'
    (o/'W4_G0_REPORT_ZH.html').write_text('<!doctype html><html lang="zh-CN"><meta charset="UTF-8"><title>W4-G0 审查报告</title><style>'+css+'</style><body>'+''.join(web)+'</body></html>',encoding='utf-8')
    environment=read(r/'environment/runtime.json') if (r/'environment/runtime.json').exists() else {'status':'NOT_RUN'}
    (o/'figure_manifest.json').write_text(json.dumps(dict(status='REAL_NATIVE_DIAGNOSTICS',paper_figure_ids=[],
        figures=['G0_cap_diagnostics','G0_refinement','G0_contrasts'],case_ids=[x['case_id'] for x in data],
        raw_root=str(r),analysis=str(a.analysis),analysis_script='code/g0_analyze.py',figure_script='code/g0_report.py',
        native_environment=environment,gate_review=review,
        selection_rule='All predefined four-wall groups and all measured refinement errors; common actual endpoints, no favourable-frame selection',
        statistical_unit='One fixed seed 42 per family; parameter diagnostics are not independent microstructure replicates',
        units=dict(time='ns',area='um2',velocity='m/s',cap_fraction='percent'),material_qualification='NOT_QUALIFIED'),indent=2),encoding='utf-8')
    print('G0_REPORT_SAVED',review['status'])


if __name__=='__main__':main()
