"""Final G1 closeout from existing evidence, including incomplete fine jobs.

Never launches ExaDiS, fills missing observations with zero, changes thresholds,
or overwrites an earlier review. Numerical completion and qualification differ.
"""
import argparse,csv,html,json,re
from pathlib import Path
import numpy as np


def load(path):return json.loads(Path(path).read_text())


def raw_audit(folder):
    """Recompute graph/tail/plastic witnesses from real saved arrays."""
    c=load(folder/'config.json');source=load(folder/'spec.json');tol=c['contact_tol_nm'];z0=source['interface']['z_nm'];tails={}
    with np.load(folder/'initial.npz') as ini,np.load(folder/'final.npz') as fin:
        ij=fin['ids'];x=fin['positions_nm'];b=fin['burgers'];r=np.zeros_like(x)
        np.add.at(r,ij[:,0],b);np.add.at(r,ij[:,1],-b);burg=float(np.linalg.norm(r,axis=1).max(initial=0));assert burg<=1e-8
        np.testing.assert_allclose(np.einsum('ij,ij->i',b,fin['planes']),0,atol=1e-8)
        lengths=np.linalg.norm(x[ij[:,1]]-x[ij[:,0]],axis=1)
        for host in ('M','T'):
            node=fin['node_host']==host;edge=fin['segment_host']==host
            if not node.any():continue
            gap=(-1 if host=='M' else 1)*(x[:,2]-z0);assert gap[node].min()>=-tol
            active=gap<=tol;tails[host]=float(lengths[edge&np.all(active[ij],axis=1)].sum())
        def area(z):
            xx=z['positions_nm']-np.array([2000.,2000.,2000.]);ii=z['ids']
            return np.einsum('si,sj->ij',z['burgers'],.5*np.cross(xx[ii[:,0]],xx[ii[:,1]]))
        beta=-c['b_m']*1e9*(area(fin)-area(ini))/np.prod(c['box_nm']);native=np.array(load(folder/'plastic.json')['native_cumulative'])
        closure=float(np.linalg.norm(.5*(beta+beta.T)-native));assert closure<1e-12
    count=0;dt=0.;max_power=0.;peak=0.
    for file in sorted((folder/'accepted_velocity').glob('chunk_*.npz')):
        with np.load(file) as z:
            count+=len(z['dt_s']);dt+=float(z['dt_s'].sum());assert np.all(z['dt_s']>0)
            for key in ('positions_nm','force_N','raw_velocity_m_s','constrained_velocity_m_s','capped_velocity_m_s','reaction_N'):assert np.isfinite(z[key]).all()
            np.testing.assert_array_equal(z['capped_velocity_m_s'],z['constrained_velocity_m_s'])
            sign=np.where(z['node_host']=='M',-1.,1.);gap=sign*(z['positions_nm'][:,2]-z0);assert gap.min(initial=0)>=-tol
            np.testing.assert_allclose(z['reaction_N'][:,2],sign*z['multipliers_N'][:,0],rtol=1e-12,atol=1e-25)
            assert z['multipliers_N'].min(initial=0)>=0
            vc=z['constrained_velocity_m_s'];vr=z['raw_velocity_m_s'];reaction=z['reaction_N']
            speed=np.maximum(np.linalg.norm(vr,axis=1),np.linalg.norm(vc,axis=1));peak=max(peak,float(np.linalg.norm(vr,axis=1).max(initial=0)))
            scale=np.maximum(np.linalg.norm(reaction,axis=1),np.linalg.norm(z['force_N'],axis=1))*speed
            power=np.abs(np.einsum('ij,ij->i',reaction,vc))/np.maximum(scale,1e-30);max_power=max(max_power,float(power.max(initial=0)))
    result=load(folder/'result.json');assert count==result['steps'] and abs(dt-result['time_ns']*1e-9)<1e-22 and max_power<=1e-8
    return dict(status='PASS',native_tail_nm=tails,accepted_samples=count,accepted_time_ns=dt*1e9,max_burgers_residual_b=burg,
        final_native_plastic_closure=closure,max_contact_power_normalized=max_power,maximum_raw_velocity_m_s=peak)


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();root=a.root;out=a.out
    if out.exists():raise RuntimeError('Refuse overwrite')
    out.mkdir(parents=True);base=root/'analysis_r1';gate=load(base/'gate_review.json');manifest=load(root/'spatial_extension_manifest.json')
    criteria=load(root/'case_specs/base_MT.json')['criteria'];limit=criteria['primary_refinement_relative'];event_limit=criteria['event_relative']
    main_jobs=load(root/'campaign_manifest.json')['jobs'];assert len(main_jobs)==26
    assert all(load(root/'cases'/j['case_id']/'result.json')['status']=='PASS' for j in main_jobs)
    fine=[]
    for job in manifest['jobs']:
        name=job['case_id'];pop=job['population'];folder=root/'cases'/name;result=folder/'result.json'
        item=dict(case_id=name,population=pop,run_status='NOT_RUN',convergence_status='NOT_RUN',time_ns=None,steps=None,wall_s=None,
            primary_relative=None,event_relative=None,tails_nm=None,tail_threshold_status='NOT_RUN',minimum_gap_nm=None,raw_velocity_max_m_s=None)
        if result.exists():
            d=load(result);item.update(run_status=d['status'],time_ns=d['time_ns'],steps=d['steps'],wall_s=d['wall_s'],
                minimum_gap_nm=d.get('minimum_gap_nm'),raw_velocity_max_m_s=d.get('diagnostics',{}).get('maximum_raw_velocity_m_s'))
            if d['status']=='PASS':
                reference=root/'cases'/('mesh2p5_timehalf_'+pop);old=load(reference/'result.json')
                fine_raw=raw_audit(folder);coarse_raw=raw_audit(reference)
                c=list(csv.DictReader((reference/'history.csv').open()))[-1];f=list(csv.DictReader((folder/'history.csv').open()))[-1]
                errors={};events={};tails={}
                for host in pop:
                    for suffix in ('line_nm','tail_nm','swept_area_nm2','area_tensor_norm_nm2'):
                        key=host+'_'+suffix;x=float(c[key]);y=float(f[key]);errors[key]=abs(x-y)/max(abs(y),1e-20)
                    t0=old['first_contact_ns'][host];t1=d['first_contact_ns'][host]
                    events[host]=abs(t1-t0)/t1 if t0 is not None and t1 is not None else None
                    tails[host]=dict(coarse_nm=float(c[host+'_tail_nm']),fine_nm=float(f[host+'_tail_nm']),relative=errors[host+'_tail_nm'])
                    assert abs(tails[host]['coarse_nm']-coarse_raw['native_tail_nm'][host])<1e-10
                    assert abs(tails[host]['fine_nm']-fine_raw['native_tail_nm'][host])<1e-10
                with np.load(reference/'initial.npz') as x,np.load(folder/'initial.npz') as y:
                    for key in ('positions_nm','tags','ids','burgers','planes','cell_nm'):np.testing.assert_array_equal(x[key],y[key])
                item.update(primary_relative=errors,event_relative=events,tails_nm=tails,initial_graph_identity='EXACT',
                    actual_array_audit=fine_raw,
                    tail_threshold_status='PASS' if all(t['relative']<=limit for t in tails.values()) else 'FAIL',
                    convergence_status='PASS' if max(errors.values())<=limit and all(v is not None and v<=event_limit for v in events.values()) else 'UNRESOLVED')
        fine.append(item)
    log=(root/'spatial_extension.log').read_text() if (root/'spatial_extension.log').exists() else ''
    receipt=load(root/'spatial_extension_receipt.json') if (root/'spatial_extension_receipt.json').exists() else None
    producer=load(root/'spatial_extension_result.json')['status'] if (root/'spatial_extension_result.json').exists() else 'CENSORED' if 'CENSORED:' in log and receipt is not None else 'UNRESOLVED'
    missing=[x['case_id'] for x in fine if x['run_status']=='NOT_RUN']
    completed=sum(x['run_status']=='PASS' for x in fine)
    strict=load(root/'g0_regression/strict_full_field_replay.json');g0=load(root/'g0_regression/result.json')
    force=load(root/'force_audit/result.json');origin=load(root/'origin_sensitivity_result.json');release=load(base/'reverse_release.json');barrier=load(base/'ctb_ablation.json')
    contrasts=load(base/'paired_contrasts.json')['mesh2p5_timehalf'];causal=load(base/'causal_resolution.json')
    audit=load(base/'case_audits.json');regional=load(base/'regional_plastic.json')
    fine_force_path=root/'force_audit_mesh1p25/result.json';fine_force=load(fine_force_path)['status'] if fine_force_path.exists() else 'NOT_RUN'
    # No inferred qualification at a missing MT state or at an untested fine dt.
    overall='FAIL' if any(x['run_status']=='FAIL' for x in fine) or force['status']=='FAIL' or g0['status']=='FAIL' else 'UNRESOLVED'
    if all(x['convergence_status']=='PASS' for x in fine) and fine_force=='PASS' and all(gate[k]=='PASS' for k in ('source_inheritance_contact_release_topology','native_force_causality','G0_regression','numerical_origin_sensitivity')):overall='PASS'
    verdict=dict(gate='W4-G1',status=overall,review_date='2026-10-30',frozen_primary_relative=limit,frozen_event_relative=event_limit,
        main_native_cases=26,completed_1p25_native_cases=completed,planned_1p25_native_cases=3,missing_1p25_cases=missing,
        tail_threshold='FAIL' if any(x['tail_threshold_status']=='FAIL' for x in fine) else 'PASS' if all(x['tail_threshold_status']=='PASS' for x in fine) else 'NOT_RUN',
        spatial_extension_resource_status=producer,source_identity_contact_release_topology=gate['source_inheritance_contact_release_topology'],
        native_force_at_2p5=force['status'],native_force_at_1p25_MT=fine_force,temporal_refinement_at_1p25='NOT_RUN',
        G0_original_scientific_contract=g0['status'],G0_auxiliary_1e10_full_field_replay=strict['status'],
        periodic_G2_admission='BLOCKED' if overall!='PASS' else 'PASS',W4_G2='NOT_RUN',W6_collective='NOT_RUN',VCoNi='NOT_QUALIFIED',
        qualification_scope='One shared planar I0 CTB, nonperiodic direct force, one affine ordinary propagation loop per side. General reactions, exact ledger restart and material strengthening are not qualified.',
        evidence_policy='All old FAIL/CENSORED/UNRESOLVED and raw trajectories retained. Missing measurements are null, not zero. No new dynamics launched by this reviewer.')
    def save(name,value):(out/name).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
    save('gate_review.json',verdict);save('spatial_convergence.json',fine)
    save('resource_review.json',dict(status=producer,original_manifest=manifest,actual_receipt=receipt,unstarted_cases=missing,
        reason='Controller requires the full frozen 1500 s allowance before starting each case; it does not start a pending case when the 3000 s envelope has less than that remaining.' if producer=='CENSORED' else None))
    with (out/'tail_convergence.csv').open('w',newline='',encoding='utf-8') as f:
        fields=['case_id','host','run_status','coarse_2p5_tail_nm','fine_1p25_tail_nm','relative_error','threshold','criterion_status'];w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for item in fine:
            for host in item['population']:
                t=None if item['tails_nm'] is None else item['tails_nm'][host]
                w.writerow(dict(case_id=item['case_id'],host=host,run_status=item['run_status'],coarse_2p5_tail_nm=None if t is None else t['coarse_nm'],fine_1p25_tail_nm=None if t is None else t['fine_nm'],relative_error=None if t is None else t['relative'],threshold=limit,criterion_status='NOT_RUN' if t is None else 'PASS' if t['relative']<=limit else 'FAIL'))
    max_balance=max(max(c['cross_force_pair_residual'],c['direct_pair_residual']) for c in force['checks'])
    max_pk=max(c['PK_tensor_relative_error'] for c in force['checks'])
    max_closure=max(c['plastic_closure'] for c in audit);region_error=max(c['closure'] for c in regional)
    lines=['# W4-G1 最终中文综合审查报告','',f'审查日期：2026-10-30。最终裁决：**{overall}**。本报告取代进行中的状态摘要，但不覆盖任何历史记录。',
        '', '## 1. 执行结果与最终裁决','',
        f'共享CTB的材料侧身份继承、可释放接触及原生阶段一致性已经实现。修复版003完成26例真实ExaDiS主矩阵，另有3例原点敏感性反控；1.25 nm细化实际完成{completed}/3例。没有重跑已经完成的动力学。本轮只审查已有输出、生成报告与归档。',
        '', f'原2%尾段门槛未改变。已完成细化中存在超限，尾段验收为{verdict["tail_threshold"]}；因此即使每例正常结束、节点运动且守恒通过，G1整体仍为{overall}。W4-G2准入BLOCKED，W6未启动。',
        '', '## 2. 1.25 nm 尾段复核','',
        '相对差采用 |L(1.25)-L(2.5)|/|L(1.25)|，与原审查一致。不是改用分母更大值、有限带宽密度或另选时刻取得通过。共同原生终点为0.025 ns，初始多边形、Burgers、滑移面和节点标签逐数组相同。',
        '', '| 人口/侧别 | 原生计算 | 2.5 nm尾段 (nm) | 1.25 nm尾段 (nm) | 相对差 (%) | 原2%裁决 |','|---|---|---:|---:|---:|---|']
    for item in fine:
        for host in item['population']:
            t=None if item['tails_nm'] is None else item['tails_nm'][host]
            lines.append(f'| {item["population"]}/{host} | {item["run_status"]} | '+('未测 | 未测 | 未测 | NOT_RUN |' if t is None else f'{t["coarse_nm"]:.8f} | {t["fine_nm"]:.8f} | {100*t["relative"]:.5f} | '+('PASS' if t['relative']<=limit else 'FAIL')+' |'))
    lines+=['', f'细化控制器状态：{producer}。未启动案例：{", ".join(missing) if missing else "无"}。原预算为1500 s/例、总3000 s，并要求启动前预留完整单例预算；没有临时增额，也没有拼接轨迹。原生案例NOT_RUN与控制器CENSORED分开记账。',
        '', '5→2.5 nm 的M+T/T尾段6.34%超限以及首版UNRESOLVED审查保留。即使某个后续侧别通过，也不能删除已经失败的侧别或用平均数掩盖它。',
        '', '## 3. 已完成的真实模型与代码','',
        'M/T满足共同[111]轴60°的FCC Σ3取向关系，CTB为同一个z=2000 nm平面。两侧采用同一宏观sigma_zz=1.1 GPa、b=0.255 nm、mu=54.6 GPa、nu=0.324和同一GLIDE迁移率；参数属于有效FCC基准，未标定为VCoNi。',
        '', '每侧是普通闭合传播环，不是Frank–Read源。M-only、T-only与M+T对应初始子网络完全相同。M+T在一个ExaDiS网络中共同求解直接弹性力；按各自真实旋转的滑移系受力，而非拼接两条独立通道。非周期计算采用无限介质直接力，没有有限自由表面或各向异性界面像力。',
        '', 'material_identity.py按真实线段Burgers、滑移法向及仿射源平面唯一恢复来源；tag只是可复用地址。host_layer、origin_region、current_region、source/lineage和canonical interface字段保留。原生分裂、粗化、跨域碰撞拒绝及回滚图均有实际证据。未知反应显式拒绝并停止，不静默猜测材料侧。',
        '', 'shared_contact.py对同一个CTB按所属侧施加相反单侧法向约束。反力通过原生GLIDE响应求得，保留合法界面切向运动；驱动力反向后解除接触。不固定所有接触节点，也不在穿透后裁回坐标。积分、Collision、TopologySerial和Remesh接受阶段均检查身份、半空间和Burgers守恒。',
        '', '7项身份/能力回归、7项原生接触/分裂/粗化/拒绝回归、4项力原点修复准入有真实日志。整数邻接向量化与610个真实构型、133747个节点观测逐项相同；它不改变浮点几何谓词或物理运动律。此资格限于当前可唯一辨识的仿射单滑移源；复杂多源反应与精确累计账本恢复未获资格。',
        '', '## 4. M-only、T-only、M+T传播与弹性消融','',
        '下面为已完成2.5 nm匹配参考，尚不能作为最终收敛后的材料预测。比较按来源分开，不用M+T总位错数量翻倍制造“强化”或活动增强结论。',
        '', '| 侧别 | M+T减单侧扫掠面积 (nm²) | 相对单侧 (%) | 首次接触时刻差 (ns) |','|---|---:|---:|---:|']
    for host in ('M','T'):
        c=contrasts[host];lines.append(f'| {host} | {c["area_difference_nm2"]:.8f} | {100*c["relative_to_single"]:.6f} | {c["first_contact_difference_ns"]:.9g} |')
    lines+=['', '本夹具中跨侧相互作用略微促进传播，而不是必然使运动更少。全部已保存M+T构型执行原生共同应力PK差分、全网络减单侧力、直接跨侧段对以及瞬时移除跨侧力后的原生速度检查。它证明这里存在真实跨侧弹性反馈，但不证明Σ3独有、四厚度趋势、屈服强化或HDI。',
        '', f'2.5 nm全构型力互反最大残差{max_balance:.6g}，原10^-8门槛下通过；PK相对误差最大{max_pk:.6g}。1.25 nm M+T真实轨迹力审查为{fine_force}，不能借用静态重网格预检代替。',
        '', '## 5. 界面阻碍和真实反向释放','',
        '| 人口/侧别 | 有CTB面积 (nm²) | 无CTB面积 (nm²) | 面积降低 (%) |','|---|---:|---:|---:|']
    for b in barrier:lines.append(f'| {b["population"]}/{b["host"]} | {b["wall_area_nm2"]:.5f} | {b["no_ctb_area_nm2"]:.5f} | {100*b["area_reduction"]:.5f} |')
    lines+=['', '无CTB是移除约束而保留原始滑移流形的参考，不是真实CTB透射。反向加载均从同一真实前向末态开始；节点、线段、Burgers和滑移面逐数组一致，没有重建源或插入子环。',
        '', '| 人口 | 各侧首次全部脱离时间 (原生ns) | 实际末态绑定 | 裁决 |','|---|---|---|---|']
    for r in release:lines.append(f'| {r["population"]} | {r["first_full_release_ns"]} | {r["restart_actual_forward_final"]} | {r["status"]} |')
    lines+=['', f'主矩阵原生全局塑性闭合最大绝对残差{max_closure:.6g}，M/T来源贡献之和恢复原生全局塑性的最大残差{region_error:.6g}。接触反力与CTB承载反力、三阶段速度、真实接受dt和阶段面积账本保存。接触反力不等于实验定义的HDI应力。',
        '', '## 6. 为什么仍不能通过G1','',
        '数值正常结束、无非法穿透、Burgers守恒和全局塑性闭合，是必要条件；它们不保证界面接触尾段已经空间收敛。本例扫掠面积和事件时刻较稳定，精确位于CTB上的线段长度仍随网格变化超过2%。不将尾段改成更容易通过的带宽密度，也不选择另一终点。',
        '', 'M-only原始速度峰值从5 nm约2.10 km/s、2.5 nm约3.30 km/s增加到1.25 nm约8.62 km/s；1.25 nm独立时间减半尚NOT_RUN。这个变化提示需要核验显式Euler时间误差、接触端点离散及Remesh局部效应，但当前证据不足以确定根因，不能将它解释成已经证实的物理局域化。',
        '', '阶段净面积变化接近零不能证明局部尾段不受重网格影响。下一次修复需要同时记录接触端点和尾段在每个原生阶段前后的变化。当前数据不能支持把界面储存差异用于定量强化或VCoNi参数拟合。',
        '', '## 7. 失败保留与G0复核','',
        '002旧版变形末态力互反约1.68e-6，直接原生段对复算也失败；平移计算原点诊断将其降至约1e-11。003仅统一平移传给原生力内核的坐标及cell origin，实际图、CTB、材料、载荷和迁移率不变，不配平或投影力。旧26例FAIL资格、180秒CENSORED、原日志和构型全部保留。',
        '', f'G0原生产源码保持不变；新A0、固定节点/原生FR、接触/三速度和四条HI/HII无壁/双壁轨迹按原科学门禁为{g0["status"]}。额外10^-10完整力/速度逐字段复现为{strict["status"]}，原失败日志保留，没有放宽该项。两种标准分别裁决。G0原80例R5不覆盖。',
        '', '早期报告对缺席相的空线段绘图发生异常，原traceback保留；仅修复渲染后生成新目录。没有更改数据或将绘图成功当成科学PASS。4000 m/s限速影响G0时序的已有发现仍有效；未标定的无上限线性GLIDE时钟不是VCoNi真实时间预测。',
        '', '## 8. 最终门禁与后续方向','',
        f'- 实现/原生阶段准入：PASS；主矩阵26例数值完成：PASS；精确接触尾段2%验收：{verdict["tail_threshold"]}；G1整体：{overall}。',
        '- 1.25 nm未测数据保持NOT_RUN/null；资源控制器CENSORED独立保存。不启动新的细化，不重复已完成计算。',
        '- 下一恢复任务应先在固定1.25 nm网格下完成独立时间步敏感性，并审计接触端点的分阶段Remesh行为；保持2%门槛和物理参数，失败则继续保留。此处仅提出方案，没有执行新的动力学。',
        '- 在尾段与新增时间资格通过前，G2周期层片、W6集体矩阵均不得由本结果放行。复杂CTB反应、VCoNi标定、Nye/GND、强化、局域化和HDI均NOT_RUN或NOT_QUALIFIED。整个W4–W9未完成。',
        '', '## 9. 可复现证据与归档','',
        f'本次实际VM根：`{root}`。原始配置、求解命令、退出码、源码快照、初末构型、关键帧、接受步字段、事件与分析脚本均在对应案例目录。',
        '', '002完整失败归档与003完整修复版归档分别保存为G1_native_evidence_20261030_002.tar和G1_native_evidence_20261030_003.tar。小型ZIP/native视图省略密集accepted_velocity块，完整TAR保留。库存与ZIP CRC核对，不使用大型SHA256重复扫描。实际归档尺寸和文件数见项目EVIDENCE_INDEX。',
        '', '下列图仅展示2.5 nm匹配参考的真实构型与全过程，明确属于尚未通过最终尾段资格的数据；不是1.25 nm缺失MT的替代，也不是论文最终资格图。']
    (out/'W4_G1_FINAL_REPORT_ZH.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    def inline(s):
        s=html.escape(s);s=re.sub(r'\*\*(.+?)\*\*',r'<strong>\1</strong>',s);return re.sub(r'`(.+?)`',r'<code>\1</code>',s)
    body='';table=False
    for line in lines:
        if line.startswith('|'):
            if not line.replace('|','').replace('-','').replace(':','').strip():continue
            cells=[c.strip() for c in line.strip('|').split('|')]
            if not table:body+='<table><thead><tr>'+''.join('<th>'+inline(c)+'</th>' for c in cells)+'</tr></thead><tbody>';table=True
            else:body+='<tr>'+''.join('<td>'+inline(c)+'</td>' for c in cells)+'</tr>'
            continue
        if table:body+='</tbody></table>';table=False
        if line.startswith('# '):body+='<h1>'+inline(line[2:])+'</h1>'
        elif line.startswith('## '):body+='<h2>'+inline(line[3:])+'</h2>'
        elif line:body+='<p>'+inline(line)+'</p>'
    if table:body+='</tbody></table>'
    for name in ('G1_native_networks','G1_actual_glide_planes','G1_propagation_and_paired_difference','G1_reverse_release'):
        body+='<figure><img src="../reference_figures/'+name+'.png"><figcaption>'+name+' — 2.5 nm真实参考，最终资格UNRESOLVED</figcaption></figure>'
    css='body{max-width:1120px;margin:35px auto;padding:20px;font:16px/1.75 "Microsoft YaHei",sans-serif;color:#243441}h1,h2{color:#165477}h2{border-bottom:1px solid #d4dee4;margin-top:2em}table{border-collapse:collapse;width:100%;font-size:14px}th,td{border:1px solid #d5dfe4;padding:8px}th{background:#eef3f6}strong{color:#9a5d16}code{overflow-wrap:anywhere;background:#eef1f3}img{width:100%}figure{margin:32px 0}p{white-space:pre-wrap}'
    (out/'W4_G1_FINAL_REPORT_ZH.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>W4-G1 最终中文综合审查</title><style>'+css+'</style><body>'+body+'</body></html>',encoding='utf-8')
    print('FINAL_G1_REVIEW',overall,'fine_completed',completed,'resource',producer,'tail',verdict['tail_threshold'])


if __name__=='__main__':main()
