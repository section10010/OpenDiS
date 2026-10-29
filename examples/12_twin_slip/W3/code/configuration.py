"""One validated source of units, geometry, loading and native numerical settings."""
from copy import deepcopy
import json
from pathlib import Path
import numpy as np


DEFAULT={
 'schema':'OpenDiS.W3.v1','name':'fr_reference','gate':'A1','kind':'FR_source',
 'material_scope':'effective FCC; not calibrated VCoNi',
 'b_m':2.55e-10,'mu_pa':54.6e9,'nu':.324,'core_a_b':6.,
 'Medge':64103.,'Mscrew':64103.,'friction_pa':0.,'vmax_m_s':4000.,
 'box_nm':[4000.,4000.,4000.],'periodic':[False,False,False],
 'stress_pa':[[0.,0.,0.],[0.,0.,0.],[0.,0.,600e6]],
 'phase':'matrix','mode':'HI','slip_index':None,'source_count':1,
 'source_gauge_nm':240.,'backbone_offset_nm':80.,'bow_nm':2.,
 'loop_radius_nm':60.,'threader_separation_nm':150.,'standoff_nm':10.,
 'emission_sign':1,'seed':42,'center_nm':None,
 'hT_nm':200.,'hM_nm':400.,'sources_twin':60,'sources_matrix':120,
 'placement_guard_nm':5.,'walls':[],
 'force_mode':'direct','grid':[48,48,48],
 'integrator':'Trapezoid','maxseg_nm':20.,'minseg_nm':4.,'rann_nm':.6,
 'rtol_nm':.1,'maxdt_ps':.5,'firstdt_ps':.05,'max_node_step_nm':.5,
 'contact_tol_nm':1e-6,'end_ns':1.,'max_steps':100000,'max_wall_s':900.,
 'threads':2,'frame_interval_ns':.02,'stop_emissions':0,
 'topology':True,'collision':True,'remesh':True,
}


def config(**updates):
    unknown=set(updates)-set(DEFAULT)
    if unknown:raise ValueError('Unknown config fields: '+str(sorted(unknown)))
    c=deepcopy(DEFAULT);c.update(deepcopy(updates));return validate(c)


def validate(c):
    if set(c)!=set(DEFAULT):raise ValueError('Configuration must have exactly the W3 schema fields')
    if c['schema']!=DEFAULT['schema']:raise ValueError('Unsupported schema')
    for k in ('b_m','mu_pa','core_a_b','Medge','Mscrew','source_gauge_nm','backbone_offset_nm',
              'loop_radius_nm','maxseg_nm','minseg_nm','rann_nm','rtol_nm','maxdt_ps','firstdt_ps',
              'max_node_step_nm','contact_tol_nm','end_ns','max_wall_s','frame_interval_ns'):
        if not np.isfinite(c[k]) or c[k]<=0:raise ValueError('Positive finite field required: '+k)
    if not 0<c['nu']<.5:raise ValueError('Invalid Poisson ratio')
    if c['friction_pa']!=0:raise ValueError('Nonzero bulk friction is not yet implemented; never silently ignored')
    if c['vmax_m_s'] is not None and c['vmax_m_s']<=0:raise ValueError('vmax must be positive or null')
    box=np.array(c['box_nm'],float);s=np.array(c['stress_pa'],float)
    if box.shape!=(3,) or not np.isfinite(box).all() or np.any(box<=0):raise ValueError('Invalid box')
    if s.shape!=(3,3) or not np.isfinite(s).all() or not np.allclose(s,s.T,rtol=0,atol=1e-6):raise ValueError('Symmetric finite stress tensor required')
    if len(c['periodic'])!=3 or any(type(v) is not bool for v in c['periodic']):raise ValueError('Three periodic booleans required')
    if c['minseg_nm']>=c['maxseg_nm']:raise ValueError('minseg must be below maxseg')
    if c['kind'] not in ('FR_source','propagation_loop','prepared_threader','multilamellar'):raise ValueError('Unknown kind')
    if c['phase'] not in ('matrix','twin') or c['mode'] not in ('HI','HII','soft'):raise ValueError('Invalid phase/mode')
    if c['slip_index'] is not None and not 0<=c['slip_index']<12:raise ValueError('Invalid slip index')
    if c['emission_sign'] not in (-1,1):raise ValueError('Emission sign must be +/-1')
    for k in ('source_count','sources_twin','sources_matrix','seed','threads','max_steps','stop_emissions'):
        if type(c[k]) is not int or c[k]<0:raise ValueError('Nonnegative integer field required: '+k)
    if c['threads']<1 or c['max_steps']<1:raise ValueError('Positive threads/steps required')
    if c['force_mode'] not in ('direct','fft'):raise ValueError('Unknown force model')
    if c['force_mode']=='fft' and not all(c['periodic']):raise ValueError('FFT benchmark requires periodic cell')
    if c['integrator'] not in ('Trapezoid','EulerForward'):raise ValueError('Unsupported integrator')
    if c['walls'] and c['integrator']!='EulerForward':raise ValueError('I0 requires the event-limited native Euler path')
    if c['walls'] and c['periodic'][2]:raise ValueError('Minimal I0 uses physical independent walls inside nonperiodic z')
    ids=[]
    for wall in c['walls']:
        if set(wall)!={'interface_id','z_nm','inward','enabled'}:raise ValueError('Invalid interface fields')
        if wall['inward'] not in (-1,1) or not 0<wall['z_nm']<box[2]:raise ValueError('Invalid wall')
        ids.append(wall['interface_id'])
    if len(ids)!=len(set(ids)):raise ValueError('Unique interface_id required')
    if c['kind']=='multilamellar':
        period=c['hT_nm']+c['hM_nm'];q=box[2]/period
        if min(c['hT_nm'],c['hM_nm'])<=0 or not np.isclose(q,round(q),rtol=0,atol=1e-9):raise ValueError('Box must contain full lamella periods')
    return c


def load(path):
    return config(**json.loads(Path(path).read_text(encoding='utf-8-sig')))


def save_json(path,value):
    def convert(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False,default=convert),encoding='utf-8')
