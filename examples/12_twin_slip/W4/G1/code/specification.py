"""Explicit capability refusal for the actually admitted G1 I0 model."""
import numpy as np


def validate_spec(spec,c):
    allowed={'schema','scope','contact_enabled','initial_geometry_updates','common_native_updates','interface','regions','criteria'}
    if not allowed<=set(spec) or set(spec)-allowed-{'force_evaluation_origin_nm'}:raise ValueError('Unknown/missing G1 specification field; no silent physics fallback')
    if spec['schema']!='W4.G1.SharedEvolution.v1':raise ValueError('Unsupported G1 schema')
    if type(spec['contact_enabled']) is not bool:raise ValueError('contact_enabled must be boolean')
    if c['kind']!='propagation_loop' or c['source_count']!=1:raise NotImplementedError('G1 admission supports one ordinary propagation loop per region, not FR or arbitrary source populations')
    if c['integrator']!='EulerForward' or any(c['periodic']) or c['walls'] or c['force_mode']!='direct':
        raise NotImplementedError('G1 admission supports nonperiodic direct-force Euler with one shared CTB only')
    if 'force_evaluation_origin_nm' in spec:
        origin=np.array(spec['force_evaluation_origin_nm'],float)
        if origin.shape!=(3,) or not np.isfinite(origin).all():raise ValueError('Invalid native force evaluation origin')
    if set(spec['initial_geometry_updates'])!={'maxseg_nm'} or spec['initial_geometry_updates']['maxseg_nm']<=0:
        raise ValueError('Explicit common initial polygon mesh required')
    interface=spec['interface']
    if set(interface)!={'id','z_nm'} or interface['id']!='CTB0' or not 0<interface['z_nm']<c['box_nm'][2]:
        raise NotImplementedError('Only one named planar CTB0 normal to sample z is admitted')
    if [(r['region'],r['phase'],r['allowed_normal']) for r in spec['regions']]!=[('M','matrix',[0,0,-1]),('T','twin',[0,0,1])]:
        raise ValueError('Incorrect M/T phase, side or normal contract')
    for r in spec['regions']:
        if set(r)!={'region','phase','center_nm','allowed_normal'} or np.array(r['center_nm']).shape!=(3,):raise ValueError('Invalid region specification')
    crit=spec['criteria']
    if not np.isclose(c['contact_tol_nm'],crit['penetration_nm'],rtol=0,atol=0):raise ValueError('Contact tolerance and frozen penetration criterion disagree')
    for key in ('source_plane_tol_nm','primary_refinement_relative','event_relative','cross_force_balance_relative','plastic_absolute','campaign_max_wall_s'):
        if not np.isfinite(crit[key]) or crit[key]<=0:raise ValueError('Invalid scientific criterion '+key)
