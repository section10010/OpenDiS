"""FCC geometry in physical nm; conversion to b units occurs only at native boundary."""
import numpy as np

AXES=np.array([[1,-1,0],[1,1,-2],[1,1,1]],float)
Q=AXES/np.linalg.norm(AXES,axis=1)[:,None]
ROT=np.array([[.5,-np.sqrt(3)/2,0],[np.sqrt(3)/2,.5,0],[0,0,1.]])
PLANES=np.array([[1,1,1],[1,-1,1],[1,1,-1],[-1,1,1]],float)
DIRS=np.array([[1,1,0],[1,-1,0],[1,0,1],[1,0,-1],[0,1,1],[0,1,-1]],float)


def systems(phase,stress):
    R=ROT@Q if phase=='twin' else Q;sigma=np.array(stress,float);out=[]
    for n0 in PLANES:
        for b0 in DIRS:
            if abs(n0@b0)>1e-12:continue
            n=R@(n0/np.linalg.norm(n0));b=R@(b0/np.linalg.norm(b0))
            mode='soft' if abs(n[2])>.999 else 'HII' if abs(b[2])<1e-9 else 'HI'
            out.append(dict(index=len(out),n=n,b=b,mode=mode,rss_pa=float(b@sigma@n)))
    return out


def pure_shear(phase,mode,tau_pa,index=None):
    ss=systems(phase,np.zeros((3,3)))
    s=ss[index] if index is not None else next(s for s in ss if s['mode']==mode)
    return tau_pa*(np.outer(s['b'],s['n'])+np.outer(s['n'],s['b']))


def choose(c,phase=None,index=None):
    ss=systems(phase or c['phase'],c['stress_pa']);index=c['slip_index'] if index is None else index
    return ss[index] if index is not None else max([s for s in ss if s['mode']==c['mode']],key=lambda s:abs(s['rss_pa']))


def source_local(c,s,sign=1,return_sign=1):
    n=s['n'];b=s['b']*(-np.sign(s['rss_pa']) if abs(s['rss_pa'])>1e-8 else 1.)
    q=b.copy() if s['mode']=='soft' else np.cross([0.,0.,1.],n)
    q/=np.linalg.norm(q);q*=sign;p=np.cross(n,q)
    L=c['source_gauge_nm'];offset=c['backbone_offset_nm'];d=c['maxseg_nm']
    u=np.linspace(0,1,max(3,int(np.ceil(L/d))+1))
    arm=(-L/2+L*u)[:,None]*q-c['bow_nm']*np.sin(np.pi*u)[:,None]*p
    points=list(arm);mobile_edges=len(arm)-1
    corners=[L/2*q,L/2*q+return_sign*offset*n,-L/2*q+return_sign*offset*n,-L/2*q]
    for x,y in zip(corners[:-1],corners[1:]):
        for a in np.linspace(0,1,int(np.ceil(np.linalg.norm(y-x)/d))+1)[1:]:points.append(x+a*(y-x))
    points=np.array(points[:-1]);constraints=np.full(len(points),7,int);constraints[1:len(arm)-1]=0
    planes=[];roles=[]
    for j,x in enumerate(points):
        tangent=points[(j+1)%len(points)]-x;normal=np.cross(b,tangent)
        normal=n if np.linalg.norm(normal)<1e-9 else normal/np.linalg.norm(normal)
        planes.append(n if j<mobile_edges else normal)
        roles.append('source_arm' if j<mobile_edges else 'fixed_continuation')
    return points,constraints,np.array(planes),b,q,p,roles


def apportion(total,slots,rng):
    values=np.full(slots,total//slots,int);values[rng.permutation(slots)[:total%slots]]+=1;return values


def build(c):
    rng=np.random.default_rng(c['seed']);nodes=[];links=[];sources=[];roles=[];parent=[]
    def add(phase,index,center=None,slab=None,sign=1,return_sign=1):
        s=choose(c,phase,index);kind=c['kind'];n=s['n']
        if kind in ('FR_source','multilamellar'):
            pts,con,pn,b,q,p,role=source_local(c,s,sign,return_sign)
        else:
            b=s['b']*(-np.sign(s['rss_pa']) if abs(s['rss_pa'])>1e-8 else 1.)
            q=np.cross([0.,0.,1.],n) if s['mode']!='soft' else b.copy();q/=np.linalg.norm(q);q*=sign;p=np.cross(n,q)
            if kind=='propagation_loop':
                r=c['loop_radius_nm'];num=max(12,int(np.ceil(2*np.pi*r/c['maxseg_nm'])))
                theta=np.arange(num)*2*np.pi/num
                pts=r*(np.cos(theta)[:,None]*q+np.sin(theta)[:,None]*p)
            else:
                enabled=c['walls']
                if len(enabled)<2:raise ValueError('Threader geometry requires two wall locations (may be disabled after building)')
                h=max(w['z_nm'] for w in enabled)-min(w['z_nm'] for w in enabled)
                extent=(h/2-c['standoff_nm'])/abs(p[2]);half=c['threader_separation_nm']/2
                if extent<=0:raise ValueError('Threader standoff too large')
                corners=[-half*q-extent*p,half*q-extent*p,half*q+extent*p,-half*q+extent*p]
                pts=[]
                for x,y in zip(corners,corners[1:]+corners[:1]):
                    for t in np.linspace(0,1,int(np.ceil(np.linalg.norm(y-x)/c['maxseg_nm']))+1)[:-1]:pts.append(x+t*(y-x))
                pts=np.array(pts)
            con=np.zeros(len(pts),int);pn=np.tile(n,(len(pts),1));role=['mobile_bulk']*len(pts)
        if center is None:
            center=np.array(c['box_nm'])/2 if c['center_nm'] is None else np.array(c['center_nm'],float)
        else:center=np.array(center,float)
        if slab is not None:
            lo,hi=slab;guard=c['placement_guard_nm'];a=lo+guard-pts[:,2].min();z=hi-guard-pts[:,2].max()
            if a>z:raise ValueError('Source envelope does not fit; no silent resizing')
            center[2]=a+center[2]*(z-a)
        pts+=center
        if np.any(pts<0) or np.any(pts>=np.array(c['box_nm'])):raise ValueError('Initial network outside cell')
        start=len(nodes);sid=len(sources)
        for j,x in enumerate(pts):
            nodes.append([*x,int(con[j])]);links.append([start+j,start+(j+1)%len(pts),*b,*pn[j]])
            roles.append(role[j]);parent.append(sid)
        source=dict(source_id=sid,kind=kind,phase=phase,mode=s['mode'],slip_index=s['index'],
                    center_nm=center.tolist(),b=b.tolist(),n=n.tolist(),q=q.tolist(),p=p.tolist(),
                    rss_pa=float(b@np.array(c['stress_pa'])@n),initial_node_indices=list(range(start,start+len(pts))),
                    anchor_indices=(start+np.flatnonzero(con==7)).tolist(),
                    mobile_arm_end_indices=(start+np.flatnonzero(con==7)[:2]).tolist() if kind in ('FR_source','multilamellar') else [],
                    return_sign=return_sign)
        sources.append(source)
    if c['kind']=='multilamellar':
        n=int(round(c['box_nm'][2]/(c['hT_nm']+c['hM_nm'])))
        for phase,total,offset,h in [('twin',c['sources_twin'],0,c['hT_nm']),('matrix',c['sources_matrix'],c['hT_nm'],c['hM_nm'])]:
            per_layer=apportion(total,n,rng);pool=np.tile(np.arange(12),total//12+1)[:total];rng.shuffle(pool);cursor=0
            for layer,num in enumerate(per_layer):
                band=(np.arange(num)+rng.uniform(.15,.85,num))/max(num,1)
                lo=layer*(c['hT_nm']+c['hM_nm'])+offset
                for k,z in enumerate(band):
                    add(phase,int(pool[cursor]),[rng.uniform(.05,.95)*c['box_nm'][0],rng.uniform(.05,.95)*c['box_nm'][1],z],
                        (lo,lo+h),1 if cursor%2==0 else -1,1 if cursor%2==0 else -1);cursor+=1
    else:
        for k in range(c['source_count']):
            center=None
            if c['source_count']>1:
                center=np.array(c['box_nm'])/2;center[0]+=(k-(c['source_count']-1)/2)*(c['source_gauge_nm']+2*c['backbone_offset_nm'])
            add(c['phase'],c['slip_index'],center,sign=c['emission_sign'])
    rn=np.array(nodes,float).reshape(-1,4);lk=np.array(links,float).reshape(-1,8)
    return rn,lk,sources,np.array(roles),np.array(parent,int)


def residual(n,ids,b):
    r=np.zeros((n,3));np.add.at(r,ids[:,0],b);np.add.at(r,ids[:,1],-b);return np.linalg.norm(r,axis=1)


def displacements(pos,ids,box,periodic):
    dr=pos[ids[:,1]]-pos[ids[:,0]];box=np.array(box);periodic=np.array(periodic)
    dr-=np.round(dr/box)*box*periodic;return dr


def components(n,ids):
    parent=np.arange(n)
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in ids:
        a,b=find(i),find(j)
        if a!=b:parent[b]=a
    groups={}
    for i in range(n):groups.setdefault(find(i),[]).append(i)
    return [np.array(v,int) for v in groups.values()]
