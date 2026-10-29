"""Uniform scientific plots from retained real ExaDiS arrays and histories."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from analyze import read,arrays,history,export_csv
from geometry import choose
from events import topology

COL={'none':'#536779','lower':'#2477ac','upper':'#c18732','both':'#ad3d42'}


def draw2(ax,folder,file='final.npz'):
    c=read(folder/'config.json');ss=choose(c);n=ss['n'];q=np.cross([0.,0.,1.],n);q/=np.linalg.norm(q);p=np.cross(n,q)
    a=arrays(folder/file);x=a['positions_nm']-np.array(c['box_nm'])/2;xy=np.stack([x@q,x@p],axis=1)/1000
    ij=a['ids'];colors=np.full(len(ij),'#246791',object)
    for w in c['walls']:
        hit=np.abs(a['positions_nm'][:,2]-w['z_nm'])<=c['contact_tol_nm']
        if w['enabled']:colors[np.all(hit[ij],axis=1)]='#b7373c'
        ax.axhline((w['z_nm']-c['box_nm'][2]/2)/p[2]/1000,color='#b18532' if w['enabled'] else '#b2bac2',lw=.8,ls='--',zorder=0)
    ini=arrays(folder/'initial.npz');xx=ini['positions_nm']-np.array(c['box_nm'])/2;yy=np.stack([xx@q,xx@p],axis=1)/1000
    ax.add_collection(LineCollection(yy[ini['ids']],colors='#b8bec5',linewidths=.65,zorder=1))
    ax.add_collection(LineCollection(xy[ij],colors=colors,linewidths=1.,zorder=2))
    ax.set(xlim=(-.62,.62),ylim=(-.62,.62),aspect='equal');ax.grid(alpha=.15,lw=.5)


def density(a,edges):
    ij=a['ids'];x=a['positions_nm'];dr=x[ij[:,1]]-x[ij[:,0]];length=np.linalg.norm(dr,axis=1)
    result=np.zeros(len(edges)-1)
    for (i,j),l in zip(ij,length):
        z0=x[i,2];z1=x[j,2]
        if abs(z1-z0)<1e-12:
            k=np.searchsorted(edges,z0,side='right')-1
            if 0<=k<len(result):result[k]+=l
        else:
            lo=np.minimum((edges[:-1]-z0)/(z1-z0),(edges[1:]-z0)/(z1-z0))
            hi=np.maximum((edges[:-1]-z0)/(z1-z0),(edges[1:]-z0)/(z1-z0))
            result+=l*np.maximum(0,np.minimum(hi,1)-np.maximum(lo,0))
    assert abs(result.sum()-length.sum())<1e-7*max(1,length.sum()),'z-slice line conservation'
    return result


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--analysis',type=Path);a=p.parse_args();r=a.root;o=a.out;o.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,'savefig.facecolor':'white','pdf.fonttype':42})
    families=[(k,m) for k in ('loop','threader') for m in ('HI','HII')];masks=('none','lower','upper','both')
    def save(fig,name):
        fig.savefig(o/(name+'.png'),dpi=200,bbox_inches='tight');fig.savefig(o/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(4,4,figsize=(11.8,11.4),sharex=True,sharey=True)
    for row,(kind,mode) in enumerate(families):
        for col,mask in enumerate(masks):
            folder=r/f'B_{kind}_{mode}_{mask}_tau300_t100';ax=axes[row,col];draw2(ax,folder)
            if row==0:ax.set_title({'none':'No wall','lower':'Lower wall','upper':'Upper wall','both':'Both walls'}[mask],fontweight='bold')
            if col==0:ax.set_ylabel(f'{kind} / {mode}\np (um)')
            if row==3:ax.set_xlabel('q: CTB trace (um)')
    fig.suptitle('Existing dislocation propagation | identical initial network per row | 0.10 ns\nFCC, pure RSS 300 MPa; I0 contact; h = 200 nm. Grey: initial; red: contact tail.',y=.99,fontsize=11)
    fig.tight_layout(rect=(0,0,1,.95));save(fig,'mechanism_comparison')
    fig=plt.figure(figsize=(12,8.4))
    for row,kind in enumerate(('loop','threader')):
        for col,mask in enumerate(masks):
            folder=r/f'B_{kind}_HI_{mask}_tau300_t100';c=read(folder/'config.json');aa=arrays(folder/'final.npz');x=(aa['positions_nm']-np.array(c['box_nm'])/2)/1000;ij=aa['ids']
            ax=fig.add_subplot(2,4,row*4+col+1,projection='3d');colors=np.full(len(ij),'#246791',object)
            for w in c['walls']:
                if w['enabled']:
                    hit=np.abs(aa['positions_nm'][:,2]-w['z_nm'])<=c['contact_tol_nm'];colors[np.all(hit[ij],axis=1)]='#b7373c'
                    gx,gy=np.meshgrid([-.6,.6],[-.6,.6]);gz=np.full_like(gx,(w['z_nm']-2000)/1000)
                    ax.plot_surface(gx,gy,gz,color='#b3a16e',alpha=.13,shade=False,linewidth=0)
            ax.add_collection3d(Line3DCollection(x[ij],colors=colors,linewidths=1.))
            ax.set(xlim=(-.6,.6),ylim=(-.6,.6),zlim=(-.6,.6),title=f'{kind} HI / {mask}',xlabel='x (um)',ylabel='y (um)',zlabel='z (um)');ax.view_init(22,36);ax.set_box_aspect((1,1,1));ax.tick_params(labelsize=6,pad=0)
    fig.suptitle('Real native final networks | common spatial scale | 0.10 ns',fontsize=12)
    fig.subplots_adjust(left=.03,right=.97,bottom=.06,top=.91,wspace=.08,hspace=.38);save(fig,'networks_3d')
    fig,axes=plt.subplots(2,4,figsize=(12,5.8),sharex=True)
    for col,(kind,mode) in enumerate(families):
        for mask in masks:
            h=history(r/f'B_{kind}_{mode}_{mask}_tau300_t100');t=[x['time_ns'] for x in h]
            axes[0,col].plot(t,[x['swept_area_nm2']/1e6 for x in h],color=COL[mask],label=mask,lw=1.2)
            axes[1,col].plot(t,[(x['lower_tail_nm']+x['upper_tail_nm'])/1000 for x in h],color=COL[mask],lw=1.2)
        axes[0,col].set_title(kind+' / '+mode,fontweight='bold');axes[1,col].set_xlabel('Time (ns)')
        for row in range(2):axes[row,col].grid(alpha=.18)
    axes[0,0].set_ylabel('Swept area (um2)');axes[1,0].set_ylabel('Total contact-tail length (um)');axes[0,3].legend(frameon=False,fontsize=8);fig.tight_layout();save(fig,'causal_metrics')
    fig,axes=plt.subplots(1,4,figsize=(12,3.3),sharey=True);profiles=[];edges=np.arange(0,4000.1,10.)
    for col,(kind,mode) in enumerate(families):
        for mask in ('none','both'):
            folder=r/f'B_{kind}_{mode}_{mask}_tau300_t100';aa=arrays(folder/'final.npz');d=density(aa,edges)/(4000*4000*10)*1e18
            axes[col].plot((edges[:-1]+5-2000),d/1e13,color=COL[mask],lw=1.,label=mask)
            profiles.extend(dict(case=folder.name,z_nm=float(z),rho_m2=float(v)) for z,v in zip(edges[:-1]+5,d))
        for z in (-100,100):axes[col].axvline(z,color='#9c8644',lw=.7,ls='--')
        axes[col].set(xlim=(-600,600),xlabel='z - layer centre (nm)',title=kind+' / '+mode);axes[col].grid(alpha=.15)
    axes[0].set_ylabel('Line density / 1e13 m-2');axes[-1].legend(frameon=False);fig.suptitle('Exact segment / z-bin intersections; 10 nm bins, common 4 x 4 um2 area',y=1.03);fig.tight_layout();save(fig,'local_storage');export_csv(o/'density_profiles.csv',profiles)
    source=r/'A1_ref240_80_tau200_wide';cfg=read(source/'config.json');ss=choose(cfg);q=np.cross([0.,0.,1.],ss['n']);q/=np.linalg.norm(q);pv=np.cross(ss['n'],q)
    lin=read(source/'source_lineage.json');births=[x for x in lin['events'] if x['event']=='native_detachment_candidate' and x['lineage_id'] in [c['lineage_id'] for c in lin['children'] if c['confirmed']]]
    files=[('initial.npz',0.)]+[(e['after_snapshot'],e['time_ns']) for e in births[:3]]
    allxy=[]
    for file,_ in files:
        xx=arrays(source/file)['positions_nm']-np.array(cfg['box_nm'])/2
        allxy.append(np.stack([xx@q,xx@pv],axis=1)/1000)
    allxy=np.concatenate(allxy);centre=(allxy.min(axis=0)+allxy.max(axis=0))/2;extent=np.ptp(allxy,axis=0).max()*.55
    fig,axes=plt.subplots(1,4,figsize=(12,4),sharex=True,sharey=True)
    for ax,(file,t) in zip(axes,files):
        aa=arrays(source/file);x=aa['positions_nm']-np.array(cfg['box_nm'])/2;xy=np.stack([x@q,x@pv],axis=1)/1000
        colors=np.full(len(aa['ids']),'#226990',object);con=aa['constraints'].ravel();colors[np.all(con[aa['ids']]==7,axis=1)]='#757d87'
        free=[g for g in topology(aa) if not g['anchored']];free.sort(key=lambda g:len(g['ix']),reverse=True)
        for g,color in zip(free,['#bf693c','#53924d','#8c6ba8','#a5a5a5']):colors[np.isin(aa['ids'][:,0],g['ix'])]=color
        ax.add_collection(LineCollection(xy[aa['ids']],colors=colors,linewidths=.9));ax.set(xlim=(centre[0]-extent,centre[0]+extent),ylim=(centre[1]-extent,centre[1]+extent),aspect='equal',title=f't = {t:.3f} ns',xlabel='q (um)');ax.grid(alpha=.15)
    axes[0].set_ylabel('p (um)');fig.suptitle('240 / 80 nm native FR source | 3 distinct large free offspring | 16 um domain\nBlue: mobile parent; grey: fixed closure; other colours: disconnected free components.',fontsize=10);fig.tight_layout(rect=(0,0,1,.88));save(fig,'source_evolution')
    fig,axes=plt.subplots(2,3,figsize=(9.5,6),sharex=True,sharey=True)
    for col,event in enumerate(births[:3]):
        for row,when in enumerate(('before','after')):
            aa=arrays(source/event[when+'_snapshot']);xx=aa['positions_nm']-np.array(cfg['box_nm'])/2;xy=np.stack([xx@q,xx@pv],axis=1)
            colors=np.full(len(aa['ids']),'#246791',object)
            for group in topology(aa):
                if not group['anchored']:colors[np.isin(aa['ids'][:,0],group['ix'])]='#b34341'
            fixed=aa['constraints'].ravel()==7;colors[np.all(fixed[aa['ids']],axis=1)]='#777e87'
            ax=axes[row,col];ax.add_collection(LineCollection(xy[aa['ids']],colors=colors,linewidths=1.1));ax.scatter(xy[fixed,0],xy[fixed,1],s=6,c='#777e87',marker='s')
            ax.set(xlim=(-180,180),ylim=(-180,180),aspect='equal',title=f'Emission {col+1}: {when} native {event["stage"]}');ax.grid(alpha=.15)
            if col==0:ax.set_ylabel('p (nm)')
            if row==1:ax.set_xlabel('q (nm)')
    fig.suptitle('Local neck / parent region only | full networks shown in source_evolution\nBlue: anchored mobile component; grey: fixed closure; red: disconnected free components.',fontsize=10);fig.tight_layout(rect=(0,0,1,.9));save(fig,'source_neck_zoom')
    data=read((a.analysis or r/'analysis')/'source_summary.json');fig,axes=plt.subplots(1,3,figsize=(11,3.4))
    scans=sorted([x for x in data if 'bracket' in x['case']],key=lambda x:x['tau_mpa'])
    for d in scans:
        t=d['first_detachment_ns'];axes[0].plot(d['tau_mpa'],2 if t is None else t,'v' if t is None else 'o',color='#ac3d42' if t is None else '#246791')
    axes[0].set(xlabel='Constant RSS (MPa)',ylabel='First detachment time (ns)',title='No emission at 50 / 75 MPa by 2 ns')
    reference=next(x for x in data if x['case']=='A1_ref240_80_tau200_wide')
    for ax,key,label in [(axes[1],'gauge_nm','Source gauge (nm)'),(axes[2],'return_nm','Fixed-return offset (nm)')]:
        dd=[reference]+[x for x in data if 'sensitivity' in x['case'] and (('source_gauge' in x['case'])==(key=='gauge_nm'))]
        dd.sort(key=lambda x:x[key]);ax.plot([x[key] for x in dd],[x['first_detachment_ns'] for x in dd],'o-',color='#246791');ax.set(xlabel=label,ylabel='First detachment time (ns)',title='200 MPa, all other inputs fixed')
    for ax in axes:ax.grid(alpha=.18)
    fig.tight_layout();save(fig,'source_activation_sensitivity');print('FIGURES_SAVED',len(list(o.glob('*.png'))))


if __name__=='__main__':main()
