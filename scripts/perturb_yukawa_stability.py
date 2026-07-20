import sys; sys.path.insert(0,".")

import numpy as np

from scripts.anarchic_bauer_s1 import DEFAULT_K_GEV, SCENARIOS, _draw_bauer_matrix, _fn_c_values
from scripts.instrument_epsK_phase import _instrument_draw, draw_nelson_barr_yukawas
from scripts.run_rs_anarchy import _load_pdg_targets
from warpConfig.wavefuncs import f_IR

targets=_load_pdg_targets(); sc=SCENARIOS["S1"]; ymax=sc["y_max"]
MKK=3000.0; xi=1.0; eps=MKK/xi/DEFAULT_K_GEV
mf,cf,jf=3.0,3.0,10.0
def evalpt(Yu,Yd,fQ,fu,fd):
    r=_instrument_draw(Yu,Yd,fQ,fu,fd,MKK,xi,targets,mf,cf,jf)
    return r["ratio_eps_K"], r["passes_pdg"], r["ratio_dm_K"]
def fit(Yu,Yd,rng):
    c_u3=float(rng.uniform(-sc["c_max"],0.5))
    cQ,cu,cd=_fn_c_values(c_u3,eps,targets,Y_u=Yu,Y_d=Yd,common_cd=False,rng=rng,c_jitter=0.0)
    return cQ,cu,cd,f_IR(cQ,eps),f_IR(cu,eps),f_IR(cd,eps)

def base_NB(rng):
    for _ in range(2000):
        Yu,Yd=draw_nelson_barr_yukawas(rng,y_max=ymax,rho_cp=1.0,eta_leak=0.0)
        cQ,cu,cd,fQ,fu,fd=fit(Yu,Yd,rng); re,ok,_=evalpt(Yu,Yd,fQ,fu,fd)
        if ok and re<=1.0: return Yu,Yd,fQ,fu,fd
def base_anarchic_phase(rng):
    best=None
    for _ in range(60000):
        Yu=_draw_bauer_matrix(rng,0.1,ymax); Yd=_draw_bauer_matrix(rng,0.1,ymax)
        cQ,cu,cd,fQ,fu,fd=fit(Yu,Yd,rng); re,ok,dm=evalpt(Yu,Yd,fQ,fu,fd)
        if ok and re<=1.0 and dm>0.05:  # survivor AND phase-aligned (large dm_K = large |C4|)
            return Yu,Yd,fQ,fu,fd
    return None

def deltacrit(Yu,Yd,fQ,fu,fd,rng,kind,mask="offdiag12",ndir=40):
    # perturb Yd: kind in {real,imag}; mask entries; scan delta; return delta at 50% pass
    idx=[(0,1),(1,0)] if mask=="offdiag12" else [(i,j) for i in range(3) for j in range(3) if i!=j]
    deltas=np.geomspace(1e-4,1.0,22); frac=[]
    for d in deltas:
        npass=0
        for _ in range(ndir):
            Z=np.zeros((3,3),complex)
            for (i,j) in idx:
                z=rng.normal()+ (1j*rng.normal() if kind=="imag" else 0.0)
                if kind=="imag": z=1j*rng.normal()
                Z[i,j]=z
            Ypp=Yd+d*Z
            re,ok,_=evalpt(Yu,Ypp,fQ,fu,fd)
            if re<=1.0: npass+=1
        frac.append(npass/ndir)
    frac=np.array(frac)
    below=np.where(frac<0.5)[0]
    return float(deltas[below[0]]) if len(below) else float("inf"), deltas, frac

rng=np.random.default_rng(7)
print("=== finding base survivors ===",flush=True)
nb=base_NB(rng); print("NB survivor found",flush=True)
an=base_anarchic_phase(rng); print("anarchic phase-aligned survivor:", "found" if an else "NOT FOUND",flush=True)

for name,base in [("NB (real down)",nb),("anarchic tuned",an)]:
    if base is None: continue
    Yu,Yd,fQ,fu,fd=base
    dcr_i,_,_=deltacrit(Yu,Yd,fQ,fu,fd,rng,"imag")
    dcr_r,_,_=deltacrit(Yu,Yd,fQ,fu,fd,rng,"real")
    print(f"\n{name}: delta_crit(IMAG/CP pert on 1-2)={dcr_i:.2e}  delta_crit(REAL pert)={dcr_r:.2e}",flush=True)
    print(f"   fractional tuning ~ 1/delta_crit(imag) = {1/dcr_i:.0f}",flush=True)
    print(f"   anisotropy real/imag = {dcr_r/dcr_i:.1f} (>>1 => protected manifold; ~1 => isolated tuned)",flush=True)
print("DONE",flush=True)
