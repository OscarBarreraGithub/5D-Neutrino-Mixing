#!/usr/bin/env python
"""Tuning-radius / stability-under-perturbation measure for RS eps_K survivors
(Opus spec, RS-FLAVOR-ALIGNMENT-2026-07). Perturb Y_d' = Y_d + delta*s_Y*(M o Z),
Z unit-Frobenius, ensembles {full, re, im}; measure delta_crit (50% pass), the CP
anisotropy A=dcrit(re)/dcrit(im), the budget exponent p (dcrit ~ eps_max^p:
0=symmetry, 1=accident), and d_sym (fraction of real dirs with eps_K identically 0).
Reuses _instrument_draw for the eps_K eval. Two budgets via ratio thresholds."""
import sys; sys.path.insert(0,".")
import numpy as np, math
from warpConfig.wavefuncs import f_IR
from scripts.anarchic_bauer_s1 import _draw_bauer_matrix, _fn_c_values, SCENARIOS, DEFAULT_K_GEV
from scripts.instrument_epsK_phase import _instrument_draw, draw_nelson_barr_yukawas, _rms
from scripts.run_rs_anarchy import _load_pdg_targets

targets=_load_pdg_targets(); sc=SCENARIOS["S1"]; ymax=sc["y_max"]
MKK=3000.0; xi=1.0; epsG=MKK/xi/DEFAULT_K_GEV; mf,cf,jf=3.0,3.0,10.0
THR_RESID=1.0; THR_OBS=2.228e-3/6.70e-5  # ratio_eps_K thresholds for the two budgets
def ratioeps(Yu,Yd,fQ,fu,fd):
    r=_instrument_draw(Yu,Yd,fQ,fu,fd,MKK,xi,targets,mf,cf,jf); return r["ratio_eps_K"],r["passes_pdg"]
def fit(Yu,Yd,rng):
    c_u3=float(rng.uniform(-sc["c_max"],0.5))
    cQ,cu,cd=_fn_c_values(c_u3,epsG,targets,Y_u=Yu,Y_d=Yd,common_cd=False,rng=rng,c_jitter=0.0)
    return fQ_fu_fd(cQ,cu,cd)
def fQ_fu_fd(cQ,cu,cd): return f_IR(cQ,epsG),f_IR(cu,epsG),f_IR(cd,epsG)
def base_NB(rng):
    for _ in range(3000):
        Yu,Yd=draw_nelson_barr_yukawas(rng,y_max=ymax,rho_cp=1.0,eta_leak=0.0)
        f=fit(Yu,Yd,rng); re,ok=ratioeps(Yu,Yd,*f)
        if ok and re<=1.0: return Yu,Yd,f
def base_anar(rng):
    for _ in range(80000):
        Yu=_draw_bauer_matrix(rng,0.1,ymax); Yd=_draw_bauer_matrix(rng,0.1,ymax)
        f=fit(Yu,Yd,rng); re,ok=ratioeps(Yu,Yd,*f)
        # phase-aligned survivor: pass AND large |C4| (need R_K sizeable). use dm proxy via 2nd eval field
        if ok and re<=1.0:
            r=_instrument_draw(Yu,Yd,*f,MKK,xi,targets,mf,cf,jf)
            if r["ratio_dm_K"]>0.05: return Yu,Yd,f
def gen_Z(rng,kind):
    if kind=="re": Z=rng.normal(size=(3,3))+0j
    elif kind=="im": Z=1j*rng.normal(size=(3,3))
    else: Z=(rng.normal(size=(3,3))+1j*rng.normal(size=(3,3)))/np.sqrt(2)
    return Z/ (np.linalg.norm(Z)+1e-30)
def dcrit(Yu,Yd,f,rng,kind,thr,ndir=96):
    sY=_rms(Yd); deltas=np.geomspace(1e-5,1.0,26); 
    for d in deltas:
        npass=0
        for _ in range(ndir):
            Ypp=Yd+d*sY*gen_Z(rng,kind); re,ok=ratioeps(Yu,Ypp,*f)
            if re<=thr and ok: npass+=1
        if npass/ndir<0.5: return d
    return float("inf")
def d_sym_frac(Yu,Yd,f,rng,ndir=64):
    # fraction of REAL directions along which eps_K stays at machine-zero at delta=0.1
    sY=_rms(Yd); c=0
    for _ in range(ndir):
        Ypp=Yd+0.1*sY*gen_Z(rng,"re"); re,_=ratioeps(Yu,Ypp,*f)
        if re<1e-8: c+=1
    return c/ndir

rng=np.random.default_rng(11)
print("finding base points...",flush=True)
pts={"NB (real Y_d)":base_NB(rng),"anarchic tuned":base_anar(rng)}
for name,base in pts.items():
    if base is None: print(f"{name}: none"); continue
    Yu,Yd,f=base
    di=dcrit(Yu,Yd,f,rng,"im",THR_RESID); dr=dcrit(Yu,Yd,f,rng,"re",THR_RESID); dfu=dcrit(Yu,Yd,f,rng,"full",THR_RESID)
    di_o=dcrit(Yu,Yd,f,rng,"im",THR_OBS)
    p = math.log(di_o/di)/math.log(THR_OBS/THR_RESID) if (np.isfinite(di) and np.isfinite(di_o) and di>0) else float("nan")
    dsym=d_sym_frac(Yu,Yd,f,rng)
    A = dr/di if di>0 else float("inf")
    print(f"\n=== {name} ===",flush=True)
    print(f"  delta_crit  full={dfu:.2e}  real={dr:.2e}  imag={di:.2e}",flush=True)
    print(f"  Delta_tune (s_Y/dcrit_full) = {_rms(Yd)/dfu:.0f}" if np.isfinite(dfu) else "  Delta_tune = inf (protected)",flush=True)
    print(f"  CP anisotropy A=dcrit(re)/dcrit(im) = {A:.1f}   (~1 tuned, >>1 protected)",flush=True)
    print(f"  budget exponent p (dcrit~eps_max^p) = {p:.2f}   (0 symmetry, 1 accident, 0.5 tangent)",flush=True)
    print(f"  d_sym: fraction of REAL dirs with eps_K==0 (machine) = {dsym:.2f}   (~1 => 9-dim real flat, ~0 => isolated)",flush=True)
print("DONE",flush=True)
