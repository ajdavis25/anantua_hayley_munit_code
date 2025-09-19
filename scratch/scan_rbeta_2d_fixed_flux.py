#!/usr/bin/env python3
# scan_rbeta_2d_fixed_flux.py  (beta_crit pinned to 1.0)
import pathlib as p, re, subprocess, tempfile, h5py, numpy as np

IP = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/ipole"
POS0_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos0.par"
POS1_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos1.par"

Rhigh_grid = [15,20,25]
sigma_grid = [0.5,1.0,1.5,2.0]
targets = {"pos0":0.530, "pos1":0.520}
goals   = {"pos0":(3.93e-3,-5.60e-3), "pos1":(1.62e-3,-2.81e-3)}

def read_param(txt,k,default=None):
    import re
    m=re.search(r"(?m)^\s*"+re.escape(k)+r"\s+(\S+)",txt)
    return m.group(1) if m else default
def set_or_replace(txt,k,v):
    import re
    line=f"{k} {v}"
    if re.search(r"(?m)^\s*"+re.escape(k)+r"\s+\S+",txt):
        return re.sub(r"(?m)^\s*"+re.escape(k)+r"\s+\S+.*$", line, txt)
    return txt.rstrip()+"\n"+line+"\n"

def run(par):
    p1 = subprocess.Popen([IP,"-par",str(par)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    out,_=p1.communicate()
    return p1.returncode,out

def Ftot(h5):
    with h5py.File(h5,"r") as H: return float(H["/Ftot"][()])

def metrics(h5):
    with h5py.File(h5,"r") as H:
        s=float(H["/header/scale"][()]); F=float(H["/Ftot"][()])
        P=np.array(H["/pol"]).transpose(1,0,2)*s
        I,Q,U,V=[P[:,:,i].sum() for i in range(4)]
        return F, (Q*Q+U*U)**0.5/F, V/F

def retune_once(txt, target):
    tmp1=p.Path(tempfile.mkstemp(suffix=".par")[1]); tmp1.write_text(txt)
    rc1,log1=run(tmp1)
    if rc1!=0: return None, None, log1, tmp1
    out=p.Path(read_param(txt,"outfile")); 
    if not out.exists(): return None, None, "outfile missing", tmp1
    F=Ftot(str(out))
    mu=float(read_param(txt,"M_unit")); mu2=mu*(target/(F if F!=0 else target))
    tuned=set_or_replace(txt,"M_unit",f"{mu2:.6E}")
    tmp2=p.Path(tempfile.mkstemp(suffix=".par")[1]); tmp2.write_text(tuned)
    rc2,log2=run(tmp2)
    if rc2!=0 or not out.exists(): return None, None, log1+"\n"+log2, tmp2
    return metrics(str(out)), tuned, None, tmp2

def scan(base_par_path, tag):
    base=p.Path(base_par_path).read_text()
    assert read_param(base,"electronModel")=="RBETA"
    tgtF=targets[tag]; P_t,V_t=goals[tag]
    base=set_or_replace(base,"beta_crit","1.0")  # pin here
    base=set_or_replace(base,"quench_output","0")
    base=set_or_replace(base,"positronRatio","0" if tag=="pos0" else "1")
    print(f"\n=== {tag} @ fixed flux {tgtF:.3f} ===")
    print(" Rhigh  sigma   |    M_unit    |   Ftot   |    P/I    |    V/I    |  |Δ|")
    best=None
    for R in Rhigh_grid:
        for S in sigma_grid:
            txt=set_or_replace(base,"trat_small","1")
            txt=set_or_replace(txt,"trat_large",str(R))
            txt=set_or_replace(txt,"sigma_cut",str(S))
            res,tuned,err,_par=retune_once(txt,tgtF)
            if err:
                print(f"  R={R:>3}  s={S:>3} : ERROR")
                # show last ~25 lines of log to pinpoint
                for line in (err.splitlines()[-25:]):
                    print("    "+line)
                continue
            (F,Pov,Vov), tuned_txt = res, tuned
            mu=float(read_param(tuned_txt,"M_unit"))
            d=abs(Pov-P_t)+abs(Vov-V_t)
            print(f"  {R:>5}  {S:>5}  | {mu:.3E} |  {F:6.3f}  | {Pov:8.3e} | {Vov:8.3e} | {d:8.3e}")
            if best is None or d<best[0]: best=(d,R,S,mu,F,Pov,Vov)
    if best:
        _,R,S,mu,F,P,V=best
        print(f"\n>>> best {tag}: Rhigh={R}  sigma_cut={S}  ->  F={F:.3f}  P/I={P:.3e}  V/I={V:.3e}  (M_unit={mu:.3E})")

if __name__=="__main__":
    scan(POS0_PAR,"pos0")
    scan(POS1_PAR,"pos1")
