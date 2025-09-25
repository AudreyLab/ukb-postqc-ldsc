#!/usr/bin/env python3
import sys, gzip, io, argparse, re, math
from typing import Dict, Tuple, Optional

def open_any(path: str):
    if path == "-" or path is None:
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def parse_info(info_str: str) -> Dict[str, str]:
    # key1=val1;key2=val2;...
    out = {}
    for kv in info_str.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out

# -------- HWE exact test (Wigginton et al., 2005), comme PLINK2 (thread-safe) --------
def hwe_pvalue(obs_hom1: int, obs_het: int, obs_hom2: int) -> float:
    """Exact test HWE (biallélique). Retourne P(mid-p)."""
    # swap pour que 'homc' = le génotype le moins fréquent
    obs_homc = min(obs_hom1, obs_hom2)
    obs_homo = max(obs_hom1, obs_hom2)
    obs_heto = obs_het
    n = obs_homc + obs_homo + obs_heto
    if n == 0:
        return 1.0
    rare = 2 * obs_homc + obs_heto
    # prob init au mode
    mid = int(rare * (2 * n - rare) / (2 * n))
    if (rare % 2) ^ (mid % 2):
        mid += 1
    # proba au mode
    prob_mid = 1.0
    prob = prob_mid
    # sum over distribution
    p_total = prob_mid
    # left tail
    i = mid
    while i > 0:
        prob *= i * (i - 1) / (4.0 * ( ( (rare - i) // 2 ) + 1 ) * ( ( ( (2*n - rare) - i) // 2 ) + 1 ))
        p_total += prob
        i -= 2
    # right tail
    prob = prob_mid
    i = mid
    while i <= rare - 2:
        prob *= ( ( (rare - i) // 2 ) * ( ( (2*n - rare) - i) // 2 ) * 4.0 ) / ( (i + 2) * (i + 1) )
        p_total += prob
        i += 2
    # mid-p
    tail = 0.0
    # recompute tail including all probabilities <= prob_mid
    prob = prob_mid
    i = mid
    tail += prob
    # left
    probL = prob_mid
    j = mid
    while j > 0:
        probL *= j * (j - 1) / (4.0 * ( ( (rare - j) // 2 ) + 1 ) * ( ( ( (2*n - rare) - j) // 2 ) + 1 ))
        if probL <= prob_mid + 1e-15:
            tail += probL
        j -= 2
    # right
    probR = prob_mid
    k = mid
    while k <= rare - 2:
        probR *= ( ( (rare - k) // 2 ) * ( ( (2*n - rare) - k) // 2 ) * 4.0 ) / ( (k + 2) * (k + 1) )
        if probR <= prob_mid + 1e-15:
            tail += probR
        k += 2
    return min(1.0, tail)

def infer_counts_from_info(info: Dict[str,str]) -> Optional[Tuple[int,int,int]]:
    """
    Essaie de déduire (hom_ref, het, hom_alt) depuis INFO si présent.
    Clés supportées (ajustables): N_HOMREF, N_HET, N_HOMALT  ou  OBS_HOM1, OBS_HET, OBS_HOM2.
    """
    keys_sets = [
        ("N_HOMREF","N_HET","N_HOMALT"),
        ("OBS_HOM1","OBS_HET","OBS_HOM2"),
        ("hom_ref","het","hom_alt"),
    ]
    for a,b,c in keys_sets:
        if a in info and b in info and c in info:
            try:
                return (int(float(info[a])), int(float(info[b])), int(float(info[c])))
            except:
                pass
    return None

def emac_from_fields(aaf: Optional[float], n_eff: Optional[float]) -> Optional[float]:
    """
    EMAC (approx.) = 2 * N_eff * min(MAF, 1-MAF)
    Si seule AAF (ALT allele freq) dispo, MAF = min(AAF,1-AAF).
    """
    if aaf is None or n_eff is None:
        return None
    maf = aaf if aaf <= 0.5 else 1.0 - aaf
    return 2.0 * n_eff * maf

def parse_float_safe(x: str) -> Optional[float]:
    try:
        return float(x)
    except:
        return None

def main():
    ap = argparse.ArgumentParser(description="Post-filter REGENIE by EMAC and HWE p-value (Python)")
    ap.add_argument("--in", dest="inp", default="-", help="Input TSV/TSV.GZ (with header), or '-' for stdin")
    ap.add_argument("--out", dest="out", default="-", help="Output TSV/TSV.GZ or '-' for stdout")
    ap.add_argument("--emac-min", type=float, default=100.0, help="Minimum EMAC (default 100)")
    ap.add_argument("--hwe-minp", type=float, default=1e-12, help="Minimum HWE p-value (default 1e-12)")
    ap.add_argument("--id-col", default="ID", help="SNP id column (default: ID)")
    ap.add_argument("--a1-col", default="A1", help="Effect allele col name if needed (optional)")
    ap.add_argument("--a2-col", default="A2", help="Other allele col name (optional)")
    ap.add_argument("--aaf-col", default="A1FREQ", help="ALT/effect allele frequency column (default: A1FREQ)")
    ap.add_argument("--n-col", default="N", help="Effective N column (or total N) (default: N)")
    ap.add_argument("--info-col", default="INFO", help="INFO key=val;... column (default: INFO)")
    args = ap.parse_args()

    inp = open_any(args.inp)
    out_handle = sys.stdout if args.out == "-" else (gzip.open(args.out, "wt") if args.out.endswith(".gz") else open(args.out, "w"))

    header = inp.readline().rstrip("\n")
    if not header:
        return
    cols = header.split("\t")
    col_idx = {c:i for i,c in enumerate(cols)}
    print(header, file=out_handle)

    id_i   = col_idx.get(args.id_col, None)
    aaf_i  = col_idx.get(args.aaf_col, None)
    n_i    = col_idx.get(args.n_col, None)
    info_i = col_idx.get(args.info_col, None)

    kept = 0
    for line in inp:
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        # parse INFO
        info = {}
        if info_i is not None and info_i < len(parts):
            info = parse_info(parts[info_i])

        # HWE : soit directement via INFO['HWE'], soit via comptes géno si présents
        hwe_p = None
        if "HWE" in info:
            hwe_p = parse_float_safe(info["HWE"])
        else:
            counts = infer_counts_from_info(info)
            if counts:
                homref, het, homalt = counts
                hwe_p = hwe_pvalue(homref, het, homalt)

        # EMAC : si INFO['EMAC'] existe, on l’utilise, sinon EMAC ≈ 2*N_eff*MAF
        emac = None
        if "EMAC" in info:
            emac = parse_float_safe(info["EMAC"])
        else:
            aaf = None
            neff = None
            if aaf_i is not None and aaf_i < len(parts):
                aaf = parse_float_safe(parts[aaf_i])
            if n_i is not None and n_i < len(parts):
                neff = parse_float_safe(parts[n_i])
            emac = emac_from_fields(aaf, neff)

        pass_hwe  = (hwe_p is None) or (hwe_p >= args.hwe_minp)   # si inconnu: ne pas filtrer
        pass_emac = (emac is None) or (emac >= args.emac_min)     # si inconnu: ne pas filtrer

        if pass_hwe and pass_emac:
            print(line, file=out_handle)
            kept += 1

    if out_handle is not sys.stdout:
        out_handle.close()

if __name__ == "__main__":
    main()

