#!/usr/bin/env python3
import sys, gzip, argparse, math

def open_any(path: str):
    if path == "-" or path is None:
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

def parse_info(s: str):
    out = {}
    for kv in s.split(";"):
        if not kv:
            continue
        if "=" in kv:
            k, v = kv.split("=", 1)
            out[k.strip()] = v.strip()
    return out

def f(x, default=None, cast=float):
    try:
        return cast(x)
    except:
        return default

def beta_se_from_info(info: dict, effect_col_val: str, is_logistic: bool):
    """
    Stratégies:
    - Si INFO a BETA et SE -> direct
    - Sinon si INFO a LOGOR et SE -> BETA=LOGOR
    - Sinon si INFO a OR et CI95L/CI95U -> BETA=ln(OR), SE=(ln(CI95U)-ln(CI95L))/3.92
    - Sinon si linear et effet déjà beta -> BETA=effect_col_val, SE dans INFO
    """
    # 1) direct
    b = f(info.get("BETA"))
    se = f(info.get("SE"))
    if b is not None and se is not None:
        return b, se

    # 2) LOGOR + SE
    logor = f(info.get("LOGOR"))
    se2   = f(info.get("SE"))
    if logor is not None and se2 is not None:
        return logor, se2

    # 3) OR + CI
    orv = f(info.get("OR"))
    lcl = f(info.get("CI95L"))
    ucl = f(info.get("CI95U"))
    if orv and lcl and ucl and lcl>0 and ucl>0:
        b3 = math.log(orv)
        se3 = (math.log(ucl) - math.log(lcl)) / 3.92
        return b3, se3

    # 4) linear: l'effet est déjà beta, SE dans INFO
    if not is_logistic:
        b4 = f(effect_col_val)
        se4 = f(info.get("SE"))
        if b4 is not None and se4 is not None:
            return b4, se4

    return None, None

def main():
    ap = argparse.ArgumentParser(description="Prepare sumstats for LDSC from REGENIE-like TSV")
    ap.add_argument("--in", dest="inp", default="-", help="Input TSV/TSV.GZ with header")
    ap.add_argument("--out", dest="out", default="-", help="Output TSV/TSV.GZ or '-'")
    ap.add_argument("--id-col", default="ID", help="SNP ID column (default: ID)")
    ap.add_argument("--a1-col", default="A1", help="Effect allele column (default: A1)")
    ap.add_argument("--a2-col", default="A2", help="Other allele column (default: A2)")
    ap.add_argument("--n-col", default="N", help="Sample size (effective) column (default: N)")
    ap.add_argument("--p-col", default="P", help="P-value column name (default: P)")
    ap.add_argument("--effect-col", default="Effect", help="Effect column (OR or BETA depending on model)")
    ap.add_argument("--info-col", default="INFO", help="INFO key=val;... column (default: INFO)")
    ap.add_argument("--logistic", action="store_true", help="Set if the model was logistic (Effect usually OR)")
    args = ap.parse_args()

    inp = open_any(args.inp)
    out = sys.stdout if args.out == "-" else (gzip.open(args.out, "wt") if args.out.endswith(".gz") else open(args.out, "w"))

    header = inp.readline().rstrip("\n")
    cols = header.split("\t")
    idx = {c:i for i,c in enumerate(cols)}

    def gi(name): 
        return idx.get(name)

    id_i   = gi(args.id_col)
    a1_i   = gi(args.a1_col)
    a2_i   = gi(args.a2_col)
    n_i    = gi(args.n_col)
    p_i    = gi(args.p_col)
    eff_i  = gi(args.effect_col)
    info_i = gi(args.info_col)

    print("SNP\tA1\tA2\tN\tZ\tP\tBETA\tSE", file=out)

    for line in inp:
        line = line.rstrip("\n")
        if not line: 
            continue
        parts = line.split("\t")
        snp  = parts[id_i]   if id_i   is not None else None
        a1   = parts[a1_i]   if a1_i   is not None else None
        a2   = parts[a2_i]   if a2_i   is not None else None
        nval = f(parts[n_i]) if n_i    is not None else None
        pval = f(parts[p_i]) if p_i    is not None else None

        info = parse_info(parts[info_i]) if (info_i is not None and info_i < len(parts)) else {}
        effect_val = parts[eff_i] if eff_i is not None else None
        beta, se = beta_se_from_info(info, effect_val, args.logistic)

        z = None
        if beta is not None and se is not None and se > 0:
            z = beta / se

        # Output even if some fields are missing; LDSC munge_sumstats fera ses propres vérifs
        print(f"{snp}\t{a1}\t{a2}\t{nval if nval is not None else ''}\t{z if z is not None else ''}\t{pval if pval is not None else ''}\t{beta if beta is not None else ''}\t{se if se is not None else ''}", file=out)

    if out is not sys.stdout:
        out.close()

if __name__ == "__main__":
    main()

