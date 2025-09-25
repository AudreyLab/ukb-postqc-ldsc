# README.md

## UK Biobank GWAS Post-QC & LDSC Pre-munge — Python Utilities

This repository provides two Python scripts that reimplement the functionality of the original C/C++ utilities (`post_filter.C` and `prep_munge.C`) for **post-GWAS QC** and **preparation of summary statistics for LDSC**.  
They are designed for outputs from **REGENIE** (or similar GWAS tools) and mimic the logic used in the PLINK2 statistical routines (HWE tests, EMAC filters).

---

## 1. Installation

Create a virtual environment and install the minimal dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -U pandas
```

Both scripts are standalone; no heavy dependencies are required.  
`pandas` is optional but recommended for downstream processing.

---

## 2. Prepare REGENIE Input

Concatenate the 22 chromosome-level REGENIE outputs into a single file, keeping only the header from the first file. Example:

```bash
(
  zcat gwasA.chr1_PHENO.regenie.gz | head -1
  for chr in $(seq 1 22); do
    zcat gwasA.chr${chr}_PHENO.regenie.gz | tail -n +2
  done
) > gwasA_PHENO.regenie.tsv
```

This produces a single table with all variants for phenotype `PHENO`.

---

## 3. Post-filtering (EMAC and HWE)

`post_filter.py` removes variants that fail **Effective Minor Allele Count (EMAC)** and **Hardy–Weinberg equilibrium (HWE)** thresholds.

```bash
python post_filter.py   --in gwasA_PHENO.regenie.tsv   --out gwasA_PHENO.filtered.tsv.gz   --emac-min 100   --hwe-minp 1e-12   --id-col ID --aaf-col A1FREQ --n-col N --info-col INFO
```

### Notes
- **EMAC**: if `INFO` contains `EMAC`, that value is used; otherwise EMAC ≈ `2 * N * min(AAF, 1-AAF)`.  
- **HWE**: if `INFO` contains `HWE`, that value is used; otherwise an exact mid-p test (Wigginton et al., 2005) is computed from genotype counts if available (`N_HOMREF`, `N_HET`, `N_HOMALT`).  
- Variants missing EMAC/HWE values are kept (no filtering).

---

## 4. Pre-munge for LDSC

`prep_munge.py` reformats the filtered GWAS results into the columns expected by LDSC (`SNP, A1, A2, N, Z, P, BETA, SE`).

```bash
python prep_munge.py   --in gwasA_PHENO.filtered.tsv.gz   --out gwasA_PHENO.pre_munged.tsv.gz   --id-col ID --a1-col A1 --a2-col A2 --n-col N --p-col P   --effect-col Effect --info-col INFO   --logistic   # include this flag if the GWAS was logistic
```

### Logic for BETA/SE extraction
- If `INFO` has `BETA` and `SE`: use directly.  
- Else if `INFO` has `LOGOR` and `SE`: use as log-odds ratio.  
- Else if `INFO` has `OR` and `CI95L/CI95U`: compute `BETA = ln(OR)` and `SE = (ln(UCL)-ln(LCL))/3.92`.  
- Else if linear model: use `Effect` column as BETA, with `SE` from `INFO`.  

Z-scores are computed as `Z = BETA / SE` where possible.

---

## 5. Workflow Example

```bash
# 1. Concatenate chromosomes
(
  zcat gwasA.chr1_PHENO.regenie.gz | head -1
  for chr in $(seq 1 22); do
    zcat gwasA.chr${chr}_PHENO.regenie.gz | tail -n +2
  done
) > gwasA_PHENO.regenie.tsv

# 2. Post-QC filter
python post_filter.py   --in gwasA_PHENO.regenie.tsv   --out gwasA_PHENO.filtered.tsv.gz   --emac-min 100 --hwe-minp 1e-12

# 3. Pre-munge for LDSC
python prep_munge.py   --in gwasA_PHENO.filtered.tsv.gz   --out gwasA_PHENO.pre_munged.tsv.gz   --logistic
```

Now the file `gwasA_PHENO.pre_munged.tsv.gz` is ready for LDSC’s `munge_sumstats.py`.

---

## 6. Column Mapping

By default, the scripts expect REGENIE-like column names:

- `ID` (variant ID)  
- `A1`, `A2` (alleles)  
- `N` (sample size)  
- `P` (p-value)  
- `Effect` (effect size, OR or BETA)  
- `INFO` (key–value annotations: EMAC, HWE, etc.)  
- `A1FREQ` (effect allele frequency)

If your files use different column names, adjust with the `--*-col` options.

---

## 7. License / Provenance

- These scripts are a Python reimplementation of your C/C++ utilities (`post_filter.C`, `prep_munge.C`) and mimic PLINK2 routines (HWE test, chi-square).  
- Original PLINK2 headers (`plink2_base`, `plink2_string`, `plink2_stats`) are under permissive licenses (Boost, BSD, LGPL). This Python code is a clean rewrite without reuse of C/C++ source.
