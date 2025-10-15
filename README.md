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

### Overview

`post_QC.py` removes variants that fail **Effective Minor Allele Count (EMAC)** and **Hardy–Weinberg equilibrium (HWE)** thresholds.

### Basic Usage

```bash
python post_QC.py \
  --in gwasA_PHENO.regenie.tsv \
  --out gwasA_PHENO.filtered.tsv.gz \
  --emac-min 100 \
  --hwe-minp 1e-12
```

### All Available Arguments

#### 📁 Input/Output Files

| Argument | Default | Description |
|----------|---------|-------------|
| `--in` | `-` (stdin) | Input TSV or TSV.GZ file |
| `--out` | `-` (stdout) | Output TSV or TSV.GZ file |

#### 🎯 Filtering Thresholds

| Argument | Default | Description |
|----------|---------|-------------|
| `--emac-min` | `100.0` | Minimum EMAC threshold |
| `--hwe-minp` | `1e-12` | Minimum HWE p-value threshold |

#### 📊 Column Names

| Argument | Default | Description |
|----------|---------|-------------|
| `--id-col` | `Name` | SNP identifier column |
| `--aaf-col` | `AAF` | Allele frequency column |
| `--n-col` | `Num_Cases` | Sample size column |
| `--info-col` | `Info` | INFO field (key=val;...) |

#### ⚙️ Calculation Options

| Argument | Type | Description |
|----------|------|-------------|
| `--use-controls` | Flag | Use control genotypes for HWE test (recommended for case-control studies) |
| `--use-mac-from-info` | Flag | Use MAC from INFO field directly instead of calculating EMAC |

### Advanced Examples

#### Standard case-control GWAS filtering
```bash
python post_QC.py \
  --in gwasA_PHENO.regenie.tsv \
  --out gwasA_PHENO.filtered.tsv.gz \
  --emac-min 100 \
  --hwe-minp 1e-12 \
  --use-controls
```

#### Custom column names
```bash
python post_QC.py \
  --in gwasA_PHENO.regenie.tsv \
  --out gwasA_PHENO.filtered.tsv.gz \
  --id-col ID \
  --aaf-col A1FREQ \
  --n-col N \
  --info-col INFO
```

#### Using MAC from INFO field
```bash
python post_QC.py \
  --in gwasA_PHENO.regenie.tsv \
  --out gwasA_PHENO.filtered.tsv.gz \
  --use-mac-from-info
```

### Filtering Logic

#### EMAC Calculation
- **If `INFO` contains `EMAC`**: use that value directly
- **Else if `--use-mac-from-info` and `INFO` contains `MAC`**: use MAC as EMAC
- **Otherwise**: calculate as `EMAC = 2 × N × min(AAF, 1-AAF)`

#### HWE p-value Calculation
- **If `INFO` contains `HWE`**: use that value directly
- **Else if genotype counts available** (`Cases_Ref/Het/Alt` or `Controls_Ref/Het/Alt`): compute exact mid-p test (Wigginton et al., 2005)
- **Else if `INFO` contains genotype counts** (`N_HOMREF`, `N_HET`, `N_HOMALT`): compute from those

⚠️ **Important**: Variants with missing EMAC or HWE values are **REJECTED** (filtered out).

A variant is **KEPT** only if:
- `EMAC ≥ emac-min` **AND**
- `HWE p-value ≥ hwe-minp`

### Output Statistics

The script outputs filtering statistics to stderr:

```
Colonnes détectées: Name, Chr, Pos, Ref, Alt, ...
Filtres appliqués: EMAC >= 100.0, HWE p-value >= 1e-12

Variant 1 (rs367896724):
  EMAC=65757.751, pass=True
  HWE p-value=0.523, pass=True

============================================================
Statistiques de filtrage:
  Total variants: 199
  Filtrés (EMAC < 100.0): 45
  Filtrés (HWE p < 1e-12): 12
  Variants conservés: 142
  Taux de rétention: 71.36%
============================================================
```

---

## 4. Pre-munge for LDSC

`prep_munge.py` reformats the filtered GWAS results into the columns expected by LDSC (`SNP, A1, A2, N, Z, P, BETA, SE`).

```bash
python prep_munge.py \
  --in gwasA_PHENO.filtered.tsv.gz \
  --out gwasA_PHENO.pre_munged.tsv.gz \
  --id-col ID \
  --a1-col A1 \
  --a2-col A2 \
  --n-col N \
  --p-col P \
  --effect-col Effect \
  --info-col INFO \
  --logistic   # include this flag if the GWAS was logistic
```

### Logic for BETA/SE extraction
- If `INFO` has `BETA` and `SE`: use directly.
- Else if `INFO` has `LOGOR` and `SE`: use as log-odds ratio.
- Else if `INFO` has `OR` and `CI95L/CI95U`: compute `BETA = ln(OR)` and `SE = (ln(UCL)-ln(LCL))/3.92`.
- Else if linear model: use `Effect` column as BETA, with `SE` from `INFO`.

Z-scores are computed as `Z = BETA / SE` where possible.

---

## 5. Complete Workflow Example

```bash
# 1. Concatenate chromosomes
(
  zcat gwasA.chr1_PHENO.regenie.gz | head -1
  for chr in $(seq 1 22); do
    zcat gwasA.chr${chr}_PHENO.regenie.gz | tail -n +2
  done
) > gwasA_PHENO.regenie.tsv

# 2. Post-QC filter with recommended settings for case-control studies
python post_QC.py \
  --in gwasA_PHENO.regenie.tsv \
  --out gwasA_PHENO.filtered.tsv.gz \
  --emac-min 100 \
  --hwe-minp 1e-12 \
  --use-controls

# 3. Pre-munge for LDSC
python prep_munge.py \
  --in gwasA_PHENO.filtered.tsv.gz \
  --out gwasA_PHENO.pre_munged.tsv.gz \
  --logistic
```

Now the file `gwasA_PHENO.pre_munged.tsv.gz` is ready for LDSC's `munge_sumstats.py`.

---

## 6. Column Mapping

### Default Column Names (REGENIE format)

The scripts expect the following column names by default:

| Script | Argument | Default Value | Description |
|--------|----------|---------------|-------------|
| `post_QC.py` | `--id-col` | `Name` | Variant identifier |
| | `--aaf-col` | `AAF` | Effect allele frequency |
| | `--n-col` | `Num_Cases` | Sample size |
| | `--info-col` | `Info` | Key-value annotations |
| `prep_munge.py` | `--id-col` | `ID` | Variant identifier |
| | `--a1-col` | `A1` | Effect allele |
| | `--a2-col` | `A2` | Reference allele |
| | `--n-col` | `N` | Sample size |
| | `--p-col` | `P` | P-value |
| | `--effect-col` | `Effect` | Effect size (OR or BETA) |
| | `--info-col` | `INFO` | Key-value annotations |

### Expected File Format

Example REGENIE output header:
```
Name	Chr	Pos	Ref	Alt	Trait	Cohort	Model	Effect	LCI_Effect	UCI_Effect	Pval	AAF	Num_Cases	Cases_Ref	Cases_Het	Cases_Alt	Num_Controls	Controls_Ref	Controls_Het	Controls_Alt	Info
```

### INFO Field Format

The INFO field should contain semicolon-separated key=value pairs:

```
REGENIE_BETA=-0.011862;REGENIE_SE=0.009675;INFO=0.466699;MAC=175413.733333
```

Common keys:
- `MAC` - Minor Allele Count
- `EMAC` - Effective Minor Allele Count
- `HWE` - Hardy-Weinberg p-value
- `BETA` - Effect size (log-odds for logistic)
- `SE` - Standard error
- `LOGOR` - Log odds ratio
- `OR` - Odds ratio
- `CI95L`, `CI95U` - 95% confidence interval bounds
- `N_HOMREF`, `N_HET`, `N_HOMALT` - Genotype counts

If your files use different column names, adjust with the `--*-col` options.

---

## 7. Recommended Thresholds

| Analysis Type | EMAC min | HWE min p | Notes |
|--------------|----------|-----------|-------|
| **Standard GWAS** | 100 | 1e-12 | Default conservative filtering |
| **Strict QC** | 200-500 | 1e-10 | For high-quality analyses |
| **Rare variant analysis** | 50 | 1e-15 | More permissive for low frequency |
| **Common variants only** | 500 | 1e-10 | Focus on well-powered variants |

### Best Practices

1. **Always use `--use-controls` for case-control studies** - HWE should be tested in controls only
2. **Check output statistics** - Review filtering percentages to ensure reasonable QC
3. **Verify column names** - Use `head -1` to check your file's header before running
4. **Use compressed files** - Save disk space with `.gz` extension
5. **Keep filtering logs** - Redirect stderr to a log file for documentation

---

## 8. Troubleshooting

### No variants are filtered

**Cause**: Column names don't match the defaults

**Solution**: Check your file header and specify correct column names:
```bash
# View header
gunzip -c file.tsv.gz | head -1

# Adjust column names
python post_QC.py --id-col <your_id_col> --aaf-col <your_freq_col> ...
```

### All variants are filtered

**Cause**: Thresholds too strict or calculation issues

**Solution**: 
1. Check the first 3 variants in stderr output
2. Try more permissive thresholds: `--emac-min 10 --hwe-minp 1e-15`
3. Verify columns are read correctly

### EMAC values are None

**Cause**: AAF or N column missing or incorrectly specified

**Solution**: Use `--use-mac-from-info` if MAC is available in INFO field

### HWE p-values are None

**Cause**: No genotype count columns or HWE field in INFO

**Solution**: If HWE data unavailable, disable HWE filtering: `--hwe-minp 0`

---

## 9. Performance Notes

- **Speed**: ~100-500K variants/second (hardware dependent)
- **Memory**: Line-by-line processing, minimal memory footprint
- **Compression**: Transparent gzip read/write support

---

## 10. License / Provenance

- These scripts are a Python reimplementation of C/C++ utilities (`post_filter.C`, `prep_munge.C`) and mimic PLINK2 routines (HWE test, chi-square).
- Original PLINK2 headers (`plink2_base`, `plink2_string`, `plink2_stats`) are under permissive licenses (Boost, BSD, LGPL). This Python code is a clean rewrite without reuse of C/C++ source.
- HWE test algorithm based on Wigginton et al. (2005) "A Note on Exact Tests of Hardy-Weinberg Equilibrium"

---

## 11. Citation

If you use these scripts in your research, please cite:

- **REGENIE**: Mbatchou et al. (2021) Nature Genetics
- **LDSC**: Bulik-Sullivan et al. (2015) Nature Genetics
- **HWE test**: Wigginton et al. (2005) AJHG
