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

#### Custom column names (for different file formats)
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

### Overview

`prep_munge.py` reformats the filtered GWAS results into the columns expected by LDSC (`SNP, A1, A2, N, Z, P, BETA, SE`).

### Basic Usage

```bash
python prep_munge.py \
  --in gwasA_PHENO.filtered.tsv.gz \
  --out gwasA_PHENO.pre_munged.tsv.gz \
  --id-col Name \
  --a1-col Alt \
  --a2-col Ref \
  --n-col Num_Cases \
  --p-col Pval \
  --effect-col Effect \
  --info-col Info \
  --logistic
```

### All Available Arguments

#### 📁 Input/Output Files

| Argument | Default | Description |
|----------|---------|-------------|
| `--in` | `-` (stdin) | Input TSV or TSV.GZ file |
| `--out` | `-` (stdout) | Output TSV or TSV.GZ file |

#### 📊 Column Names

| Argument | Default | Description |
|----------|---------|-------------|
| `--id-col` | `ID` | SNP identifier column |
| `--a1-col` | `A1` | Effect allele column |
| `--a2-col` | `A2` | Reference/other allele column |
| `--n-col` | `N` | Sample size column |
| `--p-col` | `P` | P-value column |
| `--effect-col` | `Effect` | Effect size column (OR or BETA) |
| `--info-col` | `INFO` | INFO field (key=val;...) |

#### ⚙️ Model Options

| Argument | Type | Description |
|----------|------|-------------|
| `--logistic` | Flag | Specify if GWAS was logistic regression (Effect is OR) |

### Column Mapping for REGENIE Output

For standard REGENIE output format, use these column mappings:

```bash
python prep_munge.py \
  --in gwasA_PHENO.filtered.tsv.gz \
  --out gwasA_PHENO.pre_munged.tsv.gz \
  --id-col Name \
  --a1-col Alt \
  --a2-col Ref \
  --n-col Num_Cases \
  --p-col Pval \
  --effect-col Effect \
  --info-col Info \
  --logistic
```

### Output Format

The output file contains 8 columns in LDSC-compatible format:

```
SNP	A1	A2	N	Z	P	BETA	SE
rs367896724	AC	A	82776	-2.5177	0.22018	-0.011862	0.004710
rs201106462	TA	T	82776	-11.1164	0.373242	-0.008858	0.000797
rs534229142	A	G	82776	-0.3567	0.696678	-0.051285	0.143763
```

**Column descriptions:**
- **SNP**: Variant identifier (rsID or chr:pos:ref:alt)
- **A1**: Effect allele
- **A2**: Reference allele
- **N**: Sample size
- **Z**: Z-score (calculated as BETA/SE)
- **P**: P-value
- **BETA**: Effect size (log-odds for logistic regression)
- **SE**: Standard error

### Logic for BETA/SE Extraction

The script uses multiple strategies to extract BETA and SE from the INFO field:

#### Priority 1: Direct BETA and SE
```
INFO: REGENIE_BETA=-0.011862;REGENIE_SE=0.009675;SE=0.004710;...
→ BETA = -0.011862, SE = 0.004710
```
If both `BETA` and `SE` exist in INFO, use them directly.

#### Priority 2: LOGOR and SE
```
INFO: LOGOR=-0.011862;SE=0.009675;...
→ BETA = -0.011862 (use LOGOR as BETA), SE = 0.009675
```
If `LOGOR` and `SE` exist, use LOGOR as BETA.

#### Priority 3: OR and Confidence Intervals
```
INFO: OR=0.994242;CI95L=0.985106;CI95U=1.00346;...
→ BETA = ln(0.994242) = -0.005767
→ SE = (ln(1.00346) - ln(0.985106)) / 3.92 = 0.004711
```
Calculate from odds ratio and 95% confidence intervals.

#### Priority 4: Linear Model (if --logistic not specified)
```
Effect column: 0.125
INFO: SE=0.025;...
→ BETA = 0.125, SE = 0.025
```
For linear regression, use Effect column directly as BETA.

### Z-score Calculation

Z-scores are computed as: **Z = BETA / SE**

This provides a standardized measure of effect size across the genome.

---

## 5. Complete Workflow Example

### For Case-Control GWAS (Logistic Regression)

```bash
# 1. Concatenate chromosomes
(
  zcat gwasA.chr1_PHENO.regenie.gz | head -1
  for chr in $(seq 1 22); do
    zcat gwasA.chr${chr}_PHENO.regenie.gz | tail -n +2
  done
) > gwasA_PHENO.regenie.tsv

# 2. Post-QC filter (recommended settings for case-control)
python post_QC.py \
  --in gwasA_PHENO.regenie.tsv \
  --out gwasA_PHENO.filtered.tsv.gz \
  --emac-min 100 \
  --hwe-minp 1e-12 \
  --use-controls

# 3. Pre-munge for LDSC with correct column mappings
python prep_munge.py \
  --in gwasA_PHENO.filtered.tsv.gz \
  --out gwasA_PHENO.pre_munged.tsv.gz \
  --id-col Name \
  --a1-col Alt \
  --a2-col Ref \
  --n-col Num_Cases \
  --p-col Pval \
  --effect-col Effect \
  --info-col Info \
  --logistic

# 4. Verify output format
zcat gwasA_PHENO.pre_munged.tsv.gz | head -5

# 5. Run LDSC munge_sumstats
python ~/ldsc/munge_sumstats.py \
  --sumstats gwasA_PHENO.pre_munged.tsv.gz \
  --out gwasA_PHENO.munged \
  --merge-alleles ~/ldsc/w_hm3.snplist \
  --N 82776
```

### For Quantitative Traits (Linear Regression)

```bash
# Steps 1-2: Same as above

# 3. Pre-munge WITHOUT --logistic flag
python prep_munge.py \
  --in gwasA_PHENO.filtered.tsv.gz \
  --out gwasA_PHENO.pre_munged.tsv.gz \
  --id-col Name \
  --a1-col Alt \
  --a2-col Ref \
  --n-col Num_Cases \
  --p-col Pval \
  --effect-col Effect \
  --info-col Info
```

---

## 6. Column Mapping

### Default Column Names

The scripts expect different default column names:

| Script | Argument | Default Value | Typical REGENIE Column | Description |
|--------|----------|---------------|------------------------|-------------|
| `post_QC.py` | `--id-col` | `Name` | `Name` | Variant identifier |
| | `--aaf-col` | `AAF` | `AAF` | Effect allele frequency |
| | `--n-col` | `Num_Cases` | `Num_Cases` | Sample size |
| | `--info-col` | `Info` | `Info` | Key-value annotations |
| `prep_munge.py` | `--id-col` | `ID` | `Name` | Variant identifier |
| | `--a1-col` | `A1` | `Alt` | Effect allele |
| | `--a2-col` | `A2` | `Ref` | Reference allele |
| | `--n-col` | `N` | `Num_Cases` | Sample size |
| | `--p-col` | `P` | `Pval` | P-value |
| | `--effect-col` | `Effect` | `Effect` | Effect size (OR or BETA) |
| | `--info-col` | `INFO` | `Info` | Key-value annotations |

### Expected REGENIE File Format

Standard REGENIE output header:
```
Name	Chr	Pos	Ref	Alt	Trait	Cohort	Model	Effect	LCI_Effect	UCI_Effect	Pval	AAF	Num_Cases	Cases_Ref	Cases_Het	Cases_Alt	Num_Controls	Controls_Ref	Controls_Het	Controls_Alt	Info
```

**Key columns:**
- `Name`: Variant ID (rsID or chr:pos:ref:alt format)
- `Ref`: Reference allele
- `Alt`: Alternate/effect allele  
- `Effect`: Odds ratio (logistic) or beta (linear)
- `Pval`: P-value
- `AAF`: Alternate allele frequency
- `Num_Cases`: Number of cases
- `Num_Controls`: Number of controls
- `Cases_Ref/Het/Alt`: Genotype counts in cases
- `Controls_Ref/Het/Alt`: Genotype counts in controls
- `Info`: Semicolon-separated key=value pairs

### INFO Field Format

The INFO field contains semicolon-separated key=value pairs:

```
REGENIE_BETA=-0.011862;REGENIE_SE=0.009675;SE=0.004710;INFO=0.466699;MAC=175413.733333
```

**Common INFO keys:**
- `MAC` - Minor Allele Count
- `EMAC` - Effective Minor Allele Count
- `HWE` - Hardy-Weinberg equilibrium p-value
- `REGENIE_BETA` - Beta coefficient from REGENIE
- `REGENIE_SE` - Standard error from REGENIE
- `BETA` - Effect size (log-odds for logistic)
- `SE` - Standard error
- `LOGOR` - Log odds ratio
- `OR` - Odds ratio
- `CI95L`, `CI95U` - 95% confidence interval lower/upper bounds
- `N_HOMREF`, `N_HET`, `N_HOMALT` - Genotype counts
- `INFO` - Imputation info score

### Sample Size Considerations

**For case-control studies:**
- `Num_Cases`: Number of cases (e.g., 82776)
- `Num_Controls`: Number of controls (e.g., 138035)
- **Total N**: Cases + Controls (e.g., 220811)

**Recommendation**: Use `Num_Cases` for `--n-col` in `prep_munge.py`, then specify total N in LDSC's `munge_sumstats.py` using `--N` flag.

Alternatively, create a total N column:
```bash
zcat gwasA_PHENO.filtered.tsv.gz | awk 'BEGIN{FS=OFS="\t"} 
  NR==1 {print $0, "N_total"; next} 
  {print $0, $14+$18}' | gzip > gwasA_PHENO.filtered_with_N.tsv.gz
```

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
6. **Specify correct model type** - Use `--logistic` flag for case-control GWAS
7. **Match column names to your data** - REGENIE uses `Name`, `Alt`, `Ref`, `Pval`, `Info`

---

## 8. Troubleshooting

### post_QC.py Issues

#### No variants are filtered

**Cause**: Column names don't match the defaults

**Solution**: Check your file header and specify correct column names:
```bash
# View header
gunzip -c file.tsv.gz | head -1

# Adjust column names
python post_QC.py --id-col Name --aaf-col AAF --n-col Num_Cases --info-col Info ...
```

#### All variants are filtered

**Cause**: Thresholds too strict or calculation issues

**Solution**: 
1. Check the first 3 variants in stderr output
2. Try more permissive thresholds: `--emac-min 10 --hwe-minp 1e-15`
3. Verify columns are read correctly

#### EMAC values are None

**Cause**: AAF or N column missing or incorrectly specified

**Solution**: 
- Verify column names match your file
- Use `--use-mac-from-info` if MAC is available in INFO field

#### HWE p-values are None

**Cause**: No genotype count columns or HWE field in INFO

**Solution**: If HWE data unavailable, disable HWE filtering: `--hwe-minp 0`

### prep_munge.py Issues

#### Missing BETA or SE values

**Cause**: INFO field doesn't contain required keys or wrong column specified

**Solution**:
1. Check INFO field format: `zcat file.gz | cut -f22 | head -5`
2. Verify INFO contains `REGENIE_BETA`, `SE`, or `OR` with `CI95L/CI95U`
3. Ensure `--info-col Info` matches your file (not `INFO`)

#### Wrong allele coding

**Cause**: A1/A2 columns swapped

**Solution**: 
- For REGENIE: `--a1-col Alt --a2-col Ref`
- A1 should be the effect allele (usually Alt)
- A2 should be the reference allele (usually Ref)

#### Column not found errors

**Cause**: Default column names don't match REGENIE format

**Solution**: Always specify REGENIE columns explicitly:
```bash
python prep_munge.py \
  --id-col Name \
  --a1-col Alt \
  --a2-col Ref \
  --n-col Num_Cases \
  --p-col Pval \
  --info-col Info \
  --logistic
```

---

## 9. Performance Notes

- **Speed**: ~100-500K variants/second (hardware dependent)
- **Memory**: Line-by-line processing, minimal memory footprint (~50MB)
- **Compression**: Transparent gzip read/write support
- **File sizes**: 
  - Raw REGENIE: ~5-10GB uncompressed
  - Filtered: ~3-7GB compressed
  - Pre-munged: ~1-2GB compressed

---

## 10. Validation

### Verify post_QC.py output

```bash
# Check filtering statistics
python post_QC.py --in input.tsv --out output.tsv.gz 2>&1 | grep -A 10 "Statistiques"

# Compare variant counts
echo "Original: $(zcat input.tsv.gz | wc -l)"
echo "Filtered: $(zcat output.tsv.gz | wc -l)"
```

### Verify prep_munge.py output

```bash
# Check output format
zcat gwasA_PHENO.pre_munged.tsv.gz | head -1
# Should show: SNP	A1	A2	N	Z	P	BETA	SE

# Check for missing values
zcat gwasA_PHENO.pre_munged.tsv.gz | awk -F'\t' 'NR>1 && ($7=="" || $8=="")' | wc -l
# Should be low or zero

# Verify Z-scores are reasonable
zcat gwasA_PHENO.pre_munged.tsv.gz | awk -F'\t' 'NR>1 && $5!="" {print $5}' | sort -n | head -20
```

---

## 11. License / Provenance

- These scripts are a Python reimplementation of C/C++ utilities (`post_filter.C`, `prep_munge.C`) and mimic PLINK2 routines (HWE test, chi-square).
- Original PLINK2 headers (`plink2_base`, `plink2_string`, `plink2_stats`) are under permissive licenses (Boost, BSD, LGPL). This Python code is a clean rewrite without reuse of C/C++ source.
- HWE test algorithm based on Wigginton et al. (2005) "A Note on Exact Tests of Hardy-Weinberg Equilibrium"

---

## 12. Citation

If you use these scripts in your research, please cite:

- **REGENIE**: Mbatchou et al. (2021) "Computationally efficient whole-genome regression for quantitative and binary traits" *Nature Genetics*
- **LDSC**: Bulik-Sullivan et al. (2015) "LD Score regression distinguishes confounding from polygenicity in genome-wide association studies" *Nature Genetics*
- **HWE test**: Wigginton et al. (2005) "A Note on Exact Tests of Hardy-Weinberg Equilibrium" *American Journal of Human Genetics*

---

## 13. Support

For issues, questions, or contributions:
1. Check the troubleshooting section
2. Verify your command matches the examples
3. Include sample input data and error messages when reporting issues

**Common workflow issues:**
- Column name mismatches → Use explicit `--*-col` arguments
- Missing values → Check INFO field format and extraction logic
- Performance → Use compressed files and consider chunking very large datasets
