# Yeast - BLAST: XML Analysis Tool

A Python script for parsing and summarising NCBI BLAST XML results. It works with both **protein** (`blastp`) and **RNA/nucleotide** (`blastn`) searches, extracts alignment metrics, classifies hit quality, flags yeast-species hits, and writes TSV reports plus an interactive HTML summary.

---

## Features

- **Auto-detects blastp vs blastn** from the XML — no flags needed
- **Quality classification** with thresholds tuned per BLAST type (see [Quality thresholds](#quality-thresholds))
- **Yeast-species tracking**: identifies hits from *Saccharomyces*, *Candida*, *Kluyveromyces*, and other yeasts; shows the best yeast alternative when the top hit is non-yeast
- **HTML report** with colour-coded quality rows and embedded plots
- **TSV outputs** for easy downstream analysis in R, Excel, etc.

---

## Requirements

**Required**
```
biopython
```

**Optional** (needed for PNG plots)
```
numpy
matplotlib
```

Install everything at once:
```bash
pip install biopython numpy matplotlib
```

Or the minimal install (TSV outputs only, no plots):
```bash
pip install biopython
```

---

## Usage

### Basic run

```bash
python blast_v3.py --xml results.xml --out ./output
```

### With an e-value filter

Only HSPs with e-value ≤ 1×10⁻¹⁰ are considered:

```bash
python blast_v3.py --xml results.xml --out ./output --evalue 1e-10
```

### Control how many species appear in the bar chart

```bash
python blast_v3.py --xml results.xml --out ./output --top 15
```

### Full example — protein search

```bash
python blast_v3.py \
    --xml protein_blast.xml \
    --out results/protein \
    --evalue 1e-5 \
    --top 20
```

### Full example — RNA/nucleotide search

The script detects `blastn` automatically from the XML; no extra flags are needed.

```bash
python blast_v3.py \
    --xml RNA_blast.xml \
    --out results/rna
```

---

## Arguments

| Argument | Default | Description |
|---|---|---|
| `--xml` | *(required)* | Path to the BLAST XML output file |
| `--out` | *(required)* | Directory where all output files are written (created if it does not exist) |
| `--evalue` | `None` | E-value threshold for filtering HSPs. HSPs above this value are ignored. |
| `--top` | `20` | Number of top species to include in the species bar chart |

---

## Output files

All files are written to the directory specified with `--out`.

| File | Description |
|---|---|
| `blast_report.html` | Interactive HTML table with colour-coded quality rows and embedded plots |
| `top_hits.tsv` | One row per query with all raw metrics |
| `top_hits_summary.tsv` | Simplified view including quality label and best yeast alternative |
| `top_species.png` | Bar chart of the most frequent organisms *(requires matplotlib)* |
| `distributions.png` | Histograms of e-value, % identity, and coverage *(requires matplotlib)* |
| `species_counts.tsv` | Species frequency table *(fallback when matplotlib is absent)* |
| `evalues.tsv` | Raw e-values *(fallback)* |
| `percent_id.tsv` | Raw % identity values *(fallback)* |
| `coverage.tsv` | Raw coverage values *(fallback)* |

### Column reference for `top_hits.tsv`

| Column | Description |
|---|---|
| `Query` | Query sequence name from the XML |
| `Top_hit_title` | Full title of the best hit |
| `E_value` | Best e-value (minimum across all HSPs for blastn) |
| `Bit_score` | Highest bit-score across HSPs |
| `Identities` | Total identical positions (pooled across HSPs for blastn) |
| `Align_length` | Total alignment length (summed across HSPs for blastn) |
| `Percent_identity` | `Identities / Align_length × 100` |
| `Coverage_percent` | Union of aligned query spans / query length × 100 |
| `Organism` | Organism extracted from the hit title (or `Unknown` for PDB RNA entries) |
| `PDB_id` | 4-character PDB accession code |
| `Is_Yeast` | `Yes` / `No` |

---

## Quality thresholds

Hits are classified as **Good**, **OK**, or **Poor** based on the BLAST program type detected in the XML.

### Protein (`blastp`)

| Quality | % Identity | Coverage | E-value |
|---|---|---|---|
| Good | ≥ 50% | ≥ 70% | < 1×10⁻⁵ |
| OK | ≥ 30% | ≥ 50% | < 1×10⁻³ |
| Poor | anything else | | |

### RNA / nucleotide (`blastn`)

rRNA queries routinely cover only a portion of a much longer ribosomal RNA molecule, so the coverage threshold is intentionally lower. Identity, however, is expected to be high.

| Quality | % Identity | Coverage | E-value |
|---|---|---|---|
| Good | ≥ 85% | ≥ 20% | < 1×10⁻⁵ |
| OK | ≥ 70% | ≥ 10% | < 1×10⁻³ |
| Poor | anything else | | |

> **Note on E-value = 0**: An e-value reported as exactly `0` in the XML means the score is smaller than the smallest representable float. The script treats it as the best possible value (`1×10⁻³⁰⁰`) rather than as missing data.

---

## How coverage is calculated for RNA

NCBI BLAST nucleotide searches frequently return multiple HSPs (High-scoring Segment Pairs) per hit, especially for ribosomal RNA. A naive implementation that only looks at the first HSP will underestimate coverage substantially.

This script computes coverage as the **union of all aligned query coordinate ranges**, which matches the value shown on the NCBI web interface:

```
Coverage = (union of non-overlapping query spans across all HSPs) / query length × 100
```

Example from a 25S rRNA query (length 6509 nt, 7 HSPs):

| HSP | Query range | Aligned length | Identity |
|---|---|---|---|
| 1 | 2247–2773 | 533 | 490 |
| 2 | 702–1390 | 715 | 605 |
| 3 | 1803–2221 | 425 | 376 |
| 4 | 403–564 | 165 | 147 |
| 5 | 1680–1799 | 121 | 111 |
| 6 | 208–296 | 89 | 86 |
| 7 | 588–621 | 34 | 34 |

Union span = 2033 nt → **Coverage = 31.2%**, **% Identity = 88.81%** — .

---

## Yeast tracking

The HTML report and summary TSV highlight whether the best hit comes from a yeast organism. When it does not, the best available yeast hit is shown as a secondary entry in the same row.

Recognised yeast genera:

- *Saccharomyces*
- *Candida*
- *Pichia*
- *Komagataella*
- *Kluyveromyces*
- *Schizosaccharomyces*
- *Yarrowia*
- *Debaryomyces*
- *Hansenula*

---

## HTML report colour coding

| Colour | Meaning |
|---|---|
| 🟢 Green | Good hit (yeast) |
| 🟡 Yellow | OK hit (yeast) |
| 🔴 Red | Poor hit (yeast) |
| 🟣 Purple | Top hit is non-yeast (best yeast alternative shown below if available) |

---

## Using as a module

The script can be imported and called from other Python code:

```python
from blast_v3 import analyze, parse_blast_results, classify_hit_quality

# Run the full pipeline
analyze(
    xml_path="results.xml",
    out_dir="./output",
    evalue_threshold=1e-5,
    top_n_species=20
)
```

```python
# Parse only, then inspect rows manually
rows, species_counter, evalues, pct_ids, coverages = parse_blast_results("results.xml")

for row in rows:
    quality = classify_hit_quality(
        row['evalue'], row['pct_id'], row['coverage'], row['blast_type']
    )
    print(f"{row['query']}: {row['pct_id']}% id, {row['coverage']}% cov → {quality}")
```

---

## Getting BLAST XML output from NCBI

Run a search on the [NCBI BLAST web interface](https://blast.ncbi.nlm.nih.gov/), then download results in XML format:

**Web interface:** Download → XML (hit the *Download All* button and select *XML*)

**Command-line BLAST:**
```bash
# Protein search
blastp -query proteins.fasta -db pdb -outfmt 5 -out protein_blast.xml

# Nucleotide / RNA search
blastn -query rna_sequences.fasta -db pdbnt -outfmt 5 -out RNA_blast.xml
```

`-outfmt 5` is the flag for XML output.

---

## Project structure

```
.
├── blast_v2.py       # Main script
├── README.md         # This file
└── example/
    ├── protein_blast.xml
    └── RNA_blast.xml
```
