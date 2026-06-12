# AlphaGenome Variant Analysis Tool

This tool predicts the effect of genetic variants on gene regulation using the AlphaGenome model developed by Google DeepMind. Given a list of variants (mutations), it produces scores and visualisations showing how each variant is predicted to affect gene expression, chromatin accessibility, splicing, and more.

No programming knowledge is required to use this tool.

---

## Requirements

- A computer
- An internet connection
- An AlphaGenome API key (request one at https://deepmind.google.com/science/alphagenome/)

That is all. The tool installs everything else automatically on first run.

---

## Getting started

Open a Terminal and navigate to the folder containing `run_alphagenome.sh` and `alphagenome_pipeline.py`.

Then run `chmod +x run_alphagenome.sh` to be able to run this script.

The first time you run the tool it will install a small package manager called `uv` and all required software. This takes about 2–5 minutes and only happens once.

---

## Preparing your input file

Create a spreadsheet in Excel or a plain text file (.csv) with the following columns. Save it as a CSV file (File → Save As → CSV UTF-8).

| Column | Required | Description |
|---|---|---|
| ID | Yes | A unique name for the variant (e.g. `rs80356730` or `TARDBP_M337V`) |
| CHROM | Yes | Chromosome (e.g. `chr1`) |
| POS | Yes | Genomic position of the variant (number) |
| REF | Yes* | Reference base(s). Leave empty for deletions |
| ALT | Yes* | Alternate base(s). Leave empty for deletions |
| Gene | Yes | Gene symbol (e.g. `TARDBP`) |
| deletion_length | No | For deletions: number of deleted bases |
| Strandedness | No | Strand of the gene (`+` or `-`). Optional — the strand is taken from the gene annotation. If provided, it is checked against the annotation and a warning is shown if it disagrees |
| WINDOW_BP | No | Width of the plot window in base pairs. Default: 24000 |
| CENTRE_POS | No | Centre position for the plot — must be an actual genomic coordinate on the same chromosome (e.g. `11010000`), **not** a window size. Default: variant position |
| Phenotype | No | Clinical phenotype label (e.g. `ALS`) |
| cDNA | No | cDNA notation (e.g. `c.1009A>G`) |
| HGVS | No | Protein change notation (e.g. `p.(Met337Val)`) |

*For deletions, leave REF and ALT empty and fill in `deletion_length` instead.

### Example input file

```
ID,CHROM,POS,REF,ALT,Gene,Phenotype,cDNA,HGVS
sample_1,chr1,169549811,G,AAAA,F5,Thrombophilia,c.1601G>A,p.(Arg534Gln)
sample_2,chr6,26092913,G,A,HFE,Hereditary hemochromatosis,c.845G>A,p.(Cys282Tyr)
sample_3,chr1,11796321,G,,MTHFR,Hyperhomocysteinemia,c.665C>T,p.(Ala222Val)
```

---

## Running the tool

Open a Terminal in the folder containing `run_alphagenome.sh` and run:

```bash
./run_alphagenome.sh --api-key YOUR_API_KEY --input-file variants.csv --output-dir ./results
```

Replace `YOUR_API_KEY` with your AlphaGenome API key and `variants.csv` with the name of your input file.

### Examples

**Minimal — human variants, default settings:**
```bash
./run_alphagenome.sh --api-key YOUR_API_KEY --input-file variants.csv
```

**Save results to a specific folder:**
```bash
./run_alphagenome.sh --api-key YOUR_API_KEY --input-file variants.csv --output-dir ./TARDBP_results
```

**Mouse variants:**
```bash
./run_alphagenome.sh --api-key YOUR_API_KEY --input-file mouse_variants.csv --organism mouse --output-dir ./mouse_results
```

**Restrict outputs to specific tissues (e.g. brain and neurons):**
```bash
./run_alphagenome.sh --api-key YOUR_API_KEY --input-file variants.csv --ontology-terms brain_terms.txt --output-dir ./results
```

**Find all available ontology terms before running (open the output file in Excel to browse):**
```bash
./run_alphagenome.sh --api-key YOUR_API_KEY --organism human --list-ontologies
```

**Use a custom gene annotation file:**
```bash
./run_alphagenome.sh --api-key YOUR_API_KEY --input-file variants.csv --gtf /path/to/my_annotation.gtf.gz --output-dir ./results
```

**All options together — mouse variants, custom annotation, specific tissues, custom output folder:**
```bash
./run_alphagenome.sh \
  --api-key YOUR_API_KEY \
  --input-file mouse_variants.csv \
  --organism mouse \
  --ontology-terms mouse_brain_terms.txt \
  --gtf /path/to/gencode.vM35.annotation.gtf.gz \
  --output-dir ./mouse_brain_results
```

> The `\` at the end of each line is just a line break for readability. You can also write the whole command on a single line.

### All available options

```
./run_alphagenome.sh --api-key KEY --input-file FILE [options]

Required:
  --api-key KEY          Your AlphaGenome API key
  --input-file FILE      Path to your input CSV file

Optional:
  --output-dir DIR       Where to save results (default: ./results_alphagenome)
  --organism ORGANISM    human or mouse (default: human)
  --ontology-terms FILE  File with ontology terms to restrict outputs (one term per line)
  --gtf FILE             Custom gene annotation file (.feather or .gtf/.gtf.gz)
  --list-ontologies      Write all available ontology terms to a TSV file, then exit
```

---

## Output

Results are saved in the output directory, organised per variant:

```
results/
├── all_variants_scores.tsv                        ← Scores for all variants combined
├── VARIANT_ID/
│   ├── scores_VARIANT_ID.tsv                      ← Scores for this variant
│   └── plots/
│       ├── rna_seq_overlay_VARIANT.png            ← REF vs ALT track
│       ├── rna_seq_diff_VARIANT.png               ← ALT − REF difference
│       ├── cage_overlay_VARIANT.png
│       ├── cage_diff_VARIANT.png
│       ├── splice_junctions_overlay_VARIANT.png
│       ├── splice_junctions_diff_VARIANT.png
│       └── ...                                    ← same pattern for atac, dnase, chip_histone, chip_tf, splice_sites, splice_site_usage
```

Open `.tsv` files in Excel. Open `.png` files in any image viewer.

### What the score columns mean

Each row in a `scores.tsv` file is one prediction for one variant against one output track (e.g. RNA-seq in a particular tissue). The main columns:

| Column | Meaning |
|---|---|
| variant_id | The variant, e.g. `chr1:11022418:A>G` |
| scored_interval | The genomic region that was scored |
| gene_name | Gene the score relates to (for gene-based scores), e.g. `TARDBP` |
| gene_strand | Strand of that gene (`+` or `-`) |
| output_type | What was predicted, e.g. `RNA_SEQ`, `ATAC`, `SPLICE_SITES` |
| variant_scorer | The scoring method used |
| track_name | The specific data track, e.g. `UBERON:0000955 total RNA-seq` |
| ontology_curie | Tissue/cell-type code of the track (e.g. `UBERON:0000955`) |
| biosample_name | Human-readable tissue/cell type (e.g. `brain`) |
| **raw_score** | **The predicted effect of the variant** (further from 0 = larger predicted effect) |
| quantile_score | The raw score expressed as a percentile vs a reference set (human only; blank for mouse) |

The remaining columns (`input_gene`, `CHROM`, `POS`, `cDNA`, `HGVS`, `Phenotype`, `is_snv`, `is_deletion`, …) are copied from your input file so each score row is self-contained.

To find the variants/tissues with the strongest predicted effect, sort by `raw_score` (largest positive or most negative values).

---

## Finding ontology terms

Ontology terms define which cell types and tissues to include in the output. If you do not provide any, a set of default terms is used.

To see all available terms for your organism, run:

```bash
./run_alphagenome.sh --api-key YOUR_API_KEY --organism human --list-ontologies
```

This writes a file called `available_ontologies_human.tsv` (or `_mouse.tsv`) that you can open in Excel. Pick the rows relevant to your tissue of interest and copy the `ontology_curie` values (e.g. `UBERON:0000955`) into a plain text file, one per line. Then pass that file with `--ontology-terms`.

Alternatively you can go to https://www.alphagenomedocs.com/colabs/tissue_ontology_mapping.html run it in google colab (space ship icon top right).


### Example ontology terms file (`my_terms.txt`)

```
UBERON:0000955
CL:0000540
CL:0000100
```

```bash
./run_alphagenome.sh --api-key YOUR_API_KEY --input-file variants.csv --ontology-terms my_terms.txt --output-dir ./results
```

---

## Mouse variants

For mouse variants, pass `--organism mouse` on the command line. All variants in one file are treated as the same organism, so use one file per organism:

```bash
./run_alphagenome.sh --api-key YOUR_API_KEY --input-file mouse_variants.csv --organism mouse --output-dir ./results_mouse
```

---

## Troubleshooting

**`zsh: command not found: run_alphagenome.sh`**
Add `./` before the command: `./run_alphagenome.sh ...`

**`ERROR: Missing required columns`**
Check that your CSV has at minimum: `ID`, `CHROM`, `POS`, `REF`, `ALT`, `Gene`.

**`ERROR: No variants were successfully processed`**
Your API key may be incorrect. Double-check it at https://deepmind.google.com/science/alphagenome/

**`WARNING: plot interval falls outside the prediction region`**
Your `CENTRE_POS` or `WINDOW_BP` value is outside the gene region. The tool will fall back to a gene-centred window automatically.

**First run is very slow**
The gene annotation file (~100 MB) is downloaded once and cached. Subsequent runs start much faster.
