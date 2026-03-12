# 🐉 dragonSMASH

A resumable, multi-stage bacterial genome analysis pipeline for Oxford Nanopore assemblies. Integrates assembly, QC, annotation, biosynthetic gene cluster (BGC) detection, 16S taxonomy, and whole-genome taxonomy into a single bash script with a clean HTML summary report.

**dragonSMASH** | **v1.25.0** | Requires: conda, GNU parallel

---

## Pipeline overview

```
ONT reads (FASTQ)  ──►  01 Dragonflye      Assembly
                    ──►  02 SeqKit Stats    Assembly & read QC
                    ──►  03 Bakta           Gene annotation
                    ──►  04 antiSMASH       BGC detection
                    ──►  05 16S ID          Species ID (barrnap + Kraken2/SILVA)
                    ──►  06 GTDBTk          Whole-genome taxonomy
                    ──►  07 Summary         TSV + interactive HTML report
```

Existing assemblies (FASTA) can be supplied directly, skipping Stage 1. Any stage can be independently disabled with a `--skip-*` flag. Each stage is **resume-safe** — re-running the pipeline skips completed samples automatically.

---

## Requirements

### Tested versions

| Tool | Version |
|---|---|
| Dragonflye | 1.0.13 |
| SeqKit | 2.8.2 |
| Bakta | 1.12.0 |
| antiSMASH | 8.0.4 |
| GTDBTk | 2.6.1 |
| barrnap | 0.9 |
| Kraken2 | 2.17.1 |
| GNU parallel | any recent |

Later versions may work but have not been tested. antiSMASH in particular has a history of breaking changes between major versions — stick to 8.x if possible.

### Software (conda environments)

#### MDU server users (marvin / zaphod / trillian)

All environments are available in Liam Sharkey's conda prefix. Point your `dragonsmash.conf` directly at these — no installation needed:

```
/home/lsharkey/miniconda3/envs/dragonflye
/home/lsharkey/miniconda3/envs/seqkit
/home/lsharkey/miniconda3/envs/bakta
/home/lsharkey/miniconda3/envs/antismash8
/home/lsharkey/miniconda3/envs/16S_ID_env
/home/shared/conda/envs/gtdbtk              ← shared GTDBTk env (all servers)
```

Run `setup_dragonsmash.bash` once to create your own `16S_ID_env` and verify all paths are accessible.

#### External users

Create your own environments with the pinned versions:

```bash
conda create -n dragonflye -c bioconda dragonflye=1.0.13 -y
conda create -n seqkit     -c bioconda seqkit=2.8.2 -y
conda create -n bakta      -c bioconda bakta=1.12.0 -y
conda create -n 16S_ID_env -c bioconda -c conda-forge barrnap=0.9 kraken2=2.17.1 -y
# antiSMASH 8 — see https://download.antismash.secondarymetabolites.org/
# GTDBTk — see https://ecogenomics.github.io/GTDBTk/installing/index.html
```

GNU parallel is also required:
```bash
conda install -c conda-forge parallel
```

### Databases

| Database | Used by | Notes |
|---|---|---|
| Bakta DB (full) | Stage 3 | Per-user. Download: `conda run -n bakta bakta_db download --output ~/dbs/bakta/ --type full` |
| GTDBTk release226 | Stage 6 | Shared at `/home/shared/db/gtdbtk/release226` on MDU servers |
| SILVA (Kraken2) | Stage 5 | Shared at `/home/shared/db/kraken2/silva` on MDU servers |

---

## Installation

```bash
# 1. Clone or copy dragonSMASHv3.bash and setup_dragonsmash.bash
#    into a directory of your choice, e.g.:
mkdir -p ~/Projects/dragonsmash/scripts
cd ~/Projects/dragonsmash/scripts

# 2. Run setup (creates 16S_ID_env, symlink ~/bin/dragonSMASH, checks all deps)
bash setup_dragonsmash.bash

# 3. Source your shell rc file to pick up the PATH update
source ~/.bashrc

# 4. Copy dragonsmash.conf.example to dragonsmash.conf and edit BAKTA_DB
cp dragonsmash.conf.example dragonsmash.conf
```

The setup script will check all conda environments, databases, and tools, and print a personalised `dragonsmash.conf` block at the end.

---

## Configuration

Place a `dragonsmash.conf` file in the **same directory as `dragonSMASHv3.bash`**. All settings have sensible defaults; only `BAKTA_DB` is required.

```bash
# dragonsmash.conf

# Required — path to your Bakta database
# MDU users: use the path below (or download your own)
BAKTA_DB=/home/lsharkey/dbs/bakta/v6/db

# Conda environment names OR full paths
# MDU users: use full paths to point at existing envs (no installation needed)
ENV_NAME_DRAGONFLYE=/home/lsharkey/miniconda3/envs/dragonflye
ENV_NAME_SEQKIT=/home/lsharkey/miniconda3/envs/seqkit
ENV_NAME_BAKTA=/home/lsharkey/miniconda3/envs/bakta
ENV_NAME_ANTISMASH=/home/lsharkey/miniconda3/envs/antismash8
ENV_NAME_BARRNAP=/home/lsharkey/miniconda3/envs/16S_ID_env
# External users: use short names (e.g. ENV_NAME_DRAGONFLYE=dragonflye)

# GTDBTk — shared env, same path on all MDU servers
GTDBTK_CONDA_PREFIX=/home/shared/conda/envs/gtdbtk
GTDBTK_DATA_PATH=/home/shared/db/gtdbtk/release226
GTDBTK_SKANI_SKETCH_DIR=/home/shared/db/gtdbtk/release226/skani/sketches

# SILVA database for 16S classification
SILVA_DB=/home/shared/db/kraken2/silva

# antiSMASH options
ANTISMASH_TAXON=bacteria          # or: fungi
ANTISMASH_EXTRA_ARGS=--cc-mibig --cb-general --cb-subcluster --cb-knownclusters --rre

# barrnap options
BARRNAP_KINGDOM=bac               # or: arc, euk
BARRNAP_MIN_SCORE=0.8

# Watchdog (for long antiSMASH/Bakta runs — see --watchdog flag)
WATCHDOG_INTERVAL_SECS=600
WATCHDOG_STALL_CHECKS=3
```

---

## Usage

### Standard mode (assemble from ONT reads)

Requires an `inputs.txt` file with one sample name per line (no extension). Sample names must match the FASTQ filenames in `<ont_reads_dir>`.

```
inputs.txt:
  DMG2600585
  DMG2600586
  DMG2600587
```

```bash
dragonSMASH [options] <inputs.txt> <ont_reads_dir> <output_dir>
```

### Assembly-free mode (existing FASTAs)

Sample names are inferred from FASTA filenames. Skips Stage 1.

```bash
dragonSMASH [options] --assemblies <fasta_dir> <output_dir>
```

---

## Options

```
--cpus N              CPUs to use across all stages (default: all available).
                      Bakta and antiSMASH run parallel jobs; total CPU use
                      is capped at N. Override ANTISMASH_JOBS / BAKTA_JOBS
                      in dragonsmash.conf to control parallelism.

--max-contigs N       Exclude samples with > N contigs from annotation
                      (Bakta + antiSMASH). Default: 100. Fragmented assemblies
                      yield unreliable BGC predictions. Set to 0 to disable.

--assemblies <dir>    Use existing FASTA files; skip Stage 1 (assembly).

--illumina <dir>      Provide Illumina reads for Dragonflye polishing.
                      Directory should contain per-sample subdirectories.

--skip-bakta          Skip Stage 3 (Bakta annotation). Also skips antiSMASH,
                      which requires Bakta output.
--skip-antismash      Skip Stage 4 (antiSMASH) only.
--skip-16s            Skip Stage 5 (barrnap + Kraken2/SILVA 16S ID).
--skip-gtdbtk         Skip Stage 6 (GTDBTk taxonomy).

--watchdog            Enable a hang-detection watchdog for Stages 3–4.
                      Kills stalled diamond jobs after N consecutive checks
                      with no CPU progress. Configure interval and check count
                      in dragonsmash.conf. Off by default.

--dry-run             Print all commands without executing anything.
--min-disk N          Minimum free disk space required in GB (default: 50).
--help, -h            Show usage.
```

---

## Examples

```bash
# Standard run — assemble 15 samples, all stages, 20 CPUs
dragonSMASH --cpus 20 inputs.txt reads/ out/

# Existing assemblies, skip antiSMASH, run in screen
screen -S dragonsmasher \
  dragonSMASH --assemblies fastas/ --skip-antismash --cpus 16 out/

# Large dataset — existing assemblies, watchdog enabled
screen -S nocardia \
  dragonSMASH \
    --assemblies ~/Projects/Nocardia/Best_assemblies/ \
    --cpus 4 \
    --watchdog \
    --max-contigs 500 \
    ~/Projects/Nocardia/dragonsmash_out/

# Taxonomy only (GTDBTk + 16S, no annotation)
dragonSMASH --assemblies fastas/ --skip-bakta --skip-16s --cpus 8 out/
```

---

## Outputs

```
output_dir/
  01_assembly/<sample>/
      final.contigs.fa          Assembled contigs
  02_seqkit_stats/
      assembly_stats.tsv        Per-sample assembly statistics
      read_stats.tsv            Per-sample read statistics (standard mode)
  03_bakta/<sample>/
      <sample>.gbff             GenBank flat file (input to antiSMASH)
      <sample>.tsv              Gene annotation table
  04_antismash/<sample>/
      index.html                antiSMASH interactive report
      *.region*.gbk             Per-region GenBank files
  05_16S/<sample>/
      16S.fasta                 Extracted 16S rRNA sequences
      kraken2_report.txt        Kraken2 classification report
      kraken2_output.txt        Per-sequence Kraken2 output
      classification.tsv        Species, genus, family, copy count, mixed flag
  06_gtdbtk/
      gtdbtk.bac120.summary.tsv GTDB taxonomy for all samples
  07_summary/
      summary.tsv               Master summary table (all samples × all metrics)
      summary.html              Interactive HTML report
```

### Summary TSV columns

| Column | Description |
|---|---|
| `sample` | Sample name |
| `gtdb_classification` | Full GTDB lineage string |
| `gtdb_genus` / `gtdb_species` | Parsed GTDB genus and species |
| `16S_species` / `16S_genus` / `16S_family` | Top 16S hit (majority vote) |
| `16S_copies` | Number of 16S rRNA gene copies detected |
| `16S_n_genera` | Number of distinct genera across all copies |
| `16S_mixed_flag` | `OK`, `genus_conflict`, `family_conflict`, or `unclassified` |
| `num_contigs` | Contig count |
| `assembly_size_bp` | Total assembly length |
| `N50` / `N90` / `L50` | Standard assembly statistics |
| `bakta_CDS` | Number of predicted coding sequences |
| `antismash_BGCs_total` | Total BGC regions detected |
| `antismash_BGCs_on_edge` | BGCs truncated by contig edges |
| `NRPS`, `PKS-I`, etc. | Per-type BGC counts (dynamic columns) |

### Mixed culture detection

The `16S_mixed_flag` field reports the consistency of 16S rRNA gene copies within an assembly:

| Flag | Meaning |
|---|---|
| `OK` | All copies agree at genus level — consistent with a pure culture |
| `genus_conflict` | Different genera, same family — may be intra-genomic variation or a closely related mixed culture |
| `family_conflict` | Multiple families — strong indicator of mixed culture or contamination |
| `unclassified` | No copies matched the SILVA database |

---

## Resume behaviour

dragonSMASH is designed to be safely re-run. Each stage checks for existing output before running:

- **Stages 1–5**: per-sample — a completed sample is skipped, a partial one is re-run
- **Stage 6 (GTDBTk)**: skipped if `gtdbtk.bac120.summary.tsv` already exists
- **Stage 7 (summary)**: always regenerated from existing stage outputs

To force a stage to re-run, delete its output directory:

```bash
# Re-run antiSMASH for all samples
rm -rf out/04_antismash/

# Re-run summary only
rm -rf out/07_summary/
```

---

## Troubleshooting

**antiSMASH hangs indefinitely**
Use `--watchdog`. antiSMASH runs diamond internally, which can stall on some assemblies. The watchdog monitors CPU activity and kills stalled jobs.

**"No assembly found for sample"**
Check that FASTA filenames match sample names in `inputs.txt` exactly (no extension). Supported extensions: `.fasta`, `.fa`, `.fna`, `.fas`.

**GTDBTk fails with database errors**
Ensure `GTDBTK_DATA_PATH` and `GTDBTK_SKANI_SKETCH_DIR` in `dragonsmash.conf` point to a GTDBTk release226 database. GTDBTk ≥ 2.3 is required for skani support.

**Sample filtered due to contig count**
Samples exceeding `--max-contigs` are excluded from Bakta and antiSMASH but still processed by 16S ID and GTDBTk. Increase `--max-contigs` or use `--max-contigs 0` to annotate all samples regardless.

**BGC type columns empty in summary**
Ensure you are using antiSMASH 6 or later (dragonSMASH parses the `/product` qualifier from region GBK files, which is the antiSMASH 6+ format).

---

## Citation

If you use dragonSMASH in your work, please cite the underlying tools:

- **Dragonflye**: [github.com/rpetit3/dragonflye](https://github.com/rpetit3/dragonflye)
- **SeqKit**: Shen et al. (2016) *PLOS ONE* doi:10.1371/journal.pone.0163962
- **Bakta**: Schwengers et al. (2021) *Microbial Genomics* doi:10.1099/mgen.0.000685
- **antiSMASH**: Blin et al. (2023) *Nucleic Acids Research* doi:10.1093/nar/gkad344
- **barrnap**: [github.com/tseemann/barrnap](https://github.com/tseemann/barrnap)
- **Kraken2**: Wood et al. (2019) *Genome Biology* doi:10.1186/s13059-019-1891-0
- **GTDBTk**: Chaumeil et al. (2022) *Bioinformatics* doi:10.1093/bioinformatics/btac672
- **GTDB**: Parks et al. (2022) *Nature Biotechnology* doi:10.1038/s41587-020-0501-8

---

## Acknowledgements

Developed at the [Melbourne Infectious Diseases Unit (MDU)](https://mdu.mdhs.unimelb.edu.au/), University of Melbourne.
