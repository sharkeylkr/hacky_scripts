#!/usr/bin/env bash
# ==============================================================================
# dragonSMASH v3
# Dragonflye → SeqKit Stats → Bakta → antiSMASH → GTDBTk
# ------------------------------------------------------------------------------
# Version : 1.7.0
# Date    : 11/03/2026
# Changes : Initial versioned release. ASCII progress bars, British date log
#           format, resume logic, CPU safety, final.contigs.fa symlink fix.
# ------------------------------------------------------------------------------
# Version : 1.2.0
# Date    : 07/03/2026
# Changes : Unset BAKTA_DB env var before conf load so environment cannot
#           override script default. Conf file still takes precedence.
# ------------------------------------------------------------------------------
# Version : 1.2.0
# Date    : 07/03/2026
# Changes : Fix Bakta flags to match installed version (remove --amrfinder-db).
#           Add interactive prompts for genus and gram stain (--genus, --gram).
#           override script default. Conf file still takes precedence.
# ------------------------------------------------------------------------------
# Version : 1.2.0
# Date    : 07/03/2026
# Changes : Add --skip-expert-system amrfinder to Bakta command. AMRFinderPlus
#           DB on shared server is outdated; skip until cwwalsh updates it.
#           override script default. Conf file still takes precedence.
#           format, resume logic, CPU safety (antiSMASH parallel jobs),
#           final.contigs.fa symlink fix for standard assembly mode.
# ------------------------------------------------------------------------------
# Version : 1.3.0
# Date    : 09/03/2026
# Changes : Parallelise Stage 3 (Bakta). BAKTA_JOBS controls concurrent samples
#           (default: 2); BAKTA_CPUS = CPUS / BAKTA_JOBS per job. Override in
#           dragonsmash.conf. Progress bar uses same counter-file pattern as S4.
# ------------------------------------------------------------------------------
# Version : 1.4.0
# Date    : 09/03/2026
# Changes : Add Stage 5 (GTDBTk classify_wf) for GTDB taxonomy. Runs on all
#           samples including filtered ones (no contig limit for taxonomy).
#           Adds gtdb_classification, gtdb_genus, gtdb_species columns to
#           summary TSV. Summary moved to 06_summary/summary.tsv.
# ------------------------------------------------------------------------------
# Version : 1.5.0
# Date    : 09/03/2026
# Changes : Remove all hardcoded server-specific paths. CONDA_BASE now
#           auto-detected via 'conda info --base'. BAKTA_DB required in conf
#           (no default). Conda env names configurable via ENV_NAME_* variables
#           in dragonsmash.conf. Warns if no conf file found. Adds
#           dragonsmash.conf.example with full documentation.
# ------------------------------------------------------------------------------
# Version : 1.6.0
# Date    : 11/03/2026
# Changes : Background watchdog monitors all diamond processes during Stages 3
#           and 4. Any diamond job with no output file progress for 2 hours is
#           killed and its partial sample output cleaned up so parallel can
#           restart it. Watchdog logs kills to the main log. Script renamed
#           dragonSMASHv3.bash.
# ------------------------------------------------------------------------------
# Version : 1.7.0
# Date    : 11/03/2026
# Changes : Stage 2 overhaul. Assembly stats now TSV with header-based column
#           parsing; adds N90 (-N 90) and L50 columns; drops meaningless Q20/Q30
#           from FASTA assemblies. New read_stats.tsv captures ONT (and
#           Illumina if provided) read metrics per sample (standard mode only).
#           Summary TSV gains N90 and L50 columns. All seqkit output now TSV.
# ------------------------------------------------------------------------------
# Version : 1.25.0
# Date    : 12/03/2026
# Changes : Add --skip-bakta, --skip-antismash, --skip-16s, --skip-gtdbtk flags
#           to allow any non-essential stage to be disabled. --skip-bakta
#           automatically implies --skip-antismash. Skipped stages print a
#           one-line notice and leave summary columns as NA.
# ==============================================================================

set -uo pipefail

# (banner printed after log setup below)

# Early stubs — redefined properly after log redirect is set up
tprint() { printf '%s\n' "$*"; }
eprint() { printf '%s\n' "$*" >&2; }

# ------------------------------------------------------------------------------
# Usage
# ------------------------------------------------------------------------------
usage() {
  cat << USAGE
Usage: $0 [options] <inputs.txt> <ont_reads_dir> <output_dir>
       $0 [options] --assemblies <fasta_dir> <output_dir>

Standard mode (with assembly):
  <inputs.txt>          One sample name per line (no extension).
  <ont_reads_dir>       Directory containing ONT reads (.fastq.gz or .fastq).
  <output_dir>          Root output directory for all pipeline stages.

Skip-assembly mode:
  --assemblies <dir>    Directory of existing FASTA files (.fasta, .fa, .fna, .fas).
                        Sample names are inferred from filenames. Skips Stage 1.
  <output_dir>          Root output directory for all pipeline stages.

Options:
  --cpus N              CPUs to use for all stages (default: all available via nproc).
                        Stages 1-3 use all N CPUs on one sample at a time.
                        Stage 4 (antiSMASH) runs 4 parallel jobs at N/4 CPUs each,
                        so total CPU use never exceeds N. Override ANTISMASH_JOBS in
                        dragonsmash.conf to change the number of parallel jobs.
  --max-contigs N       Skip annotation (Bakta + antiSMASH) for assemblies with more
                        than N contigs (default: 100). Highly fragmented assemblies
                        are unlikely to yield meaningful BGC predictions. Set to 0
                        to disable filtering.
  # (--force is not needed — dragonflye always overwrites partial output dirs)
  --illumina <dir>      Directory of Illumina reads in per-sample subfolders (optional).
  --dry-run             Print what would be run without executing anything.
  --min-disk N          Minimum free disk space required in GB (default: 50).
  --watchdog            Enable the diamond hang watchdog during Stages 3 and 4.
                        Kills diamond jobs with no CPU progress for 3 consecutive
                        check intervals (default: 30 min). Configure
                        WATCHDOG_INTERVAL_SECS and WATCHDOG_STALL_CHECKS in
                        dragonsmash.conf. Off by default.
  --skip-bakta          Skip Stage 3 (Bakta annotation). Also skips antiSMASH,
                        which requires Bakta output.
  --skip-antismash      Skip Stage 4 (antiSMASH BGC detection) only.
  --skip-16s            Skip Stage 5 (barrnap + Kraken2/SILVA 16S species ID).
  --skip-gtdbtk         Skip Stage 6 (GTDBTk taxonomy). Summary will show NA
                        for GTDB classification columns.
  --help, -h            Show this message.

Pipeline stages (each skipped if output already exists):
  1. Dragonflye   -> <output_dir>/01_assembly/<sample>/  [skipped with --assemblies]
  2. SeqKit Stats -> <output_dir>/02_seqkit_stats/assembly_stats.tsv
                              <output_dir>/02_seqkit_stats/read_stats.tsv  [standard mode only]
  3. Bakta        -> <output_dir>/03_bakta/<sample>/
  4. antiSMASH    -> <output_dir>/04_antismash/<sample>/
  5. 16S ID       -> <output_dir>/05_16S/<sample>/  (barrnap + Kraken2/SILVA)
  6. GTDBTk       -> <output_dir>/06_gtdbtk/

Summary TSV written to: <output_dir>/07_summary/summary.tsv
USAGE
  exit 1
}

# ------------------------------------------------------------------------------
# Parse args
# ------------------------------------------------------------------------------
CPUS=$(command -v nproc >/dev/null 2>&1 && nproc || echo 4)
ILLUMINA_READS_DIR=""
DRY_RUN=false
MIN_DISK_GB=50
EXISTING_ASSEMBLIES_DIR=""
MAX_CONTIGS=100
WATCHDOG_ENABLED=false
SKIP_BAKTA=false
SKIP_ANTISMASH=false
SKIP_16S=false
SKIP_GTDBTK=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cpus)           CPUS="$2";                      shift 2 ;;
    --illumina)       ILLUMINA_READS_DIR="$2";        shift 2 ;;
    --dry-run)        DRY_RUN=true;                   shift 1 ;;
    --min-disk)       MIN_DISK_GB="$2";               shift 2 ;;
    --assemblies)     EXISTING_ASSEMBLIES_DIR="$2";   shift 2 ;;
    --max-contigs)    MAX_CONTIGS="$2";               shift 2 ;;
    --watchdog)       WATCHDOG_ENABLED=true;          shift 1 ;;
    --skip-bakta)     SKIP_BAKTA=true;                shift 1 ;;
    --skip-antismash) SKIP_ANTISMASH=true;            shift 1 ;;
    --skip-16s)       SKIP_16S=true;                  shift 1 ;;
    --skip-gtdbtk)    SKIP_GTDBTK=true;               shift 1 ;;
    --help|-h)        usage ;;
    *)             break ;;
  esac
done

if [[ -n "$EXISTING_ASSEMBLIES_DIR" ]]; then
  # Skip-assembly mode: only output_dir is a positional arg
  if [[ $# -lt 1 ]]; then usage; fi
  INPUTS_FILE=""
  ONT_READS_DIR=""
  OUTPUT_DIR="${1%/}"
  [[ ! -d "$EXISTING_ASSEMBLIES_DIR" ]] && echo "Error: assemblies dir not found: $EXISTING_ASSEMBLIES_DIR" && exit 1
else
  # Standard mode: inputs.txt, ont_reads_dir, output_dir
  if [[ $# -lt 3 ]]; then usage; fi
  INPUTS_FILE="$1"
  ONT_READS_DIR="$2"
  OUTPUT_DIR="${3%/}"
  [[ ! -f "$INPUTS_FILE" ]]   && echo "Error: inputs file not found: $INPUTS_FILE"   && exit 1
  [[ ! -d "$ONT_READS_DIR" ]] && echo "Error: ONT reads dir not found: $ONT_READS_DIR" && exit 1
fi

# --skip-bakta implies --skip-antismash (antiSMASH requires Bakta GBFF output)
if $SKIP_BAKTA && ! $SKIP_ANTISMASH; then
  SKIP_ANTISMASH=true
  echo "Note: --skip-bakta also skips antiSMASH (requires Bakta GBFF output)"
fi

# ------------------------------------------------------------------------------
# Config — load from .conf file if present, fall back to hardcoded defaults
# ------------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF_FILE="$SCRIPT_DIR/dragonsmash.conf"

# Unset any inherited environment variables so they can't override script defaults
unset BAKTA_DB
unset BAKTA_AMRFINDER_DB
unset ANTISMASH_JOBS
unset ANTISMASH_CPUS
unset CONDA_BASE

if [[ -f "$CONF_FILE" ]]; then
  # shellcheck source=/dev/null
  source "$CONF_FILE"
  tprint "Config loaded from: $CONF_FILE"
else
  tprint "Warning: No dragonsmash.conf found at $CONF_FILE"
  tprint "         Copy dragonsmash.conf.example to dragonsmash.conf and set your paths."
fi

# CONDA_BASE: auto-detect from conda if not set in conf
if [[ -z "${CONDA_BASE:-}" ]]; then
  CONDA_BASE=$(conda info --base 2>/dev/null || echo "")
  if [[ -z "$CONDA_BASE" ]]; then
    eprint "Error: Cannot detect conda base. Set CONDA_BASE in dragonsmash.conf."
    exit 1
  fi
fi

# Conda environment names — override in dragonsmash.conf if yours differ
ENV_NAME_DRAGONFLYE="${ENV_NAME_DRAGONFLYE:-dragonflye}"
ENV_NAME_SEQKIT="${ENV_NAME_SEQKIT:-seqkit}"
ENV_NAME_BAKTA="${ENV_NAME_BAKTA:-bakta}"
ENV_NAME_ANTISMASH="${ENV_NAME_ANTISMASH:-antismash8}"
ENV_NAME_GTDBTK="${ENV_NAME_GTDBTK:-gtdbtk}"
# barrnap and kraken2 — install into 16S_ID_env:
#   conda create -n 16S_ID_env -c bioconda -c conda-forge barrnap kraken2 -y
ENV_NAME_BARRNAP="${ENV_NAME_BARRNAP:-16S_ID_env}"

# 16S Kraken2 database — SILVA by default (present on MDU servers)
SILVA_DB="${SILVA_DB:-/home/shared/db/kraken2/silva}"

# GTDBTk skani sketch directory — used with GTDBTk >= 2.3 which switched from
# fastANI to skani. Leave unset to let GTDBTk use its default. Set explicitly
# in dragonsmash.conf when using a shared database like MDU release226.
# GTDBTK_SKANI_SKETCH_DIR is unset by default — override in conf.
: "${GTDBTK_SKANI_SKETCH_DIR:=}"

# antiSMASH taxon — 'bacteria' (default) or 'fungi'
ANTISMASH_TAXON="${ANTISMASH_TAXON:-bacteria}"

# antiSMASH extra flags — appended verbatim to the antismash command.
# Default enables MIBiG comparison, ClusterBlast, and RRE detection.
# Override in dragonsmash.conf to add/remove analyses e.g.:
#   ANTISMASH_EXTRA_ARGS="--cc-mibig --cb-general --cb-knownclusters"  # lighter run
#   ANTISMASH_EXTRA_ARGS="--cc-mibig --cb-general --cb-subcluster --cb-knownclusters --rre --pfam2go --smcog-trees"  # heavier
ANTISMASH_EXTRA_ARGS="${ANTISMASH_EXTRA_ARGS:---cc-mibig --cb-general --cb-subcluster --cb-knownclusters --rre}"

# barrnap kingdom — 'bac' (default), 'arc' (archaea), or 'euk' (eukaryote)
BARRNAP_KINGDOM="${BARRNAP_KINGDOM:-bac}"

# barrnap minimum score threshold (0–1, default 0.8)
BARRNAP_MIN_SCORE="${BARRNAP_MIN_SCORE:-0.8}"

# BAKTA_DB: required — no sensible default across servers
if [[ -z "${BAKTA_DB:-}" ]]; then
  eprint "Error: BAKTA_DB is not set. Add 'BAKTA_DB=/path/to/bakta/db' to dragonsmash.conf."
  eprint "       Download with: bakta_db download --output /your/db/path"
  exit 1
fi

# antiSMASH: run multiple samples in parallel, but cap total CPU use at $CPUS.
# ANTISMASH_JOBS controls how many samples run concurrently (default: 4, or 1 if
# CPUS < 4). ANTISMASH_CPUS is the per-job thread count = CPUS / ANTISMASH_JOBS.
# Override both in dragonsmash.conf if needed.
if [[ -z "${ANTISMASH_JOBS:-}" ]]; then
  ANTISMASH_JOBS=$(( CPUS >= 4 ? 4 : 1 ))
fi
ANTISMASH_CPUS=$(( CPUS / ANTISMASH_JOBS ))
(( ANTISMASH_CPUS < 1 )) && ANTISMASH_CPUS=1

# Bakta: run multiple samples in parallel (default: 2 jobs — diamond is memory hungry).
# BAKTA_JOBS controls how many samples run concurrently.
# BAKTA_CPUS is the per-job thread count = CPUS / BAKTA_JOBS.
# Override both in dragonsmash.conf if needed.
if [[ -z "${BAKTA_JOBS:-}" ]]; then
  BAKTA_JOBS=$(( CPUS >= 2 ? 2 : 1 ))
fi
BAKTA_CPUS=$(( CPUS / BAKTA_JOBS ))
(( BAKTA_CPUS < 1 )) && BAKTA_CPUS=1

tprint "CPUs: $CPUS (Bakta: ${BAKTA_JOBS} jobs × ${BAKTA_CPUS} CPUs each | antiSMASH: ${ANTISMASH_JOBS} jobs × ${ANTISMASH_CPUS} CPUs each | GTDBTk: $CPUS CPUs)"
$DRY_RUN && tprint "*** DRY RUN MODE — no commands will be executed ***"

# ------------------------------------------------------------------------------
# Lock file — prevent concurrent runs on the same output directory
# ------------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR"
LOCK_FILE="$OUTPUT_DIR/.dragonsmash.lock"

if [[ -f "$LOCK_FILE" ]]; then
  LOCK_PID=$(cat "$LOCK_FILE" 2>/dev/null || echo "unknown")
  if [[ "$LOCK_PID" != "unknown" ]] && kill -0 "$LOCK_PID" 2>/dev/null; then
    tprint "Error: Another dragonSMASH run appears to be active on this output directory."
    tprint "  Lock file: $LOCK_FILE"
    tprint "  PID in lock: $LOCK_PID"
    tprint "  If no run is active, delete the lock file and retry:"
    tprint "    rm $LOCK_FILE"
    exit 1
  else
    tprint "Warning: Stale lock file found (PID $LOCK_PID is not running) — removing and continuing."
    rm -f "$LOCK_FILE"
  fi
fi

echo $$ > "$LOCK_FILE"
# Set trap immediately after writing lock — cleans up even if log setup hasn't run yet.
# The 'exec 1>&3' is guarded: fd 3 only exists after log redirect, so we check first.
trap 'rm -f "$LOCK_FILE"; { [[ -n "${WATCHDOG_PID:-}" ]] && kill "$WATCHDOG_PID" 2>/dev/null || true; }; { [[ -e /proc/self/fd/3 ]] && exec 1>&3 2>&1; } 2>/dev/null || true' EXIT INT TERM

# ------------------------------------------------------------------------------
# Log file setup — redirect all verbose output to log, terminal gets tprint only
# ------------------------------------------------------------------------------
LOG_FILE="$OUTPUT_DIR/dragonsmash_$(date +%d%m%Y_%H%M%S).log"
exec 3>&1                    # save terminal as fd 3
exec 1>>"$LOG_FILE" 2>&1    # redirect stdout+stderr to log

tprint() { printf '%s\n' "$*" >&3; }                              # terminal only
eprint() { printf '%s\n' "$*" >&3; printf '%s\n' "$*" >&2; }     # terminal + log
bprint() { printf '%s\n' "$*" >&3; printf '%s\n' "$*"; }         # terminal + log

# progress_bar CURRENT TOTAL LABEL
# Prints an in-place updating progress bar to the terminal.
# Call with CURRENT=0 before the loop to draw the initial bar,
# and CURRENT=TOTAL after the loop to finalise it with a newline.
_PB_WIDTH=40
progress_bar() {
  local current="$1" total="$2" label="$3"
  if (( total == 0 )); then return; fi
  local filled=$(( current * _PB_WIDTH / total ))
  local empty=$(( _PB_WIDTH - filled ))
  local bar
  bar="$(printf '%*s' "$filled" '' | tr ' ' '#')$(printf '%*s' "$empty" '' | tr ' ' '-')"
  if (( current == total )); then
    printf '\r  [%s] %d/%d %s\n' "$bar" "$current" "$total" "$label" >&3
  else
    printf '\r  [%s] %d/%d %s' "$bar" "$current" "$total" "$label" >&3
  fi
}

tprint ""
tprint "  Log file: $LOG_FILE"
tprint "  Monitor: tail -f $LOG_FILE"
tprint ""
bprint "                \\||/"
bprint "                |  @___oo"
bprint "      /\  /\   / (__,,,,|"
bprint "     ) /^\) ^\/  _)"
bprint "     )   /^\/   _)"
bprint "     )   _ /  / _)"
bprint " /\  )/\/ ||  | )_)"
bprint "<  >      |(,,) )__)"
bprint " ||      /    \)___)\\"
bprint " | \____(      )___) )___"
bprint "  \______(_______;;; __;;;"
bprint "   ---- dragonSMASH ----"
bprint "   Dragonflye → SeqKit → Bakta → antiSMASH → GTDBTk"
bprint ""

# ------------------------------------------------------------------------------
# Disk space check
# ------------------------------------------------------------------------------
check_disk_space() {
  local dir="$1"
  local min_gb="$2"
  local free_kb
  free_kb=$(df -k "$dir" | awk 'NR==2 {print $4}')
  local free_gb=$(( free_kb / 1024 / 1024 ))
  if (( free_gb < min_gb )); then
    tprint "Warning: Only ${free_gb} GB free on the filesystem containing $dir."
    tprint "         Recommended minimum is ${min_gb} GB."
    tprint "         Continue anyway? (Ctrl-C to abort, Enter to proceed)"
    read -r < /dev/tty
  else
    tprint "Disk space OK: ${free_gb} GB free (minimum: ${min_gb} GB)"
  fi
}

if ! $DRY_RUN; then
  check_disk_space "$OUTPUT_DIR" "$MIN_DISK_GB"
fi

# ------------------------------------------------------------------------------
# Bakta annotation parameters — genus and gram stain
# ------------------------------------------------------------------------------
tprint ""
if [[ -t 0 ]] || [[ -e /dev/tty ]]; then
  # Interactive — prompt the user
  printf 'Genus (e.g. Mycobacterium — leave blank to skip, auto-skips in 2 min): ' >&3
  read -r -t 120 BAKTA_GENUS < /dev/tty || { tprint "No input — genus left blank."; BAKTA_GENUS=""; }
  printf 'Gram stain — enter +, - or ? (unknown, auto-selects ? in 2 min): ' >&3
  read -r -t 120 BAKTA_GRAM < /dev/tty || { tprint "No input — gram defaulting to ?."; BAKTA_GRAM="?"; }
else
  # Non-interactive (piped/scripted) — use defaults silently
  BAKTA_GENUS=""
  BAKTA_GRAM="?"
  tprint "Non-interactive mode — genus: <not set>, gram: ?"
fi
# Validate gram input
if [[ "$BAKTA_GRAM" != "+" && "$BAKTA_GRAM" != "-" && "$BAKTA_GRAM" != "?" ]]; then
  tprint "Invalid gram stain value '$BAKTA_GRAM' — defaulting to ? (unknown)"
  BAKTA_GRAM="?"
fi
tprint "Genus: ${BAKTA_GENUS:-<not set>}  |  Gram: $BAKTA_GRAM"

# ------------------------------------------------------------------------------
# Configuration summary — printed to terminal and log
# ------------------------------------------------------------------------------
bprint ""
bprint "=============================="
bprint " Configuration"
bprint "=============================="
SCRIPT_VERSION=$(grep '^# Version' "${BASH_SOURCE[0]}" | tail -1 | awk '{print $NF}')
bprint "  Script version : $SCRIPT_VERSION"
bprint "  Output dir     : $OUTPUT_DIR"
bprint "  Bakta DB       : $BAKTA_DB"
bprint "  Bakta genus    : ${BAKTA_GENUS:-<not set>}"
bprint "  Bakta gram     : $BAKTA_GRAM"
bprint "  Max contigs    : $([ "$MAX_CONTIGS" -eq 0 ] && echo "disabled" || echo "$MAX_CONTIGS")"
bprint "  Bakta jobs     : $BAKTA_JOBS × $BAKTA_CPUS CPUs"
bprint "  antiSMASH jobs : $ANTISMASH_JOBS × $ANTISMASH_CPUS CPUs (total: $CPUS)"
bprint "  antiSMASH taxon: $ANTISMASH_TAXON"
bprint "  antiSMASH flags: $ANTISMASH_EXTRA_ARGS"
bprint "  barrnap kingdom: $BARRNAP_KINGDOM"
bprint "  barrnap min-scr: $BARRNAP_MIN_SCORE"
bprint "  GTDBTk env     : ${GTDBTK_CONDA_PREFIX:+$GTDBTK_CONDA_PREFIX (via GTDBTK_CONDA_PREFIX)}${GTDBTK_CONDA_PREFIX:-$CONDA_BASE/envs/${ENV_NAME_GTDBTK:-gtdbtk}}"
bprint "  GTDBTk DB      : ${GTDBTK_DATA_PATH:-<from conda env>}"
bprint "  GTDBTk skani   : ${GTDBTK_SKANI_SKETCH_DIR:-<default>}"
bprint "  Log file       : $LOG_FILE"
bprint "  Dry run        : $DRY_RUN"
bprint ""

# ------------------------------------------------------------------------------
# Input validation
# ------------------------------------------------------------------------------
tprint ""
tprint "=============================="
tprint " Validating inputs"
tprint "=============================="

find_ont_reads() {
  local sample="$1"
  local dir="$2"
  local path=""
  if [[ "$sample" == *.fastq.gz || "$sample" == *.fastq ]]; then
    path=$(find "$dir" -type f -name "$sample" -print -quit)
  else
    path=$(find "$dir" -type f -name "${sample}.fastq.gz" -print -quit)
    [[ -z "$path" ]] && path=$(find "$dir" -type f -name "${sample}.fastq" -print -quit)
  fi
  echo "$path"
}

find_existing_fasta() {
  local sample="$1"
  local dir="$2"
  local path=""
  for ext in fasta fa fna fas; do
    path=$(find "$dir" -maxdepth 1 -type f -name "${sample}.${ext}" -print -quit)
    [[ -n "$path" ]] && break
  done
  echo "$path"
}

MISSING_READS=()
ALL_SAMPLES=()

if [[ -n "$EXISTING_ASSEMBLIES_DIR" ]]; then
  # Infer sample names from FASTA filenames in the assemblies directory
  while IFS= read -r fasta; do
    fname=$(basename "$fasta")
    sample="${fname%.*}"
    ALL_SAMPLES+=("$sample")
  done < <(find "$EXISTING_ASSEMBLIES_DIR" -maxdepth 1 -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" -o -name "*.fas" \) | sort)

  if [[ ${#ALL_SAMPLES[@]} -eq 0 ]]; then
    eprint "Error: No FASTA files found in $EXISTING_ASSEMBLIES_DIR"
    exit 1
  fi
  tprint "Samples found from FASTA files: ${#ALL_SAMPLES[@]}"
  for s in "${ALL_SAMPLES[@]}"; do tprint "  $s"; done
else
  # Standard mode: read samples from inputs.txt and validate reads exist
  while IFS= read -r sample || [[ -n "${sample:-}" ]]; do
    [[ -z "${sample// }" || "$sample" =~ ^# ]] && continue
    sample="${sample%.fastq.gz}"; sample="${sample%.fastq}"
    sample="${sample%.fasta}";    sample="${sample%.fa}"; sample="${sample%.fna}"
    ALL_SAMPLES+=("$sample")

    reads=$(find_ont_reads "$sample" "$ONT_READS_DIR")
    if [[ -z "$reads" ]]; then
      MISSING_READS+=("$sample")
    fi
  done < "$INPUTS_FILE"

  if [[ ${#ALL_SAMPLES[@]} -eq 0 ]]; then
    eprint "Error: No samples found in $INPUTS_FILE"
    exit 1
  fi

  tprint "Samples found: ${#ALL_SAMPLES[@]}"

  if [[ ${#MISSING_READS[@]} -gt 0 ]]; then
    tprint ""
    eprint "Error: The following samples have no ONT reads in $ONT_READS_DIR:"
    for s in "${MISSING_READS[@]}"; do
      tprint "  - $s"
    done
    tprint ""
    tprint "Remove these from $INPUTS_FILE or add the missing read files, then re-run."
    exit 1
  fi

  tprint "All samples have reads — OK"
fi

# ------------------------------------------------------------------------------
# Check conda envs and tools
# ------------------------------------------------------------------------------
ENV_DRAGONFLYE="$CONDA_BASE/envs/$ENV_NAME_DRAGONFLYE"
ENV_SEQKIT="$CONDA_BASE/envs/$ENV_NAME_SEQKIT"
ENV_BAKTA="$CONDA_BASE/envs/$ENV_NAME_BAKTA"
ENV_ANTISMASH="$CONDA_BASE/envs/$ENV_NAME_ANTISMASH"
# Allow GTDBTK_CONDA_PREFIX in conf to point to an env outside CONDA_BASE
# (e.g. another user's environment). Falls back to standard path if not set.
ENV_GTDBTK="${GTDBTK_CONDA_PREFIX:-$CONDA_BASE/envs/$ENV_NAME_GTDBTK}"
# barrnap may share an env with gtdbtk or live elsewhere
ENV_BARRNAP="$CONDA_BASE/envs/$ENV_NAME_BARRNAP"

ENVS_TO_CHECK=("$ENV_SEQKIT" "$ENV_BAKTA" "$ENV_ANTISMASH" "$ENV_BARRNAP" "$ENV_GTDBTK")
[[ -z "$EXISTING_ASSEMBLIES_DIR" ]] && ENVS_TO_CHECK=("$ENV_DRAGONFLYE" "${ENVS_TO_CHECK[@]}")

for env in "${ENVS_TO_CHECK[@]}"; do
  if [[ ! -d "$env" ]]; then
    eprint "Error: conda environment not found: $env"
    exit 1
  fi
done

# SILVA Kraken2 database check
if [[ ! -f "$SILVA_DB/hash.k2d" ]]; then
  eprint "Error: SILVA Kraken2 database not found at $SILVA_DB"
  eprint "       Set SILVA_DB in dragonsmash.conf to the correct path."
  exit 1
fi

run_in_env() {
  local env_path="$1"; shift
  local cmd="$1";      shift
  env PATH="$env_path/bin:$PATH" "$env_path/bin/$cmd" "$@"
}

# ------------------------------------------------------------------------------
# Diamond watchdog — background process that monitors diamond jobs during
# Stages 3 and 4. Detects hung diamond processes using CPU jiffies from
# /proc/PID/stat rather than output file size (diamond batch-writes its output
# at the very end, so the output file stays at 0 bytes for the entire run).
#
# A diamond job is considered hung if its CPU jiffies have not increased for
# WATCHDOG_STALL_CHECKS consecutive check intervals. On detection, diamond and
# its parent bakta/antismash process are killed, and the partial sample output
# directory is removed so parallel can restart the sample cleanly.
#
# Usage:
#   start_watchdog <stage_output_dir>  # arg accepted but unused; kept for clarity
#   stop_watchdog
# The watchdog PID is stored in WATCHDOG_PID.
# Configure WATCHDOG_INTERVAL_SECS and WATCHDOG_STALL_CHECKS in dragonsmash.conf.
# ------------------------------------------------------------------------------
# Watchdog configuration.
# Detection logic: read CPU jiffies from /proc/PID/stat at each interval.
# If CPU jiffies have not increased for WATCHDOG_STALL_CHECKS consecutive
# checks, the process is considered hung and killed.
# This correctly handles diamond's batch-write pattern where output files
# stay at 0 bytes for the entire run — file size is NOT used.
#
# Defaults: check every 10 min, kill after 3 stalled checks = 30 min no CPU.
# Override WATCHDOG_INTERVAL_SECS and WATCHDOG_STALL_CHECKS in dragonsmash.conf.
WATCHDOG_INTERVAL_SECS=${WATCHDOG_INTERVAL_SECS:-600}   # check every 10 minutes
WATCHDOG_STALL_CHECKS=${WATCHDOG_STALL_CHECKS:-3}        # kill after N consecutive stalled checks
WATCHDOG_PID=""
WATCHDOG_STATE_DIR=""

start_watchdog() {
  local stage_outdir="$1"
  local watchdog_log="$LOG_FILE"   # captured at fork time

  WATCHDOG_STATE_DIR=$(mktemp -d /tmp/dragonsmash_watchdog_XXXXXX)

  (
    local state_dir="$WATCHDOG_STATE_DIR"
    local interval=$WATCHDOG_INTERVAL_SECS
    local stall_limit=$WATCHDOG_STALL_CHECKS
    local log_prefix="[watchdog]"

    wlog() {
      printf '%s\n' "$*" >&3 2>/dev/null || true
      [[ -n "$watchdog_log" ]] && printf '%s\n' "$*" >> "$watchdog_log" 2>/dev/null || true
    }

    # get_cpu_jiffies PID
    # Returns total CPU jiffies (utime+stime) from /proc/PID/stat, or -1 on error.
    get_cpu_jiffies() {
      local pid="$1"
      local stat
      stat=$(cat /proc/"$pid"/stat 2>/dev/null) || { echo -1; return; }
      # Fields 14 (utime) and 15 (stime) in /proc/PID/stat
      echo "$stat" | awk '{print $14 + $15}'
    }

    while true; do
      sleep "$interval"

      while IFS= read -r pid; do
        [[ -z "$pid" ]] && continue

        local current_jiffies
        current_jiffies=$(get_cpu_jiffies "$pid")
        [[ "$current_jiffies" == "-1" ]] && continue   # process gone

        # State file: <last_jiffies>:<stall_count>
        local state_file="$state_dir/$pid"

        if [[ -f "$state_file" ]]; then
          local last_jiffies stall_count
          last_jiffies=$(cut -d: -f1 "$state_file")
          stall_count=$(cut -d: -f2  "$state_file")

          if [[ "$current_jiffies" -gt "$last_jiffies" ]]; then
            # CPU is progressing — reset stall count
            echo "${current_jiffies}:0" > "$state_file"
          else
            # No CPU progress this interval
            stall_count=$(( stall_count + 1 ))
            wlog "$log_prefix diamond PID $pid: no CPU progress (stall ${stall_count}/${stall_limit}, jiffies=${current_jiffies})"

            if [[ "$stall_count" -ge "$stall_limit" ]]; then
              wlog "$log_prefix diamond PID $pid confirmed hung — killing"

              local parent_pid
              parent_pid=$(awk '/PPid/{print $2}' /proc/"$pid"/status 2>/dev/null || echo "")

              # Read parent cmdline BEFORE killing — it will be gone afterward
              local parent_cmdline=""
              if [[ -n "$parent_pid" ]]; then
                parent_cmdline=$(tr '\0' ' ' < /proc/"$parent_pid"/cmdline 2>/dev/null) || true
              fi

              wlog "$log_prefix killing diamond PID $pid"
              kill "$pid" 2>/dev/null || true

              if [[ -n "$parent_pid" ]] && \
                 echo "$parent_cmdline" | grep -q "bakta\|antismash"; then
                wlog "$log_prefix killing parent PID $parent_pid"
                kill "$parent_pid" 2>/dev/null || true
              fi

              # Find sample output dir from parent cmdline (--output for bakta,
              # --output-dir for antismash) and remove it so parallel can restart.
              local sample_outdir=""
              if [[ -n "$parent_cmdline" ]]; then
                sample_outdir=$(echo "$parent_cmdline" \
                  | grep -oP '(?<=--output(?:-dir)? )\S+' | head -1) || true
              fi

              if [[ -n "$sample_outdir" ]] && [[ -d "$sample_outdir" ]]; then
                wlog "$log_prefix removing partial output: $sample_outdir"
                rm -rf "$sample_outdir"
              else
                wlog "$log_prefix could not identify sample output dir — skipping cleanup"
              fi

              rm -f "$state_file"
            else
              echo "${current_jiffies}:${stall_count}" > "$state_file"
            fi
          fi
        else
          # First time seeing this PID
          echo "${current_jiffies}:0" > "$state_file"
          wlog "$log_prefix now watching diamond PID $pid (jiffies=${current_jiffies})"
        fi

      done < <(pgrep -u "$(whoami)" diamond 2>/dev/null)

    done
  ) &

  WATCHDOG_PID=$!
  tprint "Watchdog started (PID $WATCHDOG_PID, interval: $(( WATCHDOG_INTERVAL_SECS / 60 ))min, kill after ${WATCHDOG_STALL_CHECKS} stalled checks = $(( WATCHDOG_INTERVAL_SECS * WATCHDOG_STALL_CHECKS / 60 ))min no CPU)"
}

stop_watchdog() {
  if [[ -n "$WATCHDOG_PID" ]]; then
    kill "$WATCHDOG_PID" 2>/dev/null || true
    wait "$WATCHDOG_PID" 2>/dev/null || true
    WATCHDOG_PID=""
    tprint "Watchdog stopped."
  fi
  if [[ -n "$WATCHDOG_STATE_DIR" ]] && [[ -d "$WATCHDOG_STATE_DIR" ]]; then
    rm -rf "$WATCHDOG_STATE_DIR"
    WATCHDOG_STATE_DIR=""
  fi
}

[[ ! -d "$BAKTA_DB" ]] && eprint "Error: Bakta DB not found at: $BAKTA_DB" && exit 1

# Check Bakta DB has expected number of files
BAKTA_DB_FILE_COUNT=$(ls "$BAKTA_DB" | wc -l)
if [[ "$BAKTA_DB_FILE_COUNT" -lt 25 ]]; then
  eprint "Error: Bakta DB at $BAKTA_DB looks incomplete ($BAKTA_DB_FILE_COUNT files found, expected ~29)."
  eprint "       Re-run bakta_db download or check the path in dragonsmash.conf."
  exit 1
fi
tprint "Bakta DB OK: $BAKTA_DB ($BAKTA_DB_FILE_COUNT files)"

# Check AMRFinderPlus DB and warn if outdated
AMRFINDER_DB="$BAKTA_DB/amrfinderplus-db"
if [[ -d "$AMRFINDER_DB" ]]; then
  AMRFINDER_LATEST=$(readlink "$AMRFINDER_DB/latest" 2>/dev/null || echo "unknown")
  # Warn if latest DB version is older than 6 months
  if [[ "$AMRFINDER_LATEST" != "unknown" ]]; then
    DB_DATE=$(echo "$AMRFINDER_LATEST" | grep -oP '^\d{4}-\d{2}-\d{2}' || echo "")
    if [[ -n "$DB_DATE" ]]; then
      DB_EPOCH=$(date -d "$DB_DATE" +%s 2>/dev/null || echo 0)
      NOW_EPOCH=$(date +%s)
      AGE_DAYS=$(( (NOW_EPOCH - DB_EPOCH) / 86400 ))
      if [[ "$AGE_DAYS" -gt 180 ]]; then
        tprint "Warning: AMRFinderPlus DB is ${AGE_DAYS} days old (version: $AMRFINDER_LATEST)."
        tprint "         Consider running: amrfinder_update --force_update --database $AMRFINDER_DB"
      else
        tprint "AMRFinderPlus DB OK: $AMRFINDER_LATEST (${AGE_DAYS} days old)"
      fi
    fi
  fi
else
  tprint "Warning: AMRFinderPlus DB not found at $AMRFINDER_DB — AMR annotation may fail."
fi

if ! command -v parallel >/dev/null 2>&1; then
  tprint "Error: 'parallel' not found on PATH."
  exit 1
fi

# ------------------------------------------------------------------------------
# Stage output dirs
# ------------------------------------------------------------------------------
ASSEMBLY_DIR="$OUTPUT_DIR/01_assembly"
SEQKIT_DIR="$OUTPUT_DIR/02_seqkit_stats"
SEQKIT_ASSEMBLY_TSV="$SEQKIT_DIR/assembly_stats.tsv"
SEQKIT_READS_TSV="$SEQKIT_DIR/read_stats.tsv"
BAKTA_DIR="$OUTPUT_DIR/03_bakta"
ANTISMASH_DIR="$OUTPUT_DIR/04_antismash"
RRNA_DIR="$OUTPUT_DIR/05_16S"
GTDBTK_DIR="$OUTPUT_DIR/06_gtdbtk"
SUMMARY_DIR="$OUTPUT_DIR/07_summary"
SUMMARY_TSV="$SUMMARY_DIR/summary.tsv"

if ! $DRY_RUN; then
  mkdir -p "$ASSEMBLY_DIR" "$SEQKIT_DIR" "$BAKTA_DIR" "$ANTISMASH_DIR" "$RRNA_DIR" "$GTDBTK_DIR" "$SUMMARY_DIR"
fi

# Illumina check
ILLUMINA_READS_PROVIDED=false
if [[ -n "$ILLUMINA_READS_DIR" && -d "$ILLUMINA_READS_DIR" ]]; then
  ILLUMINA_READS_PROVIDED=true
else
  tprint "No Illumina reads directory provided — ONT-only assembly."
fi

# Track failed samples
FAILED_SAMPLES=()

# Pipeline start time
PIPELINE_START=$(date +%s)

# ==============================================================================
# STAGE 1: Dragonflye assembly (or symlink existing assemblies)
# ==============================================================================
tprint ""
tprint "=============================="
tprint " STAGE 1: Dragonflye Assembly"
tprint "=============================="

if [[ -n "$EXISTING_ASSEMBLIES_DIR" ]]; then
  tprint "Skipping assembly — symlinking existing FASTA files from: $EXISTING_ASSEMBLIES_DIR"
  local_total=${#ALL_SAMPLES[@]}; local_done=0
  progress_bar 0 "$local_total" "Stage 1 (linking)"
  for sample in "${ALL_SAMPLES[@]}"; do
    outdir="$ASSEMBLY_DIR/$sample"
    link="$outdir/final.contigs.fa"
    if [[ -L "$link" || -f "$link" ]]; then
      (( local_done++ ))
      progress_bar "$local_done" "$local_total" "Stage 1 (linking)"
      continue
    fi
    src=$(find_existing_fasta "$sample" "$EXISTING_ASSEMBLIES_DIR")
    if [[ -z "$src" ]]; then
      eprint "ERROR [Stage 1]: No FASTA found for '$sample' in $EXISTING_ASSEMBLIES_DIR"
      FAILED_SAMPLES+=("$sample (missing fasta)")
      (( local_done++ ))
      progress_bar "$local_done" "$local_total" "Stage 1 (linking)"
      continue
    fi
    if $DRY_RUN; then
      tprint "  [dry-run] ln -s $src $link"
    else
      mkdir -p "$outdir"
      ln -s "$(realpath "$src")" "$link"
    fi
    (( local_done++ ))
    progress_bar "$local_done" "$local_total" "Stage 1 (linking)"
  done
else
  local_total=${#ALL_SAMPLES[@]}; local_done=0
  progress_bar 0 "$local_total" "Stage 1 (assembly)"
  for sample in "${ALL_SAMPLES[@]}"; do

    outdir="$ASSEMBLY_DIR/$sample"
    assembly_fasta="$outdir/final.contigs.fa"

    if [[ -f "$assembly_fasta" ]]; then
      (( local_done++ ))
      progress_bar "$local_done" "$local_total" "Stage 1 (assembly)"
      continue
    fi

    reads_ont=$(find_ont_reads "$sample" "$ONT_READS_DIR")

    assembly_failed=false

    if $ILLUMINA_READS_PROVIDED; then
      shopt -s nullglob
      illumina_r1_files=("$ILLUMINA_READS_DIR/$sample/"*R1*.fastq*)
      shopt -u nullglob
      if [[ ${#illumina_r1_files[@]} -gt 0 && -f "${illumina_r1_files[0]}" ]]; then
        tprint "Assembling $sample (ONT + Illumina polishing)..."
        if $DRY_RUN; then
          echo "  [dry-run] dragonflye --cpus $CPUS --assembler flye --racon 1 --polypolish 1 --reads $reads_ont --outdir $outdir"
        else
          mkdir -p "$outdir"
          for r1 in "${illumina_r1_files[@]}"; do
            r2="${r1/_R1_/_R2_}"
          run_in_env "$ENV_DRAGONFLYE" dragonflye --cpus "$CPUS" --assembler flye --racon 1 --polypolish 1 --force \
              --reads "$reads_ont" --outdir "$outdir" \
              --R1 "$r1" --R2 "$r2" || { assembly_failed=true; break; }
          done
        fi
      else
        tprint "No Illumina reads found for $sample — assembling ONT only..."
        if $DRY_RUN; then
          echo "  [dry-run] dragonflye --cpus $CPUS --assembler flye --racon 1 --reads $reads_ont --outdir $outdir"
        else
          mkdir -p "$outdir"
          run_in_env "$ENV_DRAGONFLYE" dragonflye --cpus "$CPUS" --assembler flye --racon 1 --force \
            --reads "$reads_ont" --outdir "$outdir" || assembly_failed=true
        fi
      fi
    else
      tprint "Assembling $sample (ONT only)..."
      if $DRY_RUN; then
        echo "  [dry-run] dragonflye --cpus $CPUS --assembler flye --racon 1 --reads $reads_ont --outdir $outdir"
      else
        mkdir -p "$outdir"
        run_in_env "$ENV_DRAGONFLYE" dragonflye --cpus "$CPUS" --assembler flye --racon 1 --force \
          --reads "$reads_ont" --outdir "$outdir" || assembly_failed=true
      fi
    fi

    if $assembly_failed; then
      eprint "ERROR [Stage 1]: Assembly failed for '$sample' — check $outdir/dragonflye.log for details."
      FAILED_SAMPLES+=("$sample (assembly)")
      (( local_done++ ))
      progress_bar "$local_done" "$local_total" "Stage 1 (assembly)"
      continue
    fi

    # Ensure final.contigs.fa exists for downstream stages.
    # Illumina-polished runs produce final.contigs.fa directly; ONT-only runs
    # produce contigs.fa — symlink it so Stages 2–4 have a consistent target.
    if $DRY_RUN; then
      tprint "  [dry-run] ln -sf contigs.fa $outdir/final.contigs.fa (if not already present)"
    elif [[ ! -f "$outdir/final.contigs.fa" ]]; then
      ln -sf contigs.fa "$outdir/final.contigs.fa"
    fi

    $DRY_RUN || echo "Assembly done: $sample"
    (( local_done++ ))
    progress_bar "$local_done" "$local_total" "Stage 1 (assembly)"

  done
fi

# ==============================================================================
# STAGE 2: SeqKit Stats
# Two outputs:
#   assembly_stats.tsv — per-sample assembly metrics (num_contigs, N50, N90, L50, etc.)
#   read_stats.tsv     — per-sample raw read metrics (skipped in --assemblies mode)
# ==============================================================================
tprint ""
tprint "=============================="
tprint " STAGE 2: SeqKit Stats"
tprint "=============================="

# TSV headers
ASSEMBLY_TSV_HEADER=$'sample\tnum_contigs\tassembly_size_bp\tmin_len\tavg_len\tmax_len\tN50\tN90\tL50'
READS_TSV_HEADER=$'sample\tread_type\tnum_reads\ttotal_bases\tmin_len\tavg_len\tmax_len\tN50'

if $DRY_RUN; then
  tprint "  [dry-run] seqkit stats -a -N 90 on each assembly → $SEQKIT_ASSEMBLY_TSV"
  [[ -z "$EXISTING_ASSEMBLIES_DIR" ]] && \
    tprint "  [dry-run] seqkit stats -a on each read file → $SEQKIT_READS_TSV"
else
  # --------------------------------------------------------------------------
  # Assembly stats
  # --------------------------------------------------------------------------
  declare -A SEQKIT_ASM_DONE
  if [[ -f "$SEQKIT_ASSEMBLY_TSV" ]]; then
    while IFS=$'\t' read -r s _rest; do
      SEQKIT_ASM_DONE["$s"]=1
    done < <(tail -n +2 "$SEQKIT_ASSEMBLY_TSV")
  else
    printf '%s\n' "$ASSEMBLY_TSV_HEADER" > "$SEQKIT_ASSEMBLY_TSV"
  fi

  local_total=${#ALL_SAMPLES[@]}; local_done=0
  progress_bar 0 "$local_total" "Stage 2 (assembly stats)"
  for sample in "${ALL_SAMPLES[@]}"; do
    assembly_fasta="$ASSEMBLY_DIR/$sample/final.contigs.fa"

    if [[ ! -f "$assembly_fasta" ]]; then
      tprint "Warning: No assembly found for '$sample' — skipping assembly SeqKit stats."
      (( local_done++ )); progress_bar "$local_done" "$local_total" "Stage 2 (assembly stats)"
      continue
    fi

    if [[ -n "${SEQKIT_ASM_DONE[$sample]+_}" ]]; then
      (( local_done++ )); progress_bar "$local_done" "$local_total" "Stage 2 (assembly stats)"
      continue
    fi

    # -a: all stats  -N 90: also compute N90  -T: tabular (tab-separated, no padding)
    # seqkit stats -T columns (1-indexed):
    #   1:file 2:format 3:type 4:num_seqs 5:sum_len 6:min_len 7:avg_len 8:max_len
    #   9:Q20(%) 10:Q30(%) 11:GC(%) 12:sum_gap 13:N50 14:Q20_bases 15:Q30_bases
    #   16:AvgQual 17:N90 (from -N 90)  18:L50 (position in seqkit -a output)
    # With -a -N 90 -T the column order is stable; we extract by header name for safety.
    run_in_env "$ENV_SEQKIT" seqkit stats -a -N 90 -T "$assembly_fasta" \
      2>/dev/null | awk -v sample="$sample" '
      NR==1 {
        for(i=1;i<=NF;i++) col[$i]=i
        next
      }
      {
        num_seqs = $col["num_seqs"]
        sum_len  = $col["sum_len"]
        min_len  = $col["min_len"]
        avg_len  = $col["avg_len"]
        max_len  = $col["max_len"]
        n50      = (("N50"  in col) ? $col["N50"]  : "NA")
        n90      = (("N90"  in col) ? $col["N90"]  : "NA")
        l50      = (("L50"  in col) ? $col["L50"]  : "NA")
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
          sample, num_seqs, sum_len, min_len, avg_len, max_len, n50, n90, l50
      }' >> "$SEQKIT_ASSEMBLY_TSV" \
      || tprint "WARNING [Stage 2]: seqkit assembly stats failed for '$sample'."

    (( local_done++ ))
    progress_bar "$local_done" "$local_total" "Stage 2 (assembly stats)"
  done

  # --------------------------------------------------------------------------
  # Read stats (standard mode only — skipped with --assemblies)
  # --------------------------------------------------------------------------
  if [[ -z "$EXISTING_ASSEMBLIES_DIR" ]]; then
    declare -A SEQKIT_READ_DONE
    if [[ -f "$SEQKIT_READS_TSV" ]]; then
      while IFS=$'\t' read -r s _rest; do
        SEQKIT_READ_DONE["$s"]=1
      done < <(tail -n +2 "$SEQKIT_READS_TSV")
    else
      printf '%s\n' "$READS_TSV_HEADER" > "$SEQKIT_READS_TSV"
    fi

    local_total=${#ALL_SAMPLES[@]}; local_done=0
    progress_bar 0 "$local_total" "Stage 2 (read stats)"
    for sample in "${ALL_SAMPLES[@]}"; do
      if [[ -n "${SEQKIT_READ_DONE[$sample]+_}" ]]; then
        (( local_done++ )); progress_bar "$local_done" "$local_total" "Stage 2 (read stats)"
        continue
      fi

      # ONT reads
      reads_ont=$(find_ont_reads "$sample" "$ONT_READS_DIR")
      if [[ -n "$reads_ont" ]]; then
        run_in_env "$ENV_SEQKIT" seqkit stats -a -T "$reads_ont" \
          2>/dev/null | awk -v sample="$sample" -v rtype="ONT" '
          NR==1 { for(i=1;i<=NF;i++) col[$i]=i; next }
          {
            printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
              sample, rtype,
              $col["num_seqs"], $col["sum_len"],
              $col["min_len"],  $col["avg_len"],
              $col["max_len"],  (("N50" in col) ? $col["N50"] : "NA")
          }' >> "$SEQKIT_READS_TSV" \
          || tprint "WARNING [Stage 2]: seqkit read stats failed for '$sample' (ONT)."
      fi

      # Illumina reads (if provided)
      if $ILLUMINA_READS_PROVIDED; then
        illumina_r1=($(find "$ILLUMINA_READS_DIR/$sample/" -name "*R1*.fastq*" 2>/dev/null))
        illumina_r2=($(find "$ILLUMINA_READS_DIR/$sample/" -name "*R2*.fastq*" 2>/dev/null))
        for rfile in "${illumina_r1[@]}" "${illumina_r2[@]}"; do
          [[ ! -f "$rfile" ]] && continue
          local_rtype="Illumina_$(basename "$rfile" | grep -oP 'R[12]')"
          run_in_env "$ENV_SEQKIT" seqkit stats -a -T "$rfile" \
            2>/dev/null | awk -v sample="$sample" -v rtype="$local_rtype" '
            NR==1 { for(i=1;i<=NF;i++) col[$i]=i; next }
            {
              printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                sample, rtype,
                $col["num_seqs"], $col["sum_len"],
                $col["min_len"],  $col["avg_len"],
                $col["max_len"],  (("N50" in col) ? $col["N50"] : "NA")
            }' >> "$SEQKIT_READS_TSV" \
            || tprint "WARNING [Stage 2]: seqkit read stats failed for '$sample' ($local_rtype)."
        done
      fi

      (( local_done++ ))
      progress_bar "$local_done" "$local_total" "Stage 2 (read stats)"
    done
    tprint "Read stats written to:     $SEQKIT_READS_TSV"
  fi

  tprint "Assembly stats written to: $SEQKIT_ASSEMBLY_TSV"
fi

# ==============================================================================
# Contig filter — exclude highly fragmented assemblies from annotation
# ==============================================================================
FILTERED_SAMPLES=()
declare -A FILTERED_CONTIG_COUNTS
if [[ "$MAX_CONTIGS" -gt 0 ]]; then
  PASSING_SAMPLES=()
  for sample in "${ALL_SAMPLES[@]}"; do
    assembly_fasta="$ASSEMBLY_DIR/$sample/final.contigs.fa"
    if [[ -f "$assembly_fasta" ]]; then
      contig_count=$(grep -c "^>" "$assembly_fasta" 2>/dev/null || echo 0)
      if [[ "$contig_count" -gt "$MAX_CONTIGS" ]]; then
        tprint "SKIPPING $sample — $contig_count contigs exceeds --max-contigs $MAX_CONTIGS"
        FILTERED_SAMPLES+=("$sample")
        FILTERED_CONTIG_COUNTS["$sample"]="$contig_count"
      else
        PASSING_SAMPLES+=("$sample")
      fi
    else
      PASSING_SAMPLES+=("$sample")
    fi
  done
  ALL_SAMPLES=("${PASSING_SAMPLES[@]}")
  if [[ ${#FILTERED_SAMPLES[@]} -gt 0 ]]; then
    tprint "${#FILTERED_SAMPLES[@]} sample(s) excluded from annotation (use --max-contigs 0 to disable)."
  fi
fi
# STAGE 3: Bakta annotation
# ==============================================================================
tprint ""
tprint "=============================="
tprint " STAGE 3: Bakta Annotation"
tprint "=============================="

if $SKIP_BAKTA; then
  tprint "  [skipped] --skip-bakta flag set"
else

if ! $DRY_RUN && $WATCHDOG_ENABLED; then
  start_watchdog "$BAKTA_DIR"
fi

if $DRY_RUN; then
  for sample in "${ALL_SAMPLES[@]}"; do
    genus_str="${BAKTA_GENUS:+--genus $BAKTA_GENUS }"
    tprint "  [dry-run] bakta --db $BAKTA_DB --output $BAKTA_DIR/$sample --prefix $sample --threads $BAKTA_CPUS --compliant --force --skip-pseudo ${genus_str}--gram $BAKTA_GRAM $ASSEMBLY_DIR/$sample/final.contigs.fa"
  done
else
  run_bakta() {
    local sample="$1"
    local assembly_fasta="$ASSEMBLY_DIR/$sample/final.contigs.fa"
    local bakta_outdir="$BAKTA_DIR/$sample"

    if [[ ! -f "$assembly_fasta" ]]; then
      echo "Warning: No assembly for '$sample' — skipping Bakta." >&2
      return 0
    fi

    if [[ -f "$bakta_outdir/${sample}.gbff" ]]; then
      echo "Skip (already done): $sample"
      return 0
    fi

    echo "Annotating: $sample"

    local bakta_extra_args=()
    [[ -n "${BAKTA_GENUS:-}" ]] && bakta_extra_args+=(--genus "$BAKTA_GENUS")
    bakta_extra_args+=(--gram "$BAKTA_GRAM")

    if ! conda run -p "$ENV_BAKTA" --no-capture-output bakta \
      --db "$BAKTA_DB" \
      --output "$bakta_outdir" \
      --prefix "$sample" \
      --threads "$BAKTA_CPUS" \
      --compliant \
      --force \
      --skip-pseudo \
      "${bakta_extra_args[@]}" \
      "$assembly_fasta"; then
      echo "ERROR [Stage 3]: Bakta failed for '$sample'" >&2
      return 1
    fi

    echo "Bakta done: $sample"
  }

  export -f run_bakta
  export BAKTA_CPUS BAKTA_DB BAKTA_DIR ASSEMBLY_DIR ENV_BAKTA BAKTA_GENUS BAKTA_GRAM

  # Progress bar for Stage 3: counter-file pattern (same as Stage 4)
  local_total=${#ALL_SAMPLES[@]}
  _BK_COUNTER="$OUTPUT_DIR/.bakta_progress_$$"
  > "$_BK_COUNTER"
  export _BK_COUNTER local_total _PB_WIDTH

  run_bakta_tracked() {
    run_bakta "$@"
    local rc=$?
    echo x >> "$_BK_COUNTER"
    local done
    done=$(wc -l < "$_BK_COUNTER")
    local filled=$(( done * _PB_WIDTH / local_total ))
    local empty=$(( _PB_WIDTH - filled ))
    local bar
    bar="$(printf '%*s' "$filled" '' | tr ' ' '#')$(printf '%*s' "$empty" '' | tr ' ' '-')"
    if (( done >= local_total )); then
      printf '\r  [%s] %d/%d Stage 3 (bakta)\n' "$bar" "$done" "$local_total" > /dev/tty
    else
      printf '\r  [%s] %d/%d Stage 3 (bakta)' "$bar" "$done" "$local_total" > /dev/tty
    fi
    return $rc
  }
  export -f run_bakta_tracked

  progress_bar 0 "$local_total" "Stage 3 (bakta)"

  parallel -j "$BAKTA_JOBS" --linebuffer \
    run_bakta_tracked {} \
    ::: "${ALL_SAMPLES[@]}" || true   # failures detected by gbff scan below

  rm -f "$_BK_COUNTER"

  # Detect Bakta failures — parallel subshells can't write to FAILED_SAMPLES,
  # so scan for samples that have an assembly but no .gbff output.
  for sample in "${ALL_SAMPLES[@]}"; do
    if [[ -f "$ASSEMBLY_DIR/$sample/final.contigs.fa" && ! -f "$BAKTA_DIR/$sample/${sample}.gbff" ]]; then
      eprint "ERROR [Stage 3]: Bakta appears to have failed for '$sample' (no .gbff)"
      FAILED_SAMPLES+=("$sample (annotation)")
    fi
  done
fi  # end DRY_RUN check for stage 3

fi  # end skip-bakta check

# ==============================================================================
# STAGE 4: antiSMASH (using Bakta .gbff, keeping Bakta annotations)
# ==============================================================================
tprint ""
tprint "=============================="
tprint " STAGE 4: antiSMASH"
tprint "=============================="

if $SKIP_ANTISMASH; then
  tprint "  [skipped] --skip-antismash flag set"
else

if ! $DRY_RUN && $WATCHDOG_ENABLED; then
  stop_watchdog
  start_watchdog "$ANTISMASH_DIR"
fi

if $DRY_RUN; then
  for sample in "${ALL_SAMPLES[@]}"; do
    tprint "  [dry-run] antismash --genefinding-tool none --taxon $ANTISMASH_TAXON $ANTISMASH_EXTRA_ARGS --cpus $ANTISMASH_CPUS $BAKTA_DIR/$sample/${sample}.gbff → $ANTISMASH_DIR/$sample"
  done
else
  run_antismash() {
    local sample="$1"
    local bakta_outdir="$2"
    local antismash_outdir="$3"

    local gbff="$bakta_outdir/${sample}.gbff"

    if [[ ! -f "$gbff" ]]; then
      echo "Warning: No Bakta .gbff found for '$sample' — skipping antiSMASH." >&2
      return 0
    fi

    local out="$antismash_outdir/$sample"

    if [[ -d "$out" && -f "$out/index.html" ]]; then
      echo "Skip (already done): $sample"
      return 0
    fi

    mkdir -p "$out"

    echo "Running antiSMASH: $sample"
    if ! env PATH="$ENV_ANTISMASH/bin:$PATH" "$ENV_ANTISMASH/bin/antismash" \
      --genefinding-tool none \
      --taxon "$ANTISMASH_TAXON" \
      ${ANTISMASH_EXTRA_ARGS} \
      --cpus "$ANTISMASH_CPUS" \
      --output-dir "$out" \
      "$gbff"; then
      echo "ERROR [Stage 4]: antiSMASH failed on '$sample'" >&2
      return 1
    fi

    echo "antiSMASH done: $sample"
  }

  export -f run_antismash
  export CPUS ANTISMASH_CPUS BAKTA_DIR ANTISMASH_DIR ENV_ANTISMASH ANTISMASH_TAXON ANTISMASH_EXTRA_ARGS

  # Progress bar for Stage 4: each parallel job appends a line to a counter file.
  # Uses /dev/tty for direct terminal output since parallel subshells don't inherit fd 3.
  local_total=${#ALL_SAMPLES[@]}
  _AS_COUNTER="$OUTPUT_DIR/.antismash_progress_$$"
  > "$_AS_COUNTER"
  export _AS_COUNTER local_total _PB_WIDTH

  # Wrap run_antismash to tick the counter after each sample
  run_antismash_tracked() {
    run_antismash "$@"
    local rc=$?
    echo x >> "$_AS_COUNTER"
    local done
    done=$(wc -l < "$_AS_COUNTER")
    local filled=$(( done * _PB_WIDTH / local_total ))
    local empty=$(( _PB_WIDTH - filled ))
    local bar
    bar="$(printf '%*s' "$filled" '' | tr ' ' '#')$(printf '%*s' "$empty" '' | tr ' ' '-')"
    if (( done >= local_total )); then
      printf '\r  [%s] %d/%d Stage 4 (antismash)\n' "$bar" "$done" "$local_total" > /dev/tty
    else
      printf '\r  [%s] %d/%d Stage 4 (antismash)' "$bar" "$done" "$local_total" > /dev/tty
    fi
    return $rc
  }
  export -f run_antismash_tracked

  progress_bar 0 "$local_total" "Stage 4 (antismash)"

  parallel -j "$ANTISMASH_JOBS" --linebuffer \
    run_antismash_tracked {1} "$BAKTA_DIR"/{1} "$ANTISMASH_DIR" \
    ::: "${ALL_SAMPLES[@]}" || true   # failures handled by index.html scan below

  rm -f "$_AS_COUNTER"

  # Detect antiSMASH failures — parallel subshells can't write to FAILED_SAMPLES,
  # so we scan for samples that have a Bakta gbff but no antiSMASH index.html
  for sample in "${ALL_SAMPLES[@]}"; do
    gbff="$BAKTA_DIR/$sample/${sample}.gbff"
    out="$ANTISMASH_DIR/$sample"
    if [[ -f "$gbff" && ! -f "$out/index.html" ]]; then
      eprint "ERROR [Stage 4]: antiSMASH appears to have failed for '$sample' (no index.html)"
      FAILED_SAMPLES+=("$sample (antismash)")
    fi
  done
fi  # end DRY_RUN check for stage 4

fi  # end skip-antismash check

# ==============================================================================
# STAGE 5: 16S rRNA identification
# Barrnap extracts 16S sequences from each assembly, then Kraken2 classifies
# them against the SILVA database. Top hit (majority vote across copies) is
# reported as species, genus, and family in the summary.
# Runs on ALL_SAMPLES only — filtered samples lack reliable annotation.
# ==============================================================================
tprint ""
tprint "=============================="
tprint " STAGE 5: 16S Species ID"
tprint "=============================="

if $SKIP_16S; then
  tprint "  [skipped] --skip-16s flag set"
else

if $DRY_RUN; then
  for sample in "${ALL_SAMPLES[@]}"; do
    tprint "  [dry-run] barrnap + kraken2/SILVA → $RRNA_DIR/$sample"
  done
else
  # Load associative arrays for 16S results (used in summary stage)
  declare -A RRNA_SPECIES RRNA_GENUS RRNA_FAMILY RRNA_N_SEQS RRNA_N_GENERA RRNA_MIXED

  for sample in "${ALL_SAMPLES[@]}"; do
    rrna_outdir="$RRNA_DIR/$sample"
    rrna_done="$rrna_outdir/kraken2_report.txt"

    # Resume: skip if kraken2 report already exists
    if [[ -f "$rrna_done" ]]; then
      tprint "Skip 16S (already done): $sample"
      # Load existing results for summary
      if [[ -f "$rrna_outdir/classification.tsv" ]]; then
        IFS=$'\t' read -r _ sp ge fa n_seqs n_gen mflag _ \
          < <(tail -1 "$rrna_outdir/classification.tsv" 2>/dev/null || true)
        RRNA_SPECIES["$sample"]="${sp:-NA}"
        RRNA_GENUS["$sample"]="${ge:-NA}"
        RRNA_FAMILY["$sample"]="${fa:-NA}"
        RRNA_N_SEQS["$sample"]="${n_seqs:-NA}"
        RRNA_N_GENERA["$sample"]="${n_gen:-NA}"
        RRNA_MIXED["$sample"]="${mflag:-NA}"
      fi
      continue
    fi

    # Find assembly — either from assembly stage or --assemblies dir
    if [[ -n "$EXISTING_ASSEMBLIES_DIR" ]]; then
      assembly_fasta=$(find_existing_fasta "$sample" "$EXISTING_ASSEMBLIES_DIR")
    else
      assembly_fasta="$ASSEMBLY_DIR/$sample/final.contigs.fa"
    fi

    if [[ ! -f "$assembly_fasta" ]]; then
      tprint "Warning [Stage 5]: No assembly found for '$sample' — skipping 16S ID."
      RRNA_SPECIES["$sample"]="NA"; RRNA_GENUS["$sample"]="NA"; RRNA_FAMILY["$sample"]="NA"; RRNA_N_SEQS["$sample"]="NA"; RRNA_N_GENERA["$sample"]="NA"; RRNA_MIXED["$sample"]="NA"
      FAILED_SAMPLES+=("$sample (16S: no assembly)")
      continue
    fi

    mkdir -p "$rrna_outdir"

    # Step 1: barrnap — extract 16S rRNA sequences
    rrna_fasta="$rrna_outdir/rrna.fasta"
    tprint "Running barrnap: $sample"
    if ! conda run -p "$ENV_BARRNAP" --no-capture-output \
      barrnap \
        --kingdom "$BARRNAP_KINGDOM" \
        --reject "$BARRNAP_MIN_SCORE" \
        --threads "$CPUS" \
        --outseq "$rrna_fasta" \
        "$assembly_fasta" \
        > "$rrna_outdir/barrnap.gff" 2>&1; then
      tprint "Warning [Stage 5]: barrnap failed for '$sample' — skipping 16S ID."
      RRNA_SPECIES["$sample"]="NA"; RRNA_GENUS["$sample"]="NA"; RRNA_FAMILY["$sample"]="NA"; RRNA_N_SEQS["$sample"]="NA"; RRNA_N_GENERA["$sample"]="NA"; RRNA_MIXED["$sample"]="NA"
      FAILED_SAMPLES+=("$sample (16S: barrnap failed)")
      continue
    fi

    # Filter to 16S only (barrnap outputs all rRNA types)
    rrna_16s_fasta="$rrna_outdir/16S.fasta"
    grep -A1 "16S" "$rrna_fasta" | grep -v "^--$" > "$rrna_16s_fasta" 2>/dev/null || true

    if [[ ! -s "$rrna_16s_fasta" ]]; then
      tprint "Warning [Stage 5]: No 16S rRNA sequences found in '$sample' — skipping classification."
      RRNA_SPECIES["$sample"]="NA"; RRNA_GENUS["$sample"]="NA"; RRNA_FAMILY["$sample"]="NA"; RRNA_N_SEQS["$sample"]="NA"; RRNA_N_GENERA["$sample"]="NA"; RRNA_MIXED["$sample"]="NA"
      # Not added to FAILED_SAMPLES — fragmented assemblies may legitimately lack 16S
      # Write an empty report so resume logic skips this sample next time
      touch "$rrna_done"
      continue
    fi

    # Step 2: Kraken2 against SILVA
    kraken2_output="$rrna_outdir/kraken2_output.txt"
    tprint "Running Kraken2/SILVA 16S: $sample"
    if ! conda run -p "$ENV_BARRNAP" --no-capture-output \
      kraken2 \
        --db "$SILVA_DB" \
        --threads "$CPUS" \
        --use-names \
        --report "$rrna_done" \
        --output "$kraken2_output" \
        "$rrna_16s_fasta" 2>&1; then
      tprint "Warning [Stage 5]: Kraken2 failed for '$sample' — skipping 16S classification."
      RRNA_SPECIES["$sample"]="NA"; RRNA_GENUS["$sample"]="NA"; RRNA_FAMILY["$sample"]="NA"; RRNA_N_SEQS["$sample"]="NA"; RRNA_N_GENERA["$sample"]="NA"; RRNA_MIXED["$sample"]="NA"
      FAILED_SAMPLES+=("$sample (16S: kraken2 failed)")
      continue
    fi

    # Step 3: Parse top hit from Kraken2 report
    # Report format: %covered  reads_covered  reads_assigned  rank  taxid  name
    # Rank codes: S=species, G=genus, F=family
    # Take the highest-covered hit at each rank level.
    parse_kraken_rank() {
      local report="$1" rank_code="$2"
      grep $'\t'"${rank_code}"$'\t' "$report" 2>/dev/null \
        | sort -k1,1rn \
        | head -1 \
        | awk -F'\t' '{gsub(/^[[:space:]]+|[[:space:]]+$/, "", $6); print $6}'
    }

    top_species=$(parse_kraken_rank "$rrna_done" "S")
    top_genus=$(parse_kraken_rank   "$rrna_done" "G")
    top_family=$(parse_kraken_rank  "$rrna_done" "F")

    RRNA_SPECIES["$sample"]="${top_species:-NA}"
    RRNA_GENUS["$sample"]="${top_genus:-NA}"
    RRNA_FAMILY["$sample"]="${top_family:-NA}"

    # Step 4: Mixed culture detection
    # The kraken2 output file (--use-names) has one line per 16S sequence:
    #   C/U  seqid  taxonomy_string  length  kmer_hits
    # Taxonomy string format: "Nocardia farcinica (taxid 1234)"
    # We extract the genus from each classified sequence and check for conflicts.
    #
    # Conflict levels:
    #   OK             — all sequences agree at genus level (or only 1 sequence)
    #   genus_conflict — sequences disagree at genus but agree at family
    #   family_conflict — sequences disagree at family level (strong mixed signal)
    #   unclassified   — all sequences unclassified
    n_16s_seqs=$(grep -c '>' "$rrna_16s_fasta" 2>/dev/null || echo 0)
    mixed_flag="OK"
    n_genera=0

    if [[ -f "$kraken2_output" ]] && [[ "$n_16s_seqs" -gt 1 ]]; then
      # Extract genus from each classified sequence.
      # Kraken2 --use-names puts the full taxonomy path in field 3, e.g.:
      #   "Bacteria;Actinobacteria;Nocardia farcinica (taxid 1234)"
      # We want the genus — the first word of the species name is typically genus,
      # but the taxon field may contain higher ranks. The report genus lines (rank G)
      # are reliable — map each sequence's taxid back via the report is complex, so
      # instead we use the per-sequence name and extract the likely genus as the
      # first word of the deepest-level name (after the last semicolon).
      mapfile -t seq_genera < <(
        awk -F'\t' '$1=="C" {
          # field 3 is the taxonomy name string
          name = $3
          # strip taxid annotation "(taxid NNNN)"
          gsub(/[[:space:]]*\(taxid[[:space:]]+[0-9]+\)/, "", name)
          # take the last semicolon-delimited segment
          n = split(name, parts, ";")
          last = parts[n]
          # strip leading/trailing whitespace
          gsub(/^[[:space:]]+|[[:space:]]+$/, "", last)
          # first word is genus
          split(last, words, " ")
          if (words[1] != "" && words[1] != "unclassified")
            print words[1]
        }' "$kraken2_output" 2>/dev/null | sort -u
      )
      n_genera=${#seq_genera[@]}

      if [[ "$n_genera" -eq 0 ]]; then
        mixed_flag="unclassified"
      elif [[ "$n_genera" -gt 1 ]]; then
        # Check if the conflicting genera share the same family
        # Use the report's family-level lines to see how many families are represented
        n_families=$(grep $'\t'"F"$'\t' "$rrna_done" 2>/dev/null \
          | awk -F'\t' '$1 > 0' \
          | wc -l | tr -d ' ')
        if [[ "$n_families" -gt 1 ]]; then
          mixed_flag="family_conflict"
        else
          mixed_flag="genus_conflict"
        fi
      fi
    fi

    RRNA_N_SEQS["$sample"]="$n_16s_seqs"
    RRNA_N_GENERA["$sample"]="$n_genera"
    RRNA_MIXED["$sample"]="$mixed_flag"

    # Write classification TSV (extended with mixed culture columns)
    printf 'sample\tspecies\tgenus\tfamily\t16S_copies\tn_genera\tmixed_flag\tdb\n' \
      > "$rrna_outdir/classification.tsv"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tSILVA\n' \
      "$sample" \
      "${top_species:-NA}" \
      "${top_genus:-NA}" \
      "${top_family:-NA}" \
      "$n_16s_seqs" \
      "$n_genera" \
      "$mixed_flag" \
      >> "$rrna_outdir/classification.tsv"

    tprint "16S done: $sample → ${top_species:-NA} (genus: ${top_genus:-NA}) | copies: $n_16s_seqs | genera: $n_genera | mixed: $mixed_flag"
  done
fi  # end DRY_RUN check for stage 5

fi  # end skip-16s check

# ==============================================================================
# STAGE 6: GTDBTk taxonomy (classify_wf using assemblies directly)
# Runs on ALL_SAMPLES + FILTERED_SAMPLES — taxonomy doesn't care about
# fragmentation. We combine both lists into a deduplicated GTDBTK_SAMPLES array.
# ==============================================================================
if ! $DRY_RUN && $WATCHDOG_ENABLED; then
  stop_watchdog
fi

tprint ""
tprint "=============================="
tprint " STAGE 6: GTDBTk Taxonomy"
tprint "=============================="

if $SKIP_GTDBTK; then
  tprint "  [skipped] --skip-gtdbtk flag set"
else

# Build the full sample list for GTDBTk (ALL + FILTERED, deduplicated)
declare -A _gtdbtk_seen
GTDBTK_SAMPLES=()
for s in "${ALL_SAMPLES[@]}" "${FILTERED_SAMPLES[@]}"; do
  if [[ -z "${_gtdbtk_seen[$s]+x}" ]]; then
    _gtdbtk_seen[$s]=1
    GTDBTK_SAMPLES+=("$s")
  fi
done
unset _gtdbtk_seen

GTDBTK_SUMMARY="$GTDBTK_DIR/gtdbtk.bac120.summary.tsv"
GTDBTK_BATCHFILE="$GTDBTK_DIR/batchfile.tsv"

if $DRY_RUN; then
  tprint "  [dry-run] Would run GTDBTk classify_wf on ${#GTDBTK_SAMPLES[@]} samples → $GTDBTK_DIR"
  tprint "  [dry-run] skani sketch dir: ${GTDBTK_SKANI_SKETCH_DIR:-<not set, using default>}"
  for sample in "${GTDBTK_SAMPLES[@]}"; do
    tprint "    $ASSEMBLY_DIR/$sample/final.contigs.fa  $sample"
  done
else
  # Resume: skip if bacterial summary already exists
  if [[ -f "$GTDBTK_SUMMARY" ]]; then
    tprint "Skipping GTDBTk — output already exists ($GTDBTK_SUMMARY)"
  else
    # Build batchfile: <fasta_path>\t<genome_id>
    > "$GTDBTK_BATCHFILE"
    for sample in "${GTDBTK_SAMPLES[@]}"; do
      fasta="$ASSEMBLY_DIR/$sample/final.contigs.fa"
      if [[ -f "$fasta" ]]; then
        printf '%s\t%s\n' "$fasta" "$sample" >> "$GTDBTK_BATCHFILE"
      else
        tprint "Warning [Stage 5]: No assembly found for '$sample' — excluding from GTDBTk."
      fi
    done

    tprint "Running GTDBTk classify_wf on $(wc -l < "$GTDBTK_BATCHFILE") genomes..."

    # Build env override for GTDBTK_DATA_PATH if set in conf —
    # overrides whatever DB path is baked into the conda environment.
    gtdbtk_env_override=()
    if [[ -n "${GTDBTK_DATA_PATH:-}" ]]; then
      gtdbtk_env_override=(--no-capture-output env GTDBTK_DATA_PATH="$GTDBTK_DATA_PATH")
      tprint "GTDBTk: using GTDBTK_DATA_PATH=$GTDBTK_DATA_PATH"
    else
      gtdbtk_env_override=(--no-capture-output)
    fi

    if ! conda run -p "$ENV_GTDBTK" "${gtdbtk_env_override[@]}" \
      gtdbtk classify_wf \
        --batchfile "$GTDBTK_BATCHFILE" \
        --out_dir "$GTDBTK_DIR" \
        --cpus "$CPUS" \
        --force \
        ${GTDBTK_SKANI_SKETCH_DIR:+--skani_sketch_dir "$GTDBTK_SKANI_SKETCH_DIR"}; then
      eprint "ERROR [Stage 6]: GTDBTk classify_wf failed — check log for details."
      FAILED_SAMPLES+=("GTDBTk (all samples)")
    else
      tprint "GTDBTk complete."
    fi
  fi
fi  # end DRY_RUN check for stage 6

fi  # end skip-gtdbtk check
# ==============================================================================
tprint ""
tprint "=============================="
tprint " Generating summary TSV"
tprint "=============================="

if $DRY_RUN; then
  tprint "  [dry-run] Would write summary to $SUMMARY_TSV"
else
  # --------------------------------------------------------------------------
  # Helper: extract BGC types from a single antiSMASH region .gbk file.
  # antiSMASH encodes cluster type in COMMENT or structured comment lines like:
  #   ##antiSMASH-Data-START##
  #   ...
  #   ClusterType: NRPS
  # or in the /product qualifier of CDS features:
  #   /product="NRPS"
  # We parse the "ClusterType:" line from the structured comment block.
  # PKS subtypes appear as transAT-PKS, PKS-I, PKS-like etc — we normalise these.
  # --------------------------------------------------------------------------
  get_bgc_types_from_gbk() {
    local gbk="$1"
    # antiSMASH 6+ (including 8) stores BGC type in the /product qualifier
    # of the 'region' feature in each region GBK file, e.g.:
    #   region          1..45678
    #                   /product="NRPS"
    # There may be multiple /product lines if the region contains hybrid clusters.
    # We extract the value, split on commas (hybrid types), normalise PKS subtypes.
    # grep exit code 1 (no match) is normal — || true prevents pipefail abort.
    awk '
      /^     region / { in_region=1 }
      in_region && /\/product=/ {
        # Extract value between quotes
        match($0, /\/product="([^"]+)"/, arr)
        if (arr[1] != "") {
          n = split(arr[1], parts, ",")
          for (i=1; i<=n; i++) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", parts[i])
            if (parts[i] != "") print parts[i]
          }
        }
      }
      # Stop after region feature block ends (next feature keyword)
      in_region && /^     [a-zA-Z]/ && !/^     region / { in_region=0 }
    ' "$gbk" 2>/dev/null \
      | sed \
          -e 's/transAT-PKS.*/PKS-I/g' \
          -e 's/T1PKS/PKS-I/g'   -e 's/PKS1/PKS-I/g' \
          -e 's/T2PKS/PKS-II/g'  -e 's/PKS2/PKS-II/g' \
          -e 's/T3PKS\|PKS3\|PKS-III/PKS-III/g' \
      | grep -v '^$' || true
  }

  # --------------------------------------------------------------------------
  # Pass 1: collect all BGC types observed across the whole dataset
  # --------------------------------------------------------------------------
  declare -A ALL_BGC_TYPES
  declare -A SAMPLE_BGC_DATA  # sample -> space-separated list of types (with repeats)

  for sample in "${ALL_SAMPLES[@]}"; do
    antismash_out="$ANTISMASH_DIR/$sample"
    type_list=""
    if [[ -d "$antismash_out" ]]; then
      while IFS= read -r gbk; do
        while IFS= read -r bgc_type; do
          [[ -z "$bgc_type" ]] && continue
          ALL_BGC_TYPES["$bgc_type"]=1
          type_list+="$bgc_type "
        done < <(get_bgc_types_from_gbk "$gbk")
      done < <(find "$antismash_out" -maxdepth 1 -name "*.region*.gbk")
    fi
    SAMPLE_BGC_DATA["$sample"]="${type_list% }"
  done

  # Sort BGC type columns: NRPS first, then PKS-I/II/III, then rest alphabetically
  mapfile -t SORTED_BGC_TYPES < <(
    for t in "${!ALL_BGC_TYPES[@]}"; do echo "$t"; done | sort | \
    awk 'BEGIN{
      order["NRPS"]=1; order["PKS-I"]=2; order["PKS-II"]=3; order["PKS-III"]=4
    }
    {
      if ($0 in order) printf "%d\t%s\n", order[$0], $0
      else             printf "9\t%s\n", $0
    }' | sort -k1,1n -k2,2 | cut -f2
  )

  # --------------------------------------------------------------------------
  # Pass 2: write the TSV
  # --------------------------------------------------------------------------

  # Build header
  header="sample\tgtdb_classification\tgtdb_genus\tgtdb_species\t16S_species\t16S_genus\t16S_family\t16S_copies\t16S_n_genera\t16S_mixed_flag\tnum_contigs\tassembly_size_bp\tN50\tN90\tL50\tbakta_CDS\tantismash_BGCs_total\tantismash_BGCs_on_edge"
  for t in "${SORTED_BGC_TYPES[@]}"; do
    header+="\t${t}"
  done
  printf '%b\n' "$header" > "$SUMMARY_TSV"

  # Load 16S results — reload from classification.tsv files if not already in memory
  # (e.g. when Stage 5 was fully skipped on resume)
  declare -A RRNA_SPECIES RRNA_GENUS RRNA_FAMILY RRNA_N_SEQS RRNA_N_GENERA RRNA_MIXED 2>/dev/null || true
  for sample in "${ALL_SAMPLES[@]}"; do
    if [[ -z "${RRNA_SPECIES[$sample]+x}" ]]; then
      classf="$RRNA_DIR/$sample/classification.tsv"
      if [[ -f "$classf" ]]; then
        IFS=$'\t' read -r _ sp ge fa n_seqs n_gen mflag _ \
          < <(tail -1 "$classf" 2>/dev/null || true)
        RRNA_SPECIES["$sample"]="${sp:-NA}"
        RRNA_GENUS["$sample"]="${ge:-NA}"
        RRNA_FAMILY["$sample"]="${fa:-NA}"
        RRNA_N_SEQS["$sample"]="${n_seqs:-NA}"
        RRNA_N_GENERA["$sample"]="${n_gen:-NA}"
        RRNA_MIXED["$sample"]="${mflag:-NA}"
      else
        RRNA_SPECIES["$sample"]="NA"
        RRNA_GENUS["$sample"]="NA"
        RRNA_FAMILY["$sample"]="NA"
        RRNA_N_SEQS["$sample"]="NA"
        RRNA_N_GENERA["$sample"]="NA"
        RRNA_MIXED["$sample"]="NA"
      fi
    fi
  done

  # Load GTDBTk results into associative arrays keyed by sample name.
  # GTDBTk produces separate bacterial and archaeal summary files; merge both.
  declare -A GTDB_CLASSIFICATION GTDB_GENUS GTDB_SPECIES
  for gtdb_tsv in \
    "$GTDBTK_DIR/gtdbtk.bac120.summary.tsv" \
    "$GTDBTK_DIR/gtdbtk.ar53.summary.tsv"; do
    if [[ -f "$gtdb_tsv" ]]; then
      # Header: user_genome classification ...
      # Skip header line (starts with "user_genome")
      while IFS=$'\t' read -r genome classification rest; do
        [[ "$genome" == "user_genome" ]] && continue
        [[ -z "$genome" ]] && continue
        GTDB_CLASSIFICATION["$genome"]="$classification"
        # Parse genus (g__...) and species (s__...) from semicolon-delimited string
        gtdb_genus=$(echo "$classification" | grep -oP 'g__[^;]+' | sed 's/g__//' || echo "")
        gtdb_species=$(echo "$classification" | grep -oP 's__[^;]+' | sed 's/s__//' || echo "")
        GTDB_GENUS["$genome"]="${gtdb_genus:-NA}"
        GTDB_SPECIES["$genome"]="${gtdb_species:-NA}"
      done < "$gtdb_tsv"
    fi
  done

  # Load seqkit assembly TSV column indices once
  declare -A _asm_col
  if [[ -f "$SEQKIT_ASSEMBLY_TSV" ]]; then
    IFS=$'\t' read -ra _asm_hdr < "$SEQKIT_ASSEMBLY_TSV"
    for i in "${!_asm_hdr[@]}"; do _asm_col["${_asm_hdr[$i]}"]=$(( i + 1 )); done
  fi

  for sample in "${ALL_SAMPLES[@]}"; do
    # Assembly stats from SeqKit TSV
    num_contigs="NA"; assembly_size="NA"; n50="NA"; n90="NA"; l50="NA"
    if [[ -f "$SEQKIT_ASSEMBLY_TSV" ]]; then
      seqkit_row=$(grep "^${sample}"$'\t' "$SEQKIT_ASSEMBLY_TSV" | head -1)
      if [[ -n "$seqkit_row" ]]; then
        num_contigs=$(echo "$seqkit_row" | cut -f"${_asm_col[num_contigs]:-2}")
        assembly_size=$(echo "$seqkit_row" | cut -f"${_asm_col[assembly_size_bp]:-3}")
        n50=$(echo "$seqkit_row"         | cut -f"${_asm_col[N50]:-7}")
        n90=$(echo "$seqkit_row"         | cut -f"${_asm_col[N90]:-8}")
        l50=$(echo "$seqkit_row"         | cut -f"${_asm_col[L50]:-9}")
      fi
    fi

    # Bakta CDS count
    bakta_cds="NA"
    bakta_txt="$BAKTA_DIR/$sample/${sample}.txt"
    if [[ -f "$bakta_txt" ]]; then
      bakta_cds=$(grep -i "^CDS:" "$bakta_txt" | awk '{print $2}' | head -1)
      [[ -z "$bakta_cds" ]] && bakta_cds="NA"
    fi

    # antiSMASH: total BGCs, edge BGCs, per-type counts
    bgc_total="NA"; bgc_on_edge="NA"
    declare -A bgc_type_counts
    for t in "${SORTED_BGC_TYPES[@]}"; do bgc_type_counts["$t"]=0; done

    antismash_out="$ANTISMASH_DIR/$sample"
    if [[ -d "$antismash_out" ]]; then
      bgc_total=0; bgc_on_edge=0
      while IFS= read -r gbk; do
        (( bgc_total++ )) || true
        # Check if on contig edge
        if grep -q "contig_edge.*True\|Orig. start.*<\|Orig\. end.*>" "$gbk" 2>/dev/null; then
          (( bgc_on_edge++ )) || true
        fi
        # Count BGC types for this region
        while IFS= read -r bgc_type; do
          [[ -z "$bgc_type" ]] && continue
          if [[ -v bgc_type_counts["$bgc_type"] ]]; then
            (( bgc_type_counts["$bgc_type"]++ )) || true
          fi
        done < <(get_bgc_types_from_gbk "$gbk")
      done < <(find "$antismash_out" -maxdepth 1 -name "*.region*.gbk")
    fi

    # Build row — use printf '%b' only for the tab-separated structure,
    # then append a newline safely. Fields are controlled values (numbers/NA/FILTERED)
    # except gtdb strings which could contain backslashes — sanitise those.
    gtdb_class="${GTDB_CLASSIFICATION[$sample]:-NA}"; gtdb_class="${gtdb_class//\\/}"
    gtdb_gen="${GTDB_GENUS[$sample]:-NA}";            gtdb_gen="${gtdb_gen//\\/}"
    gtdb_sp="${GTDB_SPECIES[$sample]:-NA}";           gtdb_sp="${gtdb_sp//\\/}"
    rrna_sp="${RRNA_SPECIES[$sample]:-NA}";           rrna_sp="${rrna_sp//\\/}"
    rrna_ge="${RRNA_GENUS[$sample]:-NA}";             rrna_ge="${rrna_ge//\\/}"
    rrna_fa="${RRNA_FAMILY[$sample]:-NA}";            rrna_fa="${rrna_fa//\\/}"
    rrna_copies="${RRNA_N_SEQS[$sample]:-NA}"
    rrna_n_gen="${RRNA_N_GENERA[$sample]:-NA}"
    rrna_mixed="${RRNA_MIXED[$sample]:-NA}"
    row="$sample\t$gtdb_class\t$gtdb_gen\t$gtdb_sp\t$rrna_sp\t$rrna_ge\t$rrna_fa\t$rrna_copies\t$rrna_n_gen\t$rrna_mixed\t$num_contigs\t$assembly_size\t$n50\t$n90\t$l50\t$bakta_cds\t$bgc_total\t$bgc_on_edge"
    for t in "${SORTED_BGC_TYPES[@]}"; do
      row+="\t${bgc_type_counts[$t]:-0}"
    done
    printf '%b\n' "$row" >> "$SUMMARY_TSV"

    unset bgc_type_counts
    declare -A bgc_type_counts
  done

  # Write rows for filtered samples — have taxonomy but no Bakta/antiSMASH
  for sample in "${FILTERED_SAMPLES[@]}"; do
    seqkit_row=$(grep "^${sample}"$'\t' "$SEQKIT_ASSEMBLY_TSV" 2>/dev/null | head -1)
    num_contigs="NA"; assembly_size="NA"; n50="NA"; n90="NA"; l50="NA"
    if [[ -n "$seqkit_row" ]]; then
      num_contigs=$(echo "$seqkit_row"  | cut -f"${_asm_col[num_contigs]:-2}")
      assembly_size=$(echo "$seqkit_row" | cut -f"${_asm_col[assembly_size_bp]:-3}")
      n50=$(echo "$seqkit_row"          | cut -f"${_asm_col[N50]:-7}")
      n90=$(echo "$seqkit_row"          | cut -f"${_asm_col[N90]:-8}")
      l50=$(echo "$seqkit_row"          | cut -f"${_asm_col[L50]:-9}")
    fi
    gtdb_class="${GTDB_CLASSIFICATION[$sample]:-NA}"; gtdb_class="${gtdb_class//\\/}"
    gtdb_gen="${GTDB_GENUS[$sample]:-NA}";            gtdb_gen="${gtdb_gen//\\/}"
    gtdb_sp="${GTDB_SPECIES[$sample]:-NA}";           gtdb_sp="${gtdb_sp//\\/}"
    row="$sample\t$gtdb_class\t$gtdb_gen\t$gtdb_sp\tFILTERED\tFILTERED\tFILTERED\t$num_contigs\t$assembly_size\t$n50\t$n90\t$l50\tFILTERED\tFILTERED\tFILTERED"
    for t in "${SORTED_BGC_TYPES[@]}"; do row+="\tFILTERED"; done
    printf '%b\n' "$row" >> "$SUMMARY_TSV"
  done

  tprint "Summary written to: $SUMMARY_TSV"
  tprint "  BGC types found: ${SORTED_BGC_TYPES[*]:-none}"

  # --------------------------------------------------------------------------
  # Generate HTML summary
  # --------------------------------------------------------------------------
  SUMMARY_HTML="${SUMMARY_DIR}/summary.html"

  # html_cell: escape a value for safe HTML output
  html_esc() { printf '%s' "$1" | sed 's/&/\&amp;/g;s/</\&lt;/g;s/>/\&gt;/g;s/"/\&quot;/g'; }

  {
    # Relative path from 07_summary/ up to 04_antismash/
    AS_REL="../04_antismash"

    cat << 'HTMLEOF'
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>dragonSMASH Summary</title>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
         font-size: 13px; background: #f5f5f5; color: #222; padding: 16px; }
  h1 { font-size: 20px; margin-bottom: 4px; color: #1a1a2e; }
  .subtitle { color: #666; font-size: 12px; margin-bottom: 14px; }
  .controls { display: flex; gap: 10px; align-items: center; margin-bottom: 10px; flex-wrap: wrap; }
  .controls input { padding: 6px 10px; border: 1px solid #ccc; border-radius: 4px;
                    font-size: 13px; width: 260px; }
  .controls label { font-size: 12px; color: #555; display: flex; align-items: center; gap: 4px; }
  .wrap { overflow-x: auto; border-radius: 6px; box-shadow: 0 1px 4px rgba(0,0,0,.12); }
  table { border-collapse: collapse; background: #fff; min-width: 100%; white-space: nowrap; }
  th { background: #1a1a2e; color: #fff; padding: 8px 10px; text-align: left;
       cursor: pointer; user-select: none; font-weight: 600; font-size: 12px;
       position: sticky; top: 0; z-index: 1; }
  th:hover { background: #2e2e5e; }
  th.sort-asc::after  { content: " ↑"; opacity: .7; }
  th.sort-desc::after { content: " ↓"; opacity: .7; }
  td { padding: 6px 10px; border-bottom: 1px solid #ebebeb; vertical-align: middle; }
  tr:hover td { background: #f0f4ff; }
  tr.filtered-row td { background: #f8f8f8; color: #aaa; font-style: italic; }
  tr.hidden { display: none; }
  a { color: #2255cc; text-decoration: none; }
  a:hover { text-decoration: underline; }
  .badge { display: inline-block; padding: 2px 7px; border-radius: 10px;
           font-size: 11px; font-weight: 600; }
  .badge-ok     { background: #d4edda; color: #155724; }
  .badge-genus  { background: #fff3cd; color: #856404; }
  .badge-family { background: #f8d7da; color: #721c24; }
  .badge-na     { background: #e9ecef; color: #6c757d; }
  .badge-unc    { background: #e2e3e5; color: #383d41; }
  .num { text-align: right; font-variant-numeric: tabular-nums; }
  .bgc-nonzero { font-weight: 700; color: #8b0000; }
  .group-sep { border-left: 2px solid #ddd; }
  .col-toggle { font-size: 11px; }
  .stats { display: flex; gap: 16px; margin-bottom: 12px; flex-wrap: wrap; }
  .stat-box { background: #fff; border-radius: 6px; padding: 10px 16px;
              box-shadow: 0 1px 3px rgba(0,0,0,.1); min-width: 110px; }
  .stat-box .val { font-size: 22px; font-weight: 700; color: #1a1a2e; }
  .stat-box .lbl { font-size: 11px; color: #888; margin-top: 2px; }
</style>
</head>
<body>
HTMLEOF

    # ── Header ──────────────────────────────────────────────────────────────
    run_date=$(date '+%d %B %Y, %H:%M')
    n_samples=${#ALL_SAMPLES[@]}
    n_filtered=${#FILTERED_SAMPLES[@]}
    n_failed=${#FAILED_SAMPLES[@]}
    total_bgcs=0
    for s in "${ALL_SAMPLES[@]}"; do
      v=$(grep "^${s}"$'\t' "$SUMMARY_TSV" 2>/dev/null | cut -f17)
      [[ "$v" =~ ^[0-9]+$ ]] && (( total_bgcs += v )) || true
    done

    printf '<h1>dragonSMASH Summary</h1>\n'
    printf '<div class="subtitle">Generated %s &nbsp;|&nbsp; Pipeline v%s</div>\n' \
      "$run_date" "$(grep '^# Version' "${BASH_SOURCE[0]}" | tail -1 | awk '{print $NF}')"

    # Stat boxes
    printf '<div class="stats">\n'
    printf '  <div class="stat-box"><div class="val">%s</div><div class="lbl">Samples</div></div>\n' "$n_samples"
    printf '  <div class="stat-box"><div class="val">%s</div><div class="lbl">Filtered</div></div>\n' "$n_filtered"
    printf '  <div class="stat-box"><div class="val">%s</div><div class="lbl">Failed</div></div>\n' "$n_failed"
    printf '  <div class="stat-box"><div class="val">%s</div><div class="lbl">Total BGCs</div></div>\n' "$total_bgcs"
    printf '  <div class="stat-box"><div class="val">%s</div><div class="lbl">BGC types</div></div>\n' "${#SORTED_BGC_TYPES[@]}"
    printf '</div>\n'

    # Controls
    printf '<div class="controls">\n'
    printf '  <input type="text" id="search" placeholder="Search samples…" oninput="filterTable()">\n'
    printf '  <label><input type="checkbox" id="hideFiltered" onchange="filterTable()"> Hide filtered samples</label>\n'
    printf '  <label><input type="checkbox" id="hideBGCzero" onchange="filterTable()"> Hide samples with 0 BGCs</label>\n'
    printf '</div>\n'

    printf '<div class="wrap"><table id="tbl">\n<thead><tr>\n'

    # ── Column headers ───────────────────────────────────────────────────────
    # Group 1: identity
    for col in "Sample" "GTDB Species" "GTDB Genus" "GTDB Classification" \
               "16S Species" "16S Genus" "16S Family" "16S Copies" "16S Genera" "Mixed Flag"; do
      printf '  <th onclick="sortTable(this)">%s</th>\n' "$col"
    done
    # Group 2: assembly (left border separator)
    printf '  <th class="group-sep" onclick="sortTable(this)">Contigs</th>\n'
    for col in "Size (bp)" "N50" "N90" "L50"; do
      printf '  <th onclick="sortTable(this)">%s</th>\n' "$col"
    done
    # Group 3: annotation
    printf '  <th class="group-sep" onclick="sortTable(this)">CDS</th>\n'
    # Group 4: BGC summary
    printf '  <th class="group-sep" onclick="sortTable(this)">BGCs</th>\n'
    printf '  <th onclick="sortTable(this)">On-edge</th>\n'
    # Group 5: BGC types
    if [[ ${#SORTED_BGC_TYPES[@]} -gt 0 ]]; then
      first=true
      for t in "${SORTED_BGC_TYPES[@]}"; do
        if $first; then
          printf '  <th class="group-sep" onclick="sortTable(this)">%s</th>\n' "$(html_esc "$t")"
          first=false
        else
          printf '  <th onclick="sortTable(this)">%s</th>\n' "$(html_esc "$t")"
        fi
      done
    fi

    printf '</tr></thead>\n<tbody>\n'

    # ── Data rows — main samples ─────────────────────────────────────────────
    for sample in "${ALL_SAMPLES[@]}"; do
      # Re-derive values from TSV for consistency
      tsv_row=$(grep "^${sample}"$'\t' "$SUMMARY_TSV" | head -1)
      [[ -z "$tsv_row" ]] && continue
      IFS=$'\t' read -ra F <<< "$tsv_row"
      # TSV columns: 0=sample 1=gtdb_class 2=gtdb_genus 3=gtdb_species
      #              4=16S_species 5=16S_genus 6=16S_family
      #              7=16S_copies 8=16S_n_genera 9=16S_mixed_flag
      #              10=num_contigs 11=assembly_size 12=N50 13=N90 14=L50
      #              15=bakta_CDS 16=bgc_total 17=bgc_on_edge
      #              18+ = bgc types

      mixed="${F[9]:-NA}"
      case "$mixed" in
        OK)               badge_class="badge-ok";     badge_tip="All 16S rRNA gene copies in this assembly agree at the genus level — consistent with a pure culture." ;;
        genus_conflict)   badge_class="badge-genus";  badge_tip="The 16S rRNA gene copies in this assembly belong to different genera, but all fall within the same family. This may reflect natural variation between copies within one organism, or could indicate a mixed culture of closely related bacteria. Treat results with caution." ;;
        family_conflict)  badge_class="badge-family"; badge_tip="The 16S rRNA gene copies in this assembly span multiple families — a strong indicator that this sample contains more than one organism (mixed culture or contamination). Taxonomy and BGC results may be unreliable." ;;
        unclassified)     badge_class="badge-unc";    badge_tip="None of the 16S rRNA sequences could be classified against the SILVA database." ;;
        *)                badge_class="badge-na";     badge_tip="No 16S classification data available." ;;
      esac

      antismash_link="${AS_REL}/${sample}/index.html"
      bgc_total_val="${F[16]:-NA}"

      printf '<tr>\n'
      # Sample — link to antiSMASH if it exists (we emit the link; browser handles missing)
      printf '  <td><a href="%s" target="_blank">%s</a></td>\n' \
        "$(html_esc "$antismash_link")" "$(html_esc "$sample")"
      # Taxonomy
      printf '  <td>%s</td>\n'  "$(html_esc "${F[3]:-NA}")"   # GTDB species
      printf '  <td>%s</td>\n'  "$(html_esc "${F[2]:-NA}")"   # GTDB genus
      printf '  <td title="%s">%s</td>\n' \
        "$(html_esc "${F[1]:-NA}")" \
        "$(html_esc "$(echo "${F[1]:-NA}" | sed 's/.*s__//' | cut -c1-40)")"  # GTDB class truncated
      printf '  <td>%s</td>\n'  "$(html_esc "${F[4]:-NA}")"   # 16S species
      printf '  <td>%s</td>\n'  "$(html_esc "${F[5]:-NA}")"   # 16S genus
      printf '  <td>%s</td>\n'  "$(html_esc "${F[6]:-NA}")"   # 16S family
      printf '  <td class="num">%s</td>\n' "$(html_esc "${F[7]:-NA}")"  # 16S copies
      printf '  <td class="num">%s</td>\n' "$(html_esc "${F[8]:-NA}")"  # 16S n_genera
      printf '  <td><span class="badge %s" title="%s">%s</span></td>\n' "$badge_class" "$badge_tip" "$(html_esc "$mixed")"
      # Assembly
      printf '  <td class="num group-sep">%s</td>\n' "$(html_esc "${F[10]:-NA}")"
      printf '  <td class="num">%s</td>\n' "$(html_esc "${F[11]:-NA}")"
      printf '  <td class="num">%s</td>\n' "$(html_esc "${F[12]:-NA}")"
      printf '  <td class="num">%s</td>\n' "$(html_esc "${F[13]:-NA}")"
      printf '  <td class="num">%s</td>\n' "$(html_esc "${F[14]:-NA}")"
      # Annotation
      printf '  <td class="num group-sep">%s</td>\n' "$(html_esc "${F[15]:-NA}")"
      # BGC summary
      bgc_class="num group-sep"
      [[ "$bgc_total_val" =~ ^[1-9] ]] && bgc_class="num group-sep bgc-nonzero"
      printf '  <td class="%s">%s</td>\n' "$bgc_class" "$(html_esc "$bgc_total_val")"
      printf '  <td class="num">%s</td>\n' "$(html_esc "${F[17]:-NA}")"
      # BGC type columns
      col_idx=18
      first_bgc=true
      for t in "${SORTED_BGC_TYPES[@]}"; do
        val="${F[$col_idx]:-0}"
        cell_class="num"
        $first_bgc && cell_class="num group-sep" && first_bgc=false
        [[ "$val" =~ ^[1-9] ]] && cell_class+=" bgc-nonzero"
        printf '  <td class="%s">%s</td>\n' "$cell_class" "$(html_esc "$val")"
        (( col_idx++ )) || true
      done
      printf '</tr>\n'
    done

    # ── Filtered sample rows ─────────────────────────────────────────────────
    for sample in "${FILTERED_SAMPLES[@]}"; do
      tsv_row=$(grep "^${sample}"$'\t' "$SUMMARY_TSV" | head -1)
      [[ -z "$tsv_row" ]] && continue
      IFS=$'\t' read -ra F <<< "$tsv_row"
      printf '<tr class="filtered-row" data-filtered="1">\n'
      printf '  <td>%s</td>\n' "$(html_esc "$sample")"
      printf '  <td>%s</td>\n' "$(html_esc "${F[3]:-NA}")"
      printf '  <td>%s</td>\n' "$(html_esc "${F[2]:-NA}")"
      printf '  <td>%s</td>\n' "$(html_esc "${F[1]:-NA}")"
      # 16S / mixed columns — all FILTERED
      for _ in 4 5 6 7 8; do printf '  <td>—</td>\n'; done
      printf '  <td><span class="badge badge-na">filtered</span></td>\n'
      # Assembly
      printf '  <td class="num group-sep">%s</td>\n' "$(html_esc "${F[10]:-NA}")"
      for i in 11 12 13 14; do
        printf '  <td class="num">%s</td>\n' "$(html_esc "${F[$i]:-NA}")"
      done
      # Annotation + BGC — filtered
      printf '  <td class="num group-sep">—</td>\n'
      printf '  <td class="num group-sep">—</td>\n'
      printf '  <td class="num">—</td>\n'
      for _ in "${SORTED_BGC_TYPES[@]}"; do printf '  <td class="num">—</td>\n'; done
      printf '</tr>\n'
    done

    # ── JavaScript ───────────────────────────────────────────────────────────
    cat << 'JSEOF'
</tbody></table></div>
<script>
function filterTable() {
  const q = document.getElementById('search').value.toLowerCase();
  const hideFiltered = document.getElementById('hideFiltered').checked;
  const hideZero = document.getElementById('hideBGCzero').checked;
  document.querySelectorAll('#tbl tbody tr').forEach(tr => {
    const isFiltered = tr.dataset.filtered === '1';
    const text = tr.textContent.toLowerCase();
    const bgcCell = tr.querySelectorAll('td')[16];
    const bgcVal = bgcCell ? parseInt(bgcCell.textContent) || 0 : 0;
    const hide = (q && !text.includes(q))
              || (hideFiltered && isFiltered)
              || (hideZero && !isFiltered && bgcVal === 0);
    tr.classList.toggle('hidden', hide);
  });
}

function sortTable(th) {
  const table = document.getElementById('tbl');
  const tbody = table.querySelector('tbody');
  const idx = Array.from(th.parentElement.children).indexOf(th);
  const asc = th.classList.contains('sort-asc') ? -1 : 1;
  document.querySelectorAll('#tbl th').forEach(t => t.classList.remove('sort-asc','sort-desc'));
  th.classList.add(asc === 1 ? 'sort-asc' : 'sort-desc');
  const rows = Array.from(tbody.querySelectorAll('tr'));
  rows.sort((a, b) => {
    const av = a.children[idx]?.textContent.trim() ?? '';
    const bv = b.children[idx]?.textContent.trim() ?? '';
    const an = parseFloat(av), bn = parseFloat(bv);
    if (!isNaN(an) && !isNaN(bn)) return asc * (an - bn);
    return asc * av.localeCompare(bv);
  });
  rows.forEach(r => tbody.appendChild(r));
}
</script>
</body></html>
JSEOF
  } > "$SUMMARY_HTML"

  tprint "Summary HTML written to: $SUMMARY_HTML"
fi

# ==============================================================================
# Final summary
# ==============================================================================
PIPELINE_END=$(date +%s)
ELAPSED=$(( PIPELINE_END - PIPELINE_START ))
ELAPSED_FMT=$(printf '%dh %dm %ds' $((ELAPSED/3600)) $(( (ELAPSED%3600)/60 )) $((ELAPSED%60)))

tprint ""
tprint "=============================="
tprint " Pipeline complete!"
tprint "=============================="
tprint "  Samples        : ${#ALL_SAMPLES[@]}"
tprint "  Elapsed        : $ELAPSED_FMT"
tprint "  Assemblies     : $ASSEMBLY_DIR"
tprint "  Assembly stats : $SEQKIT_ASSEMBLY_TSV"
[[ -z "$EXISTING_ASSEMBLIES_DIR" ]] && tprint "  Read stats     : $SEQKIT_READS_TSV"
tprint "  Bakta          : $BAKTA_DIR"
tprint "  antiSMASH      : $ANTISMASH_DIR"
tprint "  16S ID         : $RRNA_DIR"
tprint "  GTDBTk         : $GTDBTK_DIR"
tprint "  Summary TSV    : $SUMMARY_TSV"
tprint "  Summary HTML   : ${SUMMARY_HTML:-${SUMMARY_DIR}/summary.html}"

if [[ ${#FAILED_SAMPLES[@]} -gt 0 ]]; then
  tprint ""
  tprint "=============================="
  tprint " FAILED SAMPLES (${#FAILED_SAMPLES[@]})"
  tprint "=============================="
  for f in "${FAILED_SAMPLES[@]}"; do
    tprint "  - $f"
  done
  tprint ""
  tprint "Re-run the pipeline on failed samples after resolving issues."
  tprint "Completed samples will be skipped automatically."
fi

if [[ ${#FILTERED_SAMPLES[@]} -gt 0 ]]; then
  tprint ""
  tprint "=============================="
  tprint " FILTERED SAMPLES (${#FILTERED_SAMPLES[@]}) — exceeded --max-contigs $MAX_CONTIGS"
  tprint "=============================="
  for f in "${FILTERED_SAMPLES[@]}"; do
    tprint "  - $f (${FILTERED_CONTIG_COUNTS[$f]:-?} contigs)"
  done
  tprint ""
  tprint "These samples were assembled and QC'd but skipped for annotation."
  tprint "Re-run with --max-contigs 0 to annotate them regardless."
fi

tprint ""
bprint "  \"Do not meddle in the affairs of dragons,"
bprint "   for you are crunchy and taste good with ketchup.\""
tprint ""
