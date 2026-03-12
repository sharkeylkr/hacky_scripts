#!/usr/bin/env bash
# ==============================================================================
# setup_dragonsmash.bash
# One-time setup for dragonSMASH on MDU servers (marvin / zaphod / trillian).
# Run this before your first dragonSMASH run.
#
# Usage:
#   bash setup_dragonsmash.bash
#
# What it does:
#   0. Creates ~/bin and symlink dragonSMASH → dragonSMASHv3.bash, adds to PATH
#   1. Creates the 16S_ID_env conda environment (barrnap + kraken2)
#   2. Verifies all other required conda environments exist
#   3. Checks parallel is installed (required for Bakta/antiSMASH parallelism)
#   4. Verifies shared databases are accessible (GTDBTk, SILVA)
#   5. Checks the shared GTDBTk env is functional (v2.3+, skani)
#   6. Checks for a Bakta DB and prints download instructions if missing
#
# All output is logged to ~/setup_dragonsmash_DDMMYYYY_HHMMSS.log
# ==============================================================================

# Do NOT use set -e here — we want to catch and report all errors, not abort.
set -uo pipefail

LOG_FILE="$HOME/setup_dragonsmash_$(date +%d%m%Y_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

# ------------------------------------------------------------------------------
# MDU server defaults — shared paths (same on marvin / zaphod / trillian)
# ------------------------------------------------------------------------------
SHARED_CONDA=/home/shared/conda/envs
GTDBTK_ENV="$SHARED_CONDA/gtdbtk"
GTDBTK_DB=/home/shared/db/gtdbtk/release226
GTDBTK_SKANI=$GTDBTK_DB/skani/sketches
SILVA_DB=/home/shared/db/kraken2/silva

ENV_16S=16S_ID_env

# Required per-user conda environments (excluding GTDBTk which is shared)
REQUIRED_ENVS=(dragonflye seqkit bakta antismash8)

# Bakta DB — each user needs their own copy (default location, may differ)
BAKTA_DB_DEFAULT="$HOME/dbs/bakta/db"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
timestamp() { date '+[%Y-%m-%d %H:%M:%S]'; }
log()  { echo "$(timestamp) $*"; }
ok()   { echo "$(timestamp)   ✓  $*"; }
warn() { echo "$(timestamp)   ⚠  $*"; WARNINGS+=("$*"); }
fail() { echo "$(timestamp)   ✗  $*"; ERRORS+=("$*"); }

WARNINGS=()
ERRORS=()

# Initialise conda in this shell
init_conda() {
  local conda_sh=""
  for candidate in \
      "$HOME/miniconda3/etc/profile.d/conda.sh" \
      "$HOME/anaconda3/etc/profile.d/conda.sh" \
      "/opt/conda/etc/profile.d/conda.sh"; do
    [[ -f "$candidate" ]] && { conda_sh="$candidate"; break; }
  done
  if [[ -z "$conda_sh" ]]; then
    fail "Could not find conda.sh — is conda installed and in PATH?"
    exit 1
  fi
  # shellcheck source=/dev/null
  source "$conda_sh"
}

# ------------------------------------------------------------------------------
# Header
# ------------------------------------------------------------------------------
log "=============================="
log " dragonSMASH Setup"
log " User: $(whoami)@$(hostname -s)"
log "=============================="
log "  Log file  : $LOG_FILE"
log ""

init_conda

# ==============================================================================
# Step 0: PATH and symlink setup
# ==============================================================================
log "[0/6] Setting up ~/bin and PATH..."

SCRIPT_PATH="$(realpath "${BASH_SOURCE[0]}")"
SCRIPT_DIR_REAL="$(dirname "$SCRIPT_PATH")"
DRAGONSMASHER="$SCRIPT_DIR_REAL/dragonSMASHv3.bash"
BIN_DIR="$HOME/bin"
SYMLINK="$BIN_DIR/dragonSMASH"

mkdir -p "$BIN_DIR"

if [[ -L "$SYMLINK" ]] && [[ "$(readlink "$SYMLINK")" == "$DRAGONSMASHER" ]]; then
  ok "Symlink already correct: $SYMLINK → $DRAGONSMASHER"
elif [[ -f "$DRAGONSMASHER" ]]; then
  ln -sf "$DRAGONSMASHER" "$SYMLINK"
  ok "Symlink created: $SYMLINK → $DRAGONSMASHER"
else
  fail "dragonSMASHv3.bash not found at $DRAGONSMASHER"
  fail "  Place dragonSMASHv3.bash in the same directory as this setup script."
fi

# Add ~/bin to PATH in shell rc files if not already present
for rcfile in "$HOME/.bashrc" "$HOME/.bash_profile"; do
  if [[ -f "$rcfile" ]]; then
    if grep -q 'PATH="$HOME/bin:$PATH"' "$rcfile" 2>/dev/null || \
       grep -q "PATH=\$HOME/bin:\$PATH" "$rcfile" 2>/dev/null || \
       grep -q '~/bin' "$rcfile" 2>/dev/null; then
      ok "~/bin already in PATH via $(basename "$rcfile")"
    else
      {
        echo ''
        echo '# dragonSMASH — add ~/bin to PATH'
        echo 'export PATH="$HOME/bin:$PATH"'
      } >> "$rcfile"
      ok "Added ~/bin to PATH in $rcfile"
      warn "Run 'source $rcfile' or log out/in for PATH to take effect."
    fi
  fi
done

log "[0/6] Done."
log ""

# ==============================================================================
# Step 1: Create 16S_ID_env (barrnap + kraken2)
# ==============================================================================
log "[1/6] Setting up 16S_ID_env (barrnap + kraken2)..."

if conda env list | awk '{print $1}' | grep -qx "$ENV_16S"; then
  ok "Environment '$ENV_16S' already exists."

  if conda run -n "$ENV_16S" barrnap --version &>/dev/null; then
    ok "barrnap found in $ENV_16S."
  else
    warn "barrnap not found in $ENV_16S — attempting install..."
    conda install -n "$ENV_16S" -c bioconda barrnap -y \
      && ok "barrnap installed." \
      || fail "barrnap install failed — try manually: conda install -n $ENV_16S -c bioconda barrnap"
  fi

  if conda run -n "$ENV_16S" kraken2 --version &>/dev/null; then
    ok "kraken2 found in $ENV_16S."
  else
    warn "kraken2 not found in $ENV_16S — attempting install..."
    conda install -n "$ENV_16S" -c bioconda kraken2 -y \
      && ok "kraken2 installed." \
      || fail "kraken2 install failed — try manually: conda install -n $ENV_16S -c bioconda kraken2"
  fi

else
  log "      Creating '$ENV_16S' with barrnap + kraken2 (this may take a few minutes)..."
  if conda create -n "$ENV_16S" -c bioconda -c conda-forge barrnap kraken2 -y; then
    ok "Environment '$ENV_16S' created."
  else
    fail "Failed to create '$ENV_16S' — check conda logs above."
  fi
fi

log "[1/6] Done."
log ""

# ==============================================================================
# Step 2: Check per-user conda environments exist
# ==============================================================================
log "[2/6] Checking required conda environments..."

for env in "${REQUIRED_ENVS[@]}"; do
  if conda env list | awk '{print $1}' | grep -qx "$env"; then
    ok "$env"
  else
    fail "Environment '$env' not found."
    case "$env" in
      dragonflye)  fail "  Install: conda create -n dragonflye -c bioconda dragonflye -y" ;;
      seqkit)      fail "  Install: conda create -n seqkit -c bioconda seqkit -y" ;;
      bakta)       fail "  Install: conda create -n bakta -c bioconda bakta -y" ;;
      antismash8)  fail "  Install: pip install antismash  (see https://download.antismash.secondarymetabolites.org/)" ;;
    esac
    fail "  (this stage can be skipped with --skip-bakta or --skip-antismash)"
  fi
done

# Shared GTDBTk env
if [[ -d "$GTDBTK_ENV" ]]; then
  ok "Shared GTDBTk env: $GTDBTK_ENV"
else
  fail "Shared GTDBTk env not found at $GTDBTK_ENV"
  fail "  Ask your sysadmin. (Can be skipped at runtime with --skip-gtdbtk)"
fi

log "[2/6] Done."
log ""

# ==============================================================================
# Step 3: Check GNU parallel is available
# ==============================================================================
log "[3/6] Checking for GNU parallel..."

if command -v parallel &>/dev/null; then
  PAR_VER=$(parallel --version 2>/dev/null | head -1 || echo "unknown version")
  ok "parallel found: $PAR_VER"
else
  fail "GNU parallel not found in PATH."
  fail "  Install: conda install -c conda-forge parallel"
  fail "  Or: sudo apt-get install parallel   (if you have sudo access)"
  fail "  dragonSMASH uses parallel for Bakta and antiSMASH; it is required."
fi

log "[3/6] Done."
log ""

# ==============================================================================
# Step 4: Check shared databases
# ==============================================================================
log "[4/6] Checking shared databases..."

# GTDBTk DB
if [[ -d "$GTDBTK_DB" ]]; then
  ok "GTDBTk DB: $GTDBTK_DB"
else
  fail "GTDBTk DB not found at $GTDBTK_DB"
  fail "  Ask your sysadmin. (Can be skipped at runtime with --skip-gtdbtk)"
fi

# GTDBTk skani sketches
if [[ -d "$GTDBTK_SKANI" ]]; then
  ok "GTDBTk skani sketches: $GTDBTK_SKANI"
else
  fail "GTDBTk skani sketches not found at $GTDBTK_SKANI"
  fail "  Required for GTDBTk >= 2.3. Ask your sysadmin."
fi

# SILVA DB
SILVA_OK=true
for f in hash.k2d opts.k2d taxo.k2d; do
  if [[ ! -f "$SILVA_DB/$f" ]]; then
    fail "SILVA DB incomplete — missing $SILVA_DB/$f"
    SILVA_OK=false
  fi
done
if $SILVA_OK; then
  ok "SILVA kraken2 DB: $SILVA_DB"
fi

log "[4/6] Done."
log ""

# ==============================================================================
# Step 5: Verify shared GTDBTk environment
# ==============================================================================
log "[5/6] Verifying shared GTDBTk environment..."

if [[ -d "$GTDBTK_ENV" ]]; then
  GTDBTK_VER=$(conda run -p "$GTDBTK_ENV" --no-capture-output \
    gtdbtk --version 2>/dev/null | grep -oP '[\d]+\.[\d]+\.[\d]+' | head -1 || echo "unknown")
  log "      GTDBTk version: $GTDBTK_VER"

  MAJOR=$(echo "$GTDBTK_VER" | cut -d. -f1)
  MINOR=$(echo "$GTDBTK_VER" | cut -d. -f2)
  if [[ "$MAJOR" -ge 2 && "$MINOR" -ge 3 ]] 2>/dev/null; then
    ok "GTDBTk $GTDBTK_VER — skani supported ✓"
  else
    warn "GTDBTk $GTDBTK_VER is older than 2.3 — skani may not be available."
    warn "  dragonSMASH uses --skani_sketch_dir, which requires GTDBTk >= 2.3."
  fi

  if conda run -p "$GTDBTK_ENV" --no-capture-output which skani &>/dev/null; then
    ok "skani binary found in GTDBTk env."
  else
    warn "skani binary not found in GTDBTk env — ANI step may fall back to fastANI."
  fi
else
  warn "Skipping GTDBTk verification — env not found (see Step 2)."
fi

log "[5/6] Done."
log ""

# ==============================================================================
# Step 6: Check Bakta DB
# ==============================================================================
log "[6/6] Checking Bakta database..."

# Look for a Bakta DB in common locations
BAKTA_DB_FOUND=""
for candidate in \
    "$HOME/dbs/bakta/db" \
    "$HOME/dbs/bakta/v6/db" \
    "$HOME/bakta/db" \
    "$HOME/db/bakta/db"; do
  if [[ -d "$candidate" ]]; then
    BAKTA_DB_FOUND="$candidate"
    break
  fi
done

if [[ -n "$BAKTA_DB_FOUND" ]]; then
  ok "Bakta DB found: $BAKTA_DB_FOUND"
  BAKTA_DB_DEFAULT="$BAKTA_DB_FOUND"
else
  warn "No Bakta DB found in common locations."
  warn "  Download with:"
  warn "    conda run -n bakta bakta_db download --output ~/dbs/bakta/ --type full"
  warn "  Then set BAKTA_DB in your dragonsmash.conf."
  warn "  (Can be skipped at runtime with --skip-bakta if you don't need annotation)"
  BAKTA_DB_DEFAULT="$HOME/dbs/bakta/db"
fi

log "[6/6] Done."
log ""

# ==============================================================================
# Summary
# ==============================================================================
log "=============================="
log " Summary"
log "=============================="
log ""

PASS_COUNT=$(( 6 - ${#ERRORS[@]} ))  # rough estimate

if [[ ${#ERRORS[@]} -eq 0 ]] && [[ ${#WARNINGS[@]} -eq 0 ]]; then
  log "  All checks passed! dragonSMASH is ready to run."
elif [[ ${#ERRORS[@]} -eq 0 ]]; then
  log "  Setup complete with ${#WARNINGS[@]} warning(s) — review above."
else
  log "  Setup completed with ${#ERRORS[@]} error(s) — resolve before running dragonSMASH."
  log "  Stages with missing dependencies can be skipped at runtime:"
  log "    --skip-bakta       Skip Bakta annotation (and antiSMASH)"
  log "    --skip-antismash   Skip antiSMASH BGC detection only"
  log "    --skip-16s         Skip 16S rRNA species ID"
  log "    --skip-gtdbtk      Skip GTDBTk taxonomy"
fi

if [[ ${#WARNINGS[@]} -gt 0 ]]; then
  log ""
  log "  Warnings:"
  for w in "${WARNINGS[@]}"; do log "    ⚠  $w"; done
fi

if [[ ${#ERRORS[@]} -gt 0 ]]; then
  log ""
  log "  Errors:"
  for e in "${ERRORS[@]}"; do log "    ✗  $e"; done
fi

log ""
log "  ── Suggested dragonsmash.conf for $(whoami) ──────────────────────"
log ""
log "    BAKTA_DB=$BAKTA_DB_DEFAULT"
log "    ENV_NAME_DRAGONFLYE=dragonflye"
log "    ENV_NAME_SEQKIT=seqkit"
log "    ENV_NAME_BAKTA=bakta"
log "    ENV_NAME_ANTISMASH=antismash8"
log "    GTDBTK_CONDA_PREFIX=$GTDBTK_ENV"
log "    GTDBTK_DATA_PATH=$GTDBTK_DB"
log "    GTDBTK_SKANI_SKETCH_DIR=$GTDBTK_SKANI"
log "    ENV_NAME_BARRNAP=$ENV_16S"
log "    SILVA_DB=$SILVA_DB"
log "    ANTISMASH_TAXON=bacteria"
log "    ANTISMASH_EXTRA_ARGS=--cc-mibig --cb-general --cb-subcluster --cb-knownclusters --rre"
log "    BARRNAP_KINGDOM=bac"
log "    BARRNAP_MIN_SCORE=0.8"
log ""
log "  Copy the above into a file called 'dragonsmash.conf' in the same"
log "  directory as dragonSMASHv3.bash, then edit BAKTA_DB to your actual path."
log ""
log "  Log saved to: $LOG_FILE"
log ""
