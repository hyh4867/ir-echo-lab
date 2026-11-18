# Stage 7: set up the HEASoft environment (conda prefix .heasoft/ in the project
# folder) and the local CALDB (data/cache/s7/caldb). Source this file:
#     source scripts/heasoft_env.sh
# Install: docs/PIPELINE.md section 15.
_here="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
_conda="${CONDA_EXE:-$(command -v conda)}"
eval "$("$_conda" shell.bash hook)"
conda activate "$_here/.heasoft"
# The CALDB software fails when the path to a calibration file is long
# ("cifsl3: Error constructing path"), so CALDB points at a short symbolic
# link in the user's temporary folder. Only the link is outside the project.
_link="${TMPDIR:-/tmp}"
_link="${_link%/}/irecho-s7-caldb"
ln -sfn "$_here/data/cache/s7/caldb" "$_link"
export CALDB="$_link"
export CALDBCONFIG="$CALDB/caldb.config"
export CALDBALIAS="$CALDB/alias_config.fits"
export HEADASNOQUERY=1
export HEADASPROMPT=/dev/null
export PFILES="$_here/data/cache/s7/pfiles;$HEADAS/syspfiles"
mkdir -p "$_here/data/cache/s7/pfiles"
