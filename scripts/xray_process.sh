#!/usr/bin/env bash
# Stage 7: reprocess the X-ray observations of CGCG 312-012 and extract spectra.
#
#     bash scripts/xray_process.sh xrt10|xrt22|nustar      (pipelines, one at a time)
#     bash scripts/xray_process.sh regions                 (scripts/xray_regions.py)
#     bash scripts/xray_process.sh products                (spectra, responses, grouping)
#     bash scripts/xray_process.sh products_saa_none       (stage 7b: NuSTAR spectra of the default SAA run)
#
# Settings: docs/XRAY_PLAN.md section 3. Input: data/cache/s7/{swift,nustar,bat157}
# (scripts/fetch_xray.py). Output: data/cache/s7/work/ (not in git).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
set +u
source "$ROOT/scripts/heasoft_env.sh" > /dev/null
set -u
S7="$ROOT/data/cache/s7"
WORK="$S7/work"
RA=142.4078319
DEC=62.5440626
NUOBS=60160365002

xrt_pipeline() {  # $1 = xrt10|xrt22, $2 = obsid
  mkdir -p "$WORK/$1"
  rm -rf "$WORK/$1/pipe"
  xrtpipeline indir="$S7/swift/$2" steminputs="sw$2" outdir="$WORK/$1/pipe" \
    srcra=$RA srcdec=$DEC datamode=PC createexpomap=yes cleanup=no clobber=yes \
    > "$WORK/$1/xrtpipeline.log" 2>&1
  tail -3 "$WORK/$1/xrtpipeline.log"
}

xrt_products() {  # $1 = xrt10|xrt22, $2 = obsid
  local d="$WORK/$1" evt="sw$2xpcw3po_cl.evt" expo="sw$2xpcw3po_ex.img"
  mkdir -p "$d/prod"
  rm -f "$d/prod/"*
  ( cd "$d/prod"
    cat > extract.xco <<EOF
xsel_$1
read events $evt
../pipe
yes
filter region ../src_$1.reg
extract spectrum
save spectrum src.pha
clear region
filter region ../bkg_$1.reg
extract spectrum
save spectrum bkg.pha
exit
no
EOF
    xselect @extract.xco > xselect.log 2>&1
    xrtmkarf phafile=src.pha srcx=-1 srcy=-1 outfile=src.arf psfflag=yes \
      expofile="../pipe/$expo" clobber=yes > xrtmkarf.log 2>&1
    rmf="$(grep -o "[^ ']*swxpc0to12s6[^ ']*\.rmf" xrtmkarf.log | head -1)"
    cp "$rmf" src.rmf
    basename "$rmf" > rmf_name.txt
    ftgrouppha infile=src.pha outfile=src_min1.pha grouptype=min groupscale=1 clobber=yes > group.log 2>&1
  )
  ls -l "$d/prod" | awk '{print $5, $9}'
}

case "${1:-}" in
  xrt10) xrt_pipeline xrt10 00040946001 ;;
  xrt22) xrt_pipeline xrt22 00081051001 ;;
  nustar)
    mkdir -p "$WORK/nustar"
    rm -rf "$WORK/nustar/pipe"
    nupipeline indir="$S7/nustar/$NUOBS" steminputs="nu$NUOBS" outdir="$WORK/nustar/pipe" \
      ${NU_EXTRA:-} clobber=yes > "$WORK/nustar/nupipeline.log" 2>&1
    tail -3 "$WORK/nustar/nupipeline.log"
    ;;
  regions) python "$ROOT/scripts/xray_regions.py" ;;
  products)
    if [ "${SKIP_XRT:-0}" != 1 ]; then
      xrt_products xrt10 00040946001
      xrt_products xrt22 00081051001
    fi
    for m in A B; do
      k="fpm$(echo $m | tr AB ab)"
      rm -rf "$WORK/nustar/prod_$k"
      # relative paths: nuexpomap refuses long file names
      ( cd "$WORK/nustar"
        nuproducts indir=pipe instrument="FPM$m" steminputs="nu$NUOBS" \
          outdir="prod_$k" srcregionfile="src_$k.reg" \
          bkgregionfile="bkg_$k.reg" bkgextract=yes rungrppha=no clobber=yes \
          > "nuproducts_$k.log" 2>&1 )
      ( cd "$WORK/nustar/prod_$k"
        ftgrouppha infile="nu${NUOBS}${m}01_sr.pha" outfile=src_min1.pha grouptype=min groupscale=1 clobber=yes
        ftgrouppha infile="nu${NUOBS}${m}01_sr.pha" outfile=src_min20.pha grouptype=min groupscale=20 clobber=yes )
      ls -l "$WORK/nustar/prod_$k" | awk '{print $5, $9}'
    done
    mkdir -p "$WORK/bat"
    cp "$S7/bat157/SWIFT_J0929.7+6232.pha" "$WORK/bat/bat157.pha"
    cp "$S7/bat157/BAT_157m_swiftbat_survey_full_157m.rsp" "$WORK/bat/bat157.rsp"
    nh equinox=2000 ra=$RA dec=$DEC > "$WORK/nh.log" 2>&1
    tail -4 "$WORK/nh.log"
    ;;
  products_saa_none)
    # stage 7b: spectra from the event files of the default SAA setting (same regions as the main run).
    # The first nupipeline run has to be kept as work/nustar/pipe_saamode_none.
    for m in A B; do
      k="fpm$(echo $m | tr AB ab)"
      rm -rf "$WORK/nustar/prod_none_$k"
      ( cd "$WORK/nustar"
        nuproducts indir=pipe_saamode_none instrument="FPM$m" steminputs="nu$NUOBS" \
          outdir="prod_none_$k" srcregionfile="src_$k.reg" \
          bkgregionfile="bkg_$k.reg" bkgextract=yes rungrppha=no clobber=yes \
          > "nuproducts_none_$k.log" 2>&1 )
      ( cd "$WORK/nustar/prod_none_$k"
        ftgrouppha infile="nu${NUOBS}${m}01_sr.pha" outfile=src_min1.pha grouptype=min groupscale=1 clobber=yes
        rm -f ./*_skymap*.img ./*_det2map.img )
      ls -l "$WORK/nustar/prod_none_$k" | awk '{print $5, $9}'
    done
    ;;
  *) echo "usage: $0 xrt10|xrt22|nustar|regions|products|products_saa_none"; exit 1 ;;
esac
