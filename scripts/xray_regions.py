"""Stage 7: source and background regions from the reprocessed event files.

    python scripts/xray_regions.py          (inside the HEASoft environment)

Rules: docs/XRAY_PLAN.md section 3. Reads the cleaned event files made by
scripts/xray_process.sh (data/cache/s7/work/), writes ds9 region files (fk5)
next to them, results/s7_regions.json (positions, counts, background light
curve check) and results/fig_s7_xray_images.png.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from irecho import config, xray  # noqa: E402

WORK = xray.S7_DIR / "work"
RA, DEC = xray.TARGET["ra_deg"], xray.TARGET["dec_deg"]
NU_SRC_R, NU_BKG_R = 50.0, 100.0            # arcsec
XRT_SRC_PIX, XRT_BKG_IN, XRT_BKG_OUT = 20.0, 60.0, 110.0  # pixels (2.357 arcsec)

EVENTS = {
    "fpma": WORK / "nustar" / "pipe" / "nu60160365002A01_cl.evt",
    "fpmb": WORK / "nustar" / "pipe" / "nu60160365002B01_cl.evt",
    "xrt10": WORK / "xrt10" / "pipe" / "sw00040946001xpcw3po_cl.evt",
    "xrt22": WORK / "xrt22" / "pipe" / "sw00081051001xpcw3po_cl.evt",
}


def load(path: Path) -> dict:
    with fits.open(path) as h:
        ev = h["EVENTS"]
        d, hd = ev.data, ev.header
        names = [c.name for c in ev.columns]
        ix, iy = names.index("X") + 1, names.index("Y") + 1
        w = WCS(naxis=2)
        w.wcs.crpix = [hd[f"TCRPX{ix}"], hd[f"TCRPX{iy}"]]
        w.wcs.crval = [hd[f"TCRVL{ix}"], hd[f"TCRVL{iy}"]]
        w.wcs.cdelt = [hd[f"TCDLT{ix}"], hd[f"TCDLT{iy}"]]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        out = {"x": np.array(d["X"], dtype=float), "y": np.array(d["Y"], dtype=float),
               "pi": np.array(d["PI"], dtype=float), "time": np.array(d["TIME"], dtype=float),
               "wcs": w, "scale": abs(hd[f"TCDLT{iy}"]) * 3600.0,
               "exposure": float(hd.get("EXPOSURE", np.nan)),
               "header": {k: str(hd.get(k)) for k in ("DATE-OBS", "DATE-END", "SOFTVER", "CALDBVER", "NUPRODS",
                                                      "XRTVSUB", "OBS_ID", "DATAMODE", "ONTIME", "LIVETIME")
                          if hd.get(k) is not None}}
        if "DET_ID" in names:
            out["det_id"] = np.array(d["DET_ID"], dtype=int)
        good = (out["x"] > 0) & (out["y"] > 0)
        for k in ("x", "y", "pi", "time", "det_id"):
            if k in out:
                out[k] = out[k][good]
    return out


def centroid(ev: dict, x0: float, y0: float, radius_arcsec: float = 30.0) -> tuple[float, float]:
    cx, cy = x0, y0
    for _ in range(3):
        m = np.hypot(ev["x"] - cx, ev["y"] - cy) * ev["scale"] < radius_arcsec
        if m.sum() < 5:
            break
        cx, cy = float(ev["x"][m].mean()), float(ev["y"][m].mean())
    return cx, cy


def sky(ev: dict, x: float, y: float) -> tuple[float, float]:
    ra, dec = ev["wcs"].wcs_pix2world([[x, y]], 1)[0]
    return float(ra), float(dec)


def nustar_background(ev: dict, cx: float, cy: float) -> dict:
    """Centre of a source-free circle, on the detector chip of the source if possible."""
    scale = ev["scale"]
    r = np.hypot(ev["x"] - cx, ev["y"] - cy) * scale
    src_det = int(np.bincount(ev["det_id"][r < NU_SRC_R]).argmax())
    step = 10.0  # pixels
    xe = np.arange(ev["x"].min(), ev["x"].max() + step, step)
    ye = np.arange(ev["y"].min(), ev["y"].max() + step, step)
    img, _, _ = np.histogram2d(ev["x"], ev["y"], bins=[xe, ye])
    fov = img >= 0.2 * np.median(img[img > 0])
    xc, yc = 0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:])
    gx, gy = np.meshgrid(xc, yc, indexing="ij")
    for radius, min_same in ((NU_BKG_R, 0.95), (80.0, 0.95), (NU_BKG_R, 0.0), (80.0, 0.0)):
        best = None
        for i in range(len(xc)):
            for j in range(len(yc)):
                sep = np.hypot(xc[i] - cx, yc[j] - cy) * scale
                if sep < NU_SRC_R + radius + 40.0:
                    continue
                inside = np.hypot(gx - xc[i], gy - yc[j]) * scale < radius + 0.7 * step * scale
                if not fov[inside].all():
                    continue
                m = np.hypot(ev["x"] - xc[i], ev["y"] - yc[j]) * scale < radius
                same = float(np.mean(ev["det_id"][m] == src_det)) if m.any() else 0.0
                if same < min_same:
                    continue
                key = (round(same, 2), sep)
                if best is None or key > best[0]:
                    best = (key, xc[i], yc[j], int(m.sum()), same, sep)
        if best is not None:
            return {"x": float(best[1]), "y": float(best[2]), "radius_arcsec": radius, "counts": best[3],
                    "same_chip_fraction": best[4], "separation_arcsec": float(best[5]), "source_det_id": src_det}
    raise SystemExit("no background circle found")


def xrt_contaminants(ev: dict, cx: float, cy: float) -> list[dict]:
    """Bright spots inside the background annulus (sliding 10-pixel cell)."""
    from scipy.stats import poisson

    r = np.hypot(ev["x"] - cx, ev["y"] - cy)
    ann = (r > XRT_BKG_IN) & (r < XRT_BKG_OUT)
    area = np.pi * (XRT_BKG_OUT ** 2 - XRT_BKG_IN ** 2)
    cell = 10.0
    expect = ann.sum() / area * np.pi * cell ** 2
    limit = poisson.isf(1e-5, expect)
    found = []
    x, y = ev["x"][ann], ev["y"][ann]
    for gx in np.arange(cx - XRT_BKG_OUT, cx + XRT_BKG_OUT + 1, 5.0):
        for gy in np.arange(cy - XRT_BKG_OUT, cy + XRT_BKG_OUT + 1, 5.0):
            if not XRT_BKG_IN + cell <= np.hypot(gx - cx, gy - cy) <= XRT_BKG_OUT - cell:
                continue
            n = int((np.hypot(x - gx, y - gy) < cell).sum())
            if n > limit and not any(np.hypot(gx - f["x"], gy - f["y"]) < 25 for f in found):
                found.append({"x": float(gx), "y": float(gy), "counts": n, "expected": float(expect)})
    return found


def background_lightcurve(ev: dict, cx: float, cy: float, binsize: float = 200.0) -> dict:
    """Count rate outside 120 arcsec of the source in time bins (looks for SAA-related rises)."""
    m = np.hypot(ev["x"] - cx, ev["y"] - cy) * ev["scale"] > 120.0
    t = np.sort(ev["time"][m])
    gaps = np.where(np.diff(t) > 300.0)[0]
    starts = np.concatenate([[t[0]], t[gaps + 1]])
    stops = np.concatenate([t[gaps], [t[-1]]])
    rates, errs = [], []
    for a, b in zip(starts, stops):
        nb = int((b - a) // binsize)
        for k in range(nb):
            n = int(((t >= a + k * binsize) & (t < a + (k + 1) * binsize)).sum())
            rates.append(n / binsize)
            errs.append(np.sqrt(max(n, 1)) / binsize)
    rates, errs = np.array(rates), np.array(errs)
    med = float(np.median(rates))
    dev = (rates - med) / np.sqrt(med / binsize)
    return {"binsize_s": binsize, "n_bins": int(rates.size), "n_orbits": int(starts.size), "median_rate": med,
            "max_rate": float(rates.max()), "max_deviation_sigma": float(dev.max()),
            "n_bins_above_5_sigma": int((dev > 5).sum()), "chi2_constant": float(np.sum(((rates - rates.mean()) / errs) ** 2))}


def main() -> None:
    out, panels = {}, []
    for key, path in EVENTS.items():
        ev = load(path)
        x0, y0 = ev["wcs"].wcs_world2pix([[RA, DEC]], 1)[0]
        cx, cy = centroid(ev, x0, y0)
        offset = float(np.hypot(cx - x0, cy - y0) * ev["scale"])
        ra, dec = sky(ev, cx, cy)
        rec = {"event_file": str(path.relative_to(config.PROJECT_ROOT)), "exposure_s": ev["exposure"],
               "n_events": int(ev["x"].size), "pixel_scale_arcsec": ev["scale"], "header": ev["header"],
               "target_xy": [float(x0), float(y0)], "centroid_xy": [cx, cy],
               "centroid_ra_dec": [ra, dec], "centroid_offset_arcsec": offset}
        r_as = np.hypot(ev["x"] - cx, ev["y"] - cy) * ev["scale"]
        regdir = path.parent.parent
        if key.startswith("fpm"):
            bkg = nustar_background(ev, cx, cy)
            bra, bdec = sky(ev, bkg["x"], bkg["y"])
            rec.update(source_radius_arcsec=NU_SRC_R, source_counts=int((r_as < NU_SRC_R).sum()),
                       background=dict(bkg, ra=bra, dec=bdec),
                       background_lightcurve=background_lightcurve(ev, cx, cy))
            rec["expected_background_in_source"] = bkg["counts"] * (NU_SRC_R / bkg["radius_arcsec"]) ** 2
            (regdir / f"src_{key}.reg").write_text(f"fk5\ncircle({ra:.6f},{dec:.6f},{NU_SRC_R}\")\n")
            (regdir / f"bkg_{key}.reg").write_text(f"fk5\ncircle({bra:.6f},{bdec:.6f},{bkg['radius_arcsec']}\")\n")
            panels.append((key, ev, cx, cy, [("c", cx, cy, NU_SRC_R / ev["scale"]),
                                              ("c", bkg["x"], bkg["y"], bkg["radius_arcsec"] / ev["scale"])]))
        else:
            r_pix = np.hypot(ev["x"] - cx, ev["y"] - cy)
            cont = xrt_contaminants(ev, cx, cy)
            ann = (r_pix > XRT_BKG_IN) & (r_pix < XRT_BKG_OUT)
            for c in cont:
                ann &= np.hypot(ev["x"] - c["x"], ev["y"] - c["y"]) >= 25.0
            sc = ev["scale"]
            lines = [f"fk5\nannulus({ra:.6f},{dec:.6f},{XRT_BKG_IN * sc:.2f}\",{XRT_BKG_OUT * sc:.2f}\")"]
            for c in cont:
                c["ra"], c["dec"] = sky(ev, c["x"], c["y"])
                lines.append(f"-circle({c['ra']:.6f},{c['dec']:.6f},{25.0 * sc:.2f}\")")
            (regdir / f"src_{key}.reg").write_text(f"fk5\ncircle({ra:.6f},{dec:.6f},{XRT_SRC_PIX * sc:.2f}\")\n")
            (regdir / f"bkg_{key}.reg").write_text("\n".join(lines) + "\n")
            src = r_pix < XRT_SRC_PIX
            e_kev = ev["pi"] * 0.01
            rec.update(source_radius_pix=XRT_SRC_PIX, source_counts=int(src.sum()),
                       source_counts_0p3_2=int((src & (e_kev >= 0.3) & (e_kev < 2)).sum()),
                       source_counts_2_10=int((src & (e_kev >= 2) & (e_kev <= 10)).sum()),
                       background_annulus_pix=[XRT_BKG_IN, XRT_BKG_OUT], background_counts=int(ann.sum()),
                       background_counts_0p3_2=int((ann & (e_kev >= 0.3) & (e_kev < 2)).sum()),
                       background_counts_2_10=int((ann & (e_kev >= 2) & (e_kev <= 10)).sum()),
                       contaminants_masked=cont)
            shapes = [("c", cx, cy, XRT_SRC_PIX), ("c", cx, cy, XRT_BKG_IN), ("c", cx, cy, XRT_BKG_OUT)]
            shapes += [("x", c["x"], c["y"], 25.0) for c in cont]
            panels.append((key, ev, cx, cy, shapes))
        out[key] = rec
        print(key, json.dumps({k: v for k, v in rec.items() if k not in ("header",)}, ensure_ascii=False)[:900], flush=True)
    (config.RESULTS_DIR / "s7_regions.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    titles = {"fpma": "NuSTAR FPMA, 2022-05-01 (3-79 keV)", "fpmb": "NuSTAR FPMB, 2022-05-01 (3-79 keV)",
              "xrt10": "Swift XRT, 2010-12-29 (0.3-10 keV)", "xrt22": "Swift XRT, 2022-05-01 (0.3-10 keV)"}
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.2))
    for ax, (key, ev, cx, cy, shapes) in zip(axes, panels):
        half = 150 if key.startswith("fpm") else 130
        binw = 3 if key.startswith("fpm") else 4
        edges_x = np.arange(cx - half, cx + half + binw, binw)
        edges_y = np.arange(cy - half, cy + half + binw, binw)
        # energy band of the title: NuSTAR PI = (E - 1.6) / 0.04, XRT PI = E / 0.01
        band = (ev["pi"] >= 35) & (ev["pi"] <= 1935) if key.startswith("fpm") else (ev["pi"] >= 30) & (ev["pi"] <= 1000)
        img, _, _ = np.histogram2d(ev["x"][band], ev["y"][band], bins=[edges_x, edges_y])
        ax.imshow(np.sqrt(img.T), origin="lower", cmap="gray_r",
                  extent=[edges_x[0], edges_x[-1], edges_y[0], edges_y[-1]], interpolation="nearest")
        for kind, x, y, rad in shapes:
            ax.add_patch(plt.Circle((x, y), rad, fill=False, lw=1.0, color="#eb6834" if kind == "c" else "#2a78d6",
                                    ls="-" if kind == "c" else "--"))
        ax.set_title(titles[key], fontsize=9)
        ax.set_xlabel("sky X (pixel)", fontsize=8)
        ax.tick_params(labelsize=7)
    axes[0].set_ylabel("sky Y (pixel)", fontsize=8)
    fig.suptitle("CGCG 312-012: event images (square-root scale) with source and background regions", fontsize=10)
    fig.tight_layout()
    fig.savefig(config.RESULTS_DIR / "fig_s7_xray_images.png", dpi=130)


if __name__ == "__main__":
    main()
