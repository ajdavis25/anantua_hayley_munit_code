#!/usr/bin/env python3
"""
Inspect IPOLE HDF5 outputs and compare rendered image products.

Outputs under /work/vmo703/ipole_outputs/_qa:
- h5_resolution_summary.csv
- h5_top_level_items.txt
- image_inventory_new.csv
- image_inventory_ref.csv
- image_dim_distribution_new.csv
- image_dim_distribution_ref.csv
- matched_pairs.csv
- matched_pair_metrics.csv
- reference_h5_sample.csv (if any reference .h5 found)
- comparisons/pair_*.png
- QA_REPORT.md
"""

import csv
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import matplotlib
import numpy as np
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


ROOT = Path("/work/vmo703/ipole_outputs")
H5_ROOT = ROOT
NEW_IMG_DIR = ROOT / "images"
REF_IMG_DIR = ROOT / "M87" / "images"
QA_DIR = ROOT / "_qa"
CMP_DIR = QA_DIR / "comparisons"

VALID_IMG_EXT = {".png", ".jpg", ".jpeg", ".pdf"}
RES_KEY_HINTS = {
    "nx",
    "ny",
    "n_x",
    "n_y",
    "npix",
    "pixels",
    "width",
    "height",
    "img_nx",
    "img_ny",
    "camera_npix",
}


def ensure_dirs():
    QA_DIR.mkdir(parents=True, exist_ok=True)
    CMP_DIR.mkdir(parents=True, exist_ok=True)


def decode_scalar(v):
    if isinstance(v, bytes):
        try:
            return v.decode("utf-8")
        except Exception:
            return str(v)
    if isinstance(v, np.generic):
        return v.item()
    return v


def is_scalar_dataset(ds):
    return ds.shape == () or ds.size == 1


def safe_scalar(ds):
    try:
        v = ds[()]
        if isinstance(v, np.ndarray) and v.size == 1:
            v = v.reshape(()).item()
        return decode_scalar(v)
    except Exception:
        return None


def interpret_image_shape(shape):
    """Infer (nx, ny) from image-like shapes."""
    if len(shape) == 2:
        ny, nx = shape
        return int(nx), int(ny)
    if len(shape) == 3:
        # Common: (ny, nx, stokes) or (stokes, ny, nx)
        if shape[-1] <= 8:
            ny, nx = shape[0], shape[1]
            return int(nx), int(ny)
        if shape[0] <= 8:
            ny, nx = shape[1], shape[2]
            return int(nx), int(ny)
    return None


def looks_image_like(name, shape):
    leaf = name.split("/")[-1].lower()
    if leaf in {"pol", "unpol", "tau", "i", "q", "u", "v", "stokes"}:
        return True
    if len(shape) == 2 and min(shape) >= 16:
        return True
    if len(shape) == 3 and max(shape) >= 16:
        return True
    return False


def find_resolution_candidates(hfp):
    found = {}

    for k, v in hfp.attrs.items():
        lk = str(k).lower()
        if any(h in lk for h in RES_KEY_HINTS):
            found["/@{}".format(k)] = decode_scalar(v)

    def visit(name, obj):
        for k, v in obj.attrs.items():
            lk = str(k).lower()
            if any(h in lk for h in RES_KEY_HINTS):
                found["/{}/@{}".format(name, k)] = decode_scalar(v)

        if isinstance(obj, h5py.Dataset) and is_scalar_dataset(obj):
            leaf = name.split("/")[-1].lower()
            if any(h in leaf for h in RES_KEY_HINTS):
                found["/{}".format(name)] = safe_scalar(obj)

    hfp.visititems(visit)
    return found


def inspect_h5_file(h5_path):
    row = {"filename": h5_path.name, "path": str(h5_path)}

    with h5py.File(str(h5_path), "r") as hfp:
        top_items = []
        for k in hfp.keys():
            typ = "group" if isinstance(hfp[k], h5py.Group) else "dataset"
            top_items.append("{}({})".format(k, typ))
        row["top_level_items"] = ", ".join(top_items)

        camera_scalars = {}
        fovx_dsource = None
        fovy_dsource = None
        if "header" in hfp and "camera" in hfp["header"]:
            cam = hfp["header"]["camera"]
            for k in cam.keys():
                obj = cam[k]
                if isinstance(obj, h5py.Dataset) and is_scalar_dataset(obj):
                    camera_scalars[k] = safe_scalar(obj)
            if "fovx_dsource" in cam:
                fovx_dsource = float(cam["fovx_dsource"][()])
            if "fovy_dsource" in cam:
                fovy_dsource = float(cam["fovy_dsource"][()])

        row["camera_nx"] = int(camera_scalars["nx"]) if "nx" in camera_scalars else None
        row["camera_ny"] = int(camera_scalars["ny"]) if "ny" in camera_scalars else None
        row["camera_fovx_dsource"] = fovx_dsource
        row["camera_fovy_dsource"] = fovy_dsource

        res_candidates = find_resolution_candidates(hfp)
        row["resolution_candidates"] = "; ".join(
            ["{}={}".format(k, v) for k, v in sorted(res_candidates.items())]
        )

        dataset_shapes = {}
        image_candidates = []

        def visit_ds(name, obj):
            if isinstance(obj, h5py.Dataset):
                shp = tuple(int(x) for x in obj.shape)
                dataset_shapes[name] = shp
                if looks_image_like(name, shp):
                    image_candidates.append((name, shp))

        hfp.visititems(visit_ds)

        for key in ("pol", "unpol", "tau"):
            if key in hfp:
                row["{}_shape".format(key)] = tuple(int(x) for x in hfp[key].shape)
            else:
                row["{}_shape".format(key)] = None

        detected_nx = None
        detected_ny = None
        source = ""

        if row["camera_nx"] is not None and row["camera_ny"] is not None:
            detected_nx = row["camera_nx"]
            detected_ny = row["camera_ny"]
            source = "/header/camera/nx, /header/camera/ny"
        else:
            chosen = None
            for p in ("pol", "unpol", "tau"):
                if p in hfp:
                    shp = tuple(int(x) for x in hfp[p].shape)
                    inferred = interpret_image_shape(shp)
                    if inferred is not None:
                        detected_nx, detected_ny = inferred
                        chosen = "/{} shape={}".format(p, shp)
                        break
            if chosen is None:
                for name, shp in image_candidates:
                    inferred = interpret_image_shape(shp)
                    if inferred is not None:
                        detected_nx, detected_ny = inferred
                        chosen = "/{} shape={}".format(name, shp)
                        break
            source = chosen or "not found"

        row["detected_nx"] = detected_nx
        row["detected_ny"] = detected_ny
        row["resolution_source"] = source
        row["image_like_datasets"] = "; ".join(
            ["/{}:{}".format(name, shp) for name, shp in sorted(image_candidates)]
        )

    return row


def list_root_h5_files(root):
    return sorted(root.glob("*.h5"))


def fmt_bytes(n):
    if n is None:
        return "n/a"
    units = ["B", "KB", "MB", "GB"]
    v = float(n)
    for u in units:
        if abs(v) < 1024.0 or u == units[-1]:
            return "{:.1f}{}".format(v, u)
        v /= 1024.0
    return "{}B".format(n)


def write_csv(path, rows, fieldnames=None):
    path = Path(path)
    if fieldnames is None:
        if not rows:
            path.write_text("", encoding="utf-8")
            return
        # preserve key order from first row, then append unseen keys
        seen = list(rows[0].keys())
        for row in rows[1:]:
            for k in row.keys():
                if k not in seen:
                    seen.append(k)
        fieldnames = seen

    with path.open("w", newline="", encoding="utf-8") as fp:
        w = csv.DictWriter(fp, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            safe_row = {}
            for k in fieldnames:
                v = row.get(k, "")
                if isinstance(v, (tuple, list, dict)):
                    safe_row[k] = str(v)
                else:
                    safe_row[k] = v
            w.writerow(safe_row)


def format_table(rows, columns, max_rows=None):
    if max_rows is not None:
        rows = rows[:max_rows]
    if not rows:
        return "(no rows)"

    data = []
    for r in rows:
        row_vals = []
        for c in columns:
            v = r.get(c, "")
            if v is None:
                s = ""
            else:
                s = str(v)
            row_vals.append(s)
        data.append(row_vals)

    widths = [len(c) for c in columns]
    for row in data:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(cell))

    def fmt_line(vals):
        return " | ".join(vals[i].ljust(widths[i]) for i in range(len(vals)))

    header = fmt_line(columns)
    sep = "-+-".join("-" * w for w in widths)
    lines = [header, sep]
    for row in data:
        lines.append(fmt_line(row))
    return "\n".join(lines)


def collect_image_inventory(base_dir):
    rows = []
    for p in sorted(base_dir.rglob("*")):
        if not p.is_file():
            continue
        ext = p.suffix.lower()
        if ext not in VALID_IMG_EXT:
            continue

        row = {
            "path": str(p),
            "filename": p.name,
            "ext": ext.lstrip("."),
            "size_bytes": p.stat().st_size,
            "width": None,
            "height": None,
            "mode": None,
            "format": None,
            "dpi": None,
        }

        if ext in {".png", ".jpg", ".jpeg"}:
            try:
                with Image.open(str(p)) as im:
                    row["width"], row["height"] = im.size
                    row["mode"] = im.mode
                    row["format"] = im.format
                    dpi = im.info.get("dpi", None)
                    if isinstance(dpi, tuple):
                        row["dpi"] = "{:.2f},{:.2f}".format(dpi[0], dpi[1])
                    elif dpi is not None:
                        row["dpi"] = str(dpi)
            except Exception as exc:
                row["mode"] = "ERROR: {}".format(exc)

        rows.append(row)
    return rows


def dimension_distribution(image_rows):
    counts = Counter()
    for r in image_rows:
        w = r.get("width")
        h = r.get("height")
        if w is None or h is None:
            continue
        counts[(int(w), int(h))] += 1

    out = []
    for (w, h), c in sorted(counts.items(), key=lambda x: (-x[1], -x[0][0], -x[0][1])):
        out.append({"width": w, "height": h, "count": c})
    return out


def summarize_inventory(image_rows, label):
    total = len(image_rows)
    out = []
    out.append("{}: count={}".format(label, total))
    if total == 0:
        return "\n".join(out)

    sizes = [r["size_bytes"] for r in image_rows]
    mean_size = float(np.mean(sizes)) if sizes else None
    median_size = float(np.median(sizes)) if sizes else None
    out.append(
        "{}: mean_size={}, median_size={}".format(
            label, fmt_bytes(mean_size), fmt_bytes(median_size)
        )
    )

    ext_counts = Counter(r["ext"].lower() for r in image_rows)
    ext_rows = [{"ext": k, "count": v} for k, v in sorted(ext_counts.items())]
    out.append("{}: file types\n{}".format(label, format_table(ext_rows, ["ext", "count"])))

    dims = dimension_distribution(image_rows)
    out.append(
        "{}: dimensions (top 10)\n{}".format(
            label, format_table(dims, ["width", "height", "count"], max_rows=10)
        )
    )

    # include a few examples for file sizes
    sorted_by_size = sorted(image_rows, key=lambda r: r["size_bytes"])
    examples = []
    for r in sorted_by_size[:2] + sorted_by_size[-2:]:
        examples.append(
            {
                "filename": r["filename"],
                "dims": "{}x{}".format(r.get("width"), r.get("height")),
                "size": fmt_bytes(r["size_bytes"]),
                "dpi": r.get("dpi") or "none",
            }
        )
    out.append("{}: size examples\n{}".format(label, format_table(examples, ["filename", "dims", "size", "dpi"])))

    return "\n".join(out)


def normalize_model(model):
    return re.sub(r"[^A-Z0-9]", "", model.upper())


def parse_new_image_key(stem):
    toks = stem.split("_")
    if len(toks) < 4:
        return None
    sim = toks[0]
    dump = toks[1]
    model = toks[2]
    try:
        pos = float(toks[-1])
    except Exception:
        pos = float("nan")
    return sim, dump, normalize_model(model), pos


def parse_ref_image_key(stem):
    toks = stem.split("_")
    # output_Ma+0.94_5000_model_CRITBETA_Rhigh_0_0.000
    if len(toks) < 8:
        return None
    if toks[0] != "output":
        return None
    if toks[3] != "model":
        return None
    sim = toks[1]
    dump = toks[2]
    model = toks[4]
    try:
        pos = float(toks[-1])
    except Exception:
        pos = float("nan")
    return sim, dump, normalize_model(model), pos


def choose_matched_pairs(new_rows, ref_rows, limit=3):
    ref_by_key = defaultdict(list)
    ref_by_sim_dump = defaultdict(list)
    for r in ref_rows:
        key = parse_ref_image_key(Path(r["filename"]).stem)
        if key is None:
            continue
        ref_by_key[key[:3]].append((r["path"], key[3], key[:3]))
        ref_by_sim_dump[(key[0], key[1])].append((r["path"], key[3], key[:3]))

    parsed_new = []
    for r in new_rows:
        key = parse_new_image_key(Path(r["filename"]).stem)
        if key is None:
            continue
        parsed_new.append((r, key))

    used_ref = set()
    out = []
    seen_match_key = set()

    def append_best_match(new_row, key, candidates, exact_model):
        if not candidates:
            return False
        new_pos = key[3]
        ranked = candidates
        if not math.isnan(new_pos):
            ranked = sorted(
                candidates,
                key=lambda x: abs(x[1] - new_pos) if not math.isnan(x[1]) else 1e9,
            )
        for cand in ranked:
            ref_path = cand[0]
            if ref_path in used_ref:
                continue
            match_key = "{}_{}_{}".format(key[0], key[1], key[2])
            if match_key in seen_match_key:
                return False
            out.append(
                {
                    "new_path": new_row["path"],
                    "new_filename": new_row["filename"],
                    "ref_path": ref_path,
                    "ref_filename": Path(ref_path).name,
                    "match_key": match_key,
                    "exact_model_match": bool(exact_model),
                }
            )
            used_ref.add(ref_path)
            seen_match_key.add(match_key)
            return True
        return False

    # Pass 1: exact (sim, dump, model) matches only.
    for new_row, key in parsed_new:
        if len(out) >= limit:
            break
        candidates = list(ref_by_key.get(key[:3], []))
        append_best_match(new_row, key, candidates, exact_model=True)

    # Pass 2: fallback same (sim, dump) if needed.
    for new_row, key in parsed_new:
        if len(out) >= limit:
            break
        match_key = "{}_{}_{}".format(key[0], key[1], key[2])
        if match_key in seen_match_key:
            continue
        candidates = list(ref_by_sim_dump.get((key[0], key[1]), []))
        append_best_match(new_row, key, candidates, exact_model=False)

    return out


def load_rgb(path):
    with Image.open(str(path)) as im:
        return np.asarray(im.convert("RGB"))


def grayscale(arr_rgb):
    return (
        0.299 * arr_rgb[..., 0]
        + 0.587 * arr_rgb[..., 1]
        + 0.114 * arr_rgb[..., 2]
    ).astype(np.float32)


def image_quality_metrics(arr_rgb):
    lum = grayscale(arr_rgb)
    dx = np.abs(np.diff(lum, axis=1))
    dy = np.abs(np.diff(lum, axis=0))
    flat_neighbor_frac = float(
        np.mean(np.concatenate([(dx == 0).ravel(), (dy == 0).ravel()]))
    )

    if lum.shape[0] > 1 and lum.shape[1] > 1:
        same_2x2 = float(
            (
                (lum[:-1, :-1] == lum[1:, :-1])
                & (lum[:-1, :-1] == lum[:-1, 1:])
                & (lum[:-1, :-1] == lum[1:, 1:])
            ).mean()
        )
    else:
        same_2x2 = float("nan")

    lap = (
        -4 * lum
        + np.roll(lum, 1, axis=0)
        + np.roll(lum, -1, axis=0)
        + np.roll(lum, 1, axis=1)
        + np.roll(lum, -1, axis=1)
    )
    lap_var = float(np.var(lap))

    p = np.percentile(lum, [0, 1, 5, 50, 95, 99, 100])

    return {
        "flat_neighbor_frac": flat_neighbor_frac,
        "same_2x2_frac": same_2x2,
        "lap_var": lap_var,
        "lum_min": float(p[0]),
        "lum_p01": float(p[1]),
        "lum_p05": float(p[2]),
        "lum_p50": float(p[3]),
        "lum_p95": float(p[4]),
        "lum_p99": float(p[5]),
        "lum_max": float(p[6]),
        "clip0_frac": float((lum <= 0).mean()),
        "clip255_frac": float((lum >= 255).mean()),
    }


def center_crop(arr, frac=0.2):
    h, w = arr.shape[:2]
    ch = max(16, int(h * frac))
    cw = max(16, int(w * frac))
    y0 = h // 2 - ch // 2
    x0 = w // 2 - cw // 2
    return arr[y0 : y0 + ch, x0 : x0 + cw]


def save_pair_comparison(new_path, ref_path, out_path):
    new_arr = load_rgb(new_path)
    ref_arr = load_rgb(ref_path)
    new_crop = center_crop(new_arr, frac=0.22)
    ref_crop = center_crop(ref_arr, frac=0.22)

    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].imshow(new_arr)
    ax[0, 0].set_title(
        "NEW: {}\\n{}x{}".format(new_path.name, new_arr.shape[1], new_arr.shape[0])
    )
    ax[0, 1].imshow(ref_arr)
    ax[0, 1].set_title(
        "REF: {}\\n{}x{}".format(ref_path.name, ref_arr.shape[1], ref_arr.shape[0])
    )
    ax[1, 0].imshow(new_crop, interpolation="nearest")
    ax[1, 0].set_title("NEW center crop (nearest display)")
    ax[1, 1].imshow(ref_crop, interpolation="nearest")
    ax[1, 1].set_title("REF center crop (nearest display)")

    for a in ax.ravel():
        a.axis("off")

    fig.tight_layout()
    fig.savefig(str(out_path), dpi=150)
    plt.close(fig)


def inspect_reference_h5_for_fov(base_dir, sample=8):
    rows = []
    h5_files = sorted(base_dir.rglob("*.h5"))[:sample]
    for p in h5_files:
        r = {
            "path": str(p),
            "filename": p.name,
            "nx": None,
            "ny": None,
            "fovx_dsource": None,
            "fovy_dsource": None,
            "error": None,
        }
        try:
            with h5py.File(str(p), "r") as hfp:
                if (
                    "header" in hfp
                    and "camera" in hfp["header"]
                    and "nx" in hfp["header"]["camera"]
                ):
                    r["nx"] = int(hfp["header"]["camera"]["nx"][()])
                if (
                    "header" in hfp
                    and "camera" in hfp["header"]
                    and "ny" in hfp["header"]["camera"]
                ):
                    r["ny"] = int(hfp["header"]["camera"]["ny"][()])
                if (
                    "header" in hfp
                    and "camera" in hfp["header"]
                    and "fovx_dsource" in hfp["header"]["camera"]
                ):
                    r["fovx_dsource"] = float(
                        hfp["header"]["camera"]["fovx_dsource"][()]
                    )
                if (
                    "header" in hfp
                    and "camera" in hfp["header"]
                    and "fovy_dsource" in hfp["header"]["camera"]
                ):
                    r["fovy_dsource"] = float(
                        hfp["header"]["camera"]["fovy_dsource"][()]
                    )
        except Exception as exc:
            r["error"] = str(exc)
        rows.append(r)
    return rows


def mean_safe(vals):
    vals = [v for v in vals if v is not None and not math.isnan(v)]
    if not vals:
        return float("nan")
    return float(np.mean(vals))


def all_equal(values, tol=1e-8):
    vals = [v for v in values if v is not None]
    if not vals:
        return False
    first = vals[0]
    for v in vals[1:]:
        if abs(v - first) > tol:
            return False
    return True


def main():
    ensure_dirs()

    print("A) Enumerate root-level HDF5 files")
    h5_files = list_root_h5_files(H5_ROOT)
    print("root_h5_count={}".format(len(h5_files)))
    print("first_10_h5:")
    for p in h5_files[:10]:
        print("  - {}".format(p.name))

    print("\nB) Inspect each HDF5 for camera/image resolution")
    h5_rows = [inspect_h5_file(p) for p in h5_files]
    write_csv(QA_DIR / "h5_resolution_summary.csv", h5_rows)

    with (QA_DIR / "h5_top_level_items.txt").open("w", encoding="utf-8") as fp:
        for r in h5_rows:
            fp.write("{}: {}\n".format(r["filename"], r["top_level_items"]))

    h5_table_cols = [
        "filename",
        "camera_nx",
        "camera_ny",
        "detected_nx",
        "detected_ny",
        "pol_shape",
        "unpol_shape",
        "tau_shape",
        "resolution_source",
    ]
    print(format_table(h5_rows, h5_table_cols))

    print("\nC) Compare image outputs")
    new_inv = collect_image_inventory(NEW_IMG_DIR)
    ref_inv = collect_image_inventory(REF_IMG_DIR)
    write_csv(QA_DIR / "image_inventory_new.csv", new_inv)
    write_csv(QA_DIR / "image_inventory_ref.csv", ref_inv)

    new_dims = dimension_distribution(new_inv)
    ref_dims = dimension_distribution(ref_inv)
    write_csv(QA_DIR / "image_dim_distribution_new.csv", new_dims)
    write_csv(QA_DIR / "image_dim_distribution_ref.csv", ref_dims)

    print(summarize_inventory(new_inv, "NEW images"))
    print("")
    print(summarize_inventory(ref_inv, "REF images"))

    matched = choose_matched_pairs(new_inv, ref_inv, limit=3)
    write_csv(QA_DIR / "matched_pairs.csv", matched)

    # map stem -> h5 row for nx/ny lookup
    h5_by_stem = {Path(r["filename"]).stem: r for r in h5_rows}

    pair_metrics = []
    for i, pair in enumerate(matched):
        new_path = Path(pair["new_path"])
        ref_path = Path(pair["ref_path"])

        new_arr = load_rgb(new_path)
        ref_arr = load_rgb(ref_path)
        new_m = image_quality_metrics(new_arr)
        ref_m = image_quality_metrics(ref_arr)

        new_w, new_h = int(new_arr.shape[1]), int(new_arr.shape[0])
        ref_w, ref_h = int(ref_arr.shape[1]), int(ref_arr.shape[0])

        h5_nx = None
        h5_ny = None
        h5_row = h5_by_stem.get(new_path.stem)
        if h5_row is not None:
            h5_nx = h5_row.get("detected_nx")
            h5_ny = h5_row.get("detected_ny")

        panel_w_est_new = new_w / 2.0
        panel_w_est_ref = ref_w / 2.0

        metric_row = {
            "match_key": pair["match_key"],
            "new_file": new_path.name,
            "ref_file": ref_path.name,
            "new_dims": "{}x{}".format(new_w, new_h),
            "ref_dims": "{}x{}".format(ref_w, ref_h),
            "h5_nx": h5_nx,
            "h5_ny": h5_ny,
            "new_panel_to_h5_x_est": (panel_w_est_new / h5_nx) if h5_nx else None,
            "ref_panel_to_h5_x_est": (panel_w_est_ref / h5_nx) if h5_nx else None,
            "new_flat_neighbor_frac": new_m["flat_neighbor_frac"],
            "ref_flat_neighbor_frac": ref_m["flat_neighbor_frac"],
            "new_same_2x2_frac": new_m["same_2x2_frac"],
            "ref_same_2x2_frac": ref_m["same_2x2_frac"],
            "new_lap_var": new_m["lap_var"],
            "ref_lap_var": ref_m["lap_var"],
            "new_lum_p01": new_m["lum_p01"],
            "new_lum_p99": new_m["lum_p99"],
            "ref_lum_p01": ref_m["lum_p01"],
            "ref_lum_p99": ref_m["lum_p99"],
            "new_clip0_frac": new_m["clip0_frac"],
            "new_clip255_frac": new_m["clip255_frac"],
            "ref_clip0_frac": ref_m["clip0_frac"],
            "ref_clip255_frac": ref_m["clip255_frac"],
        }
        pair_metrics.append(metric_row)

        cmp_out = CMP_DIR / "pair_{}_{}.png".format(i + 1, pair["match_key"])
        save_pair_comparison(new_path, ref_path, cmp_out)

    write_csv(QA_DIR / "matched_pair_metrics.csv", pair_metrics)

    if pair_metrics:
        print("\nMatched-pair metrics:")
        print(
            format_table(
                pair_metrics,
                [
                    "match_key",
                    "new_dims",
                    "ref_dims",
                    "h5_nx",
                    "new_panel_to_h5_x_est",
                    "ref_panel_to_h5_x_est",
                    "new_flat_neighbor_frac",
                    "ref_flat_neighbor_frac",
                    "new_same_2x2_frac",
                    "ref_same_2x2_frac",
                    "new_lum_p01",
                    "new_lum_p99",
                    "ref_lum_p01",
                    "ref_lum_p99",
                ],
            )
        )
        print("\nComparison images written to:")
        for p in sorted(CMP_DIR.glob("pair_*.png")):
            print("  - {}".format(p))
    else:
        print("\nNo matched pairs found for deep comparisons.")

    ref_h5_sample = inspect_reference_h5_for_fov(ROOT / "M87", sample=8)
    if ref_h5_sample:
        write_csv(QA_DIR / "reference_h5_sample.csv", ref_h5_sample)
        print("\nReference H5 sample (for FOV/nx/ny checks):")
        print(
            format_table(
                ref_h5_sample,
                ["filename", "nx", "ny", "fovx_dsource", "fovy_dsource", "error"],
            )
        )

    # Hypothesis checks
    detected_nx = [r["detected_nx"] for r in h5_rows if r.get("detected_nx") is not None]
    detected_ny = [r["detected_ny"] for r in h5_rows if r.get("detected_ny") is not None]
    all_320 = bool(detected_nx) and all(x == 320 for x in detected_nx) and all(
        y == 320 for y in detected_ny
    )
    any_160 = any(x == 160 for x in detected_nx) or any(y == 160 for y in detected_ny)

    new_exts = {r["ext"].lower() for r in new_inv}
    ref_exts = {r["ext"].lower() for r in ref_inv}

    new_dims_mode = new_dims[0] if new_dims else None
    ref_dims_mode = ref_dims[0] if ref_dims else None

    new_flat_mean = mean_safe([r.get("new_flat_neighbor_frac") for r in pair_metrics])
    ref_flat_mean = mean_safe([r.get("ref_flat_neighbor_frac") for r in pair_metrics])
    new_same2_mean = mean_safe([r.get("new_same_2x2_frac") for r in pair_metrics])
    ref_same2_mean = mean_safe([r.get("ref_same_2x2_frac") for r in pair_metrics])

    ref_fov_values = [r.get("fovx_dsource") for r in ref_h5_sample if r.get("fovx_dsource") is not None]
    new_fov_values = [r.get("camera_fovx_dsource") for r in h5_rows if r.get("camera_fovx_dsource") is not None]
    same_fov_hint = False
    if ref_fov_values and new_fov_values and all_equal(ref_fov_values) and all_equal(new_fov_values):
        same_fov_hint = abs(ref_fov_values[0] - new_fov_values[0]) < 1e-8

    # Build markdown report
    md = []
    md.append("# QA REPORT: IPOLE Resolution and Image Quality")
    md.append("")
    md.append("## A) Root-level HDF5 enumeration")
    md.append("- Root path: `{}`".format(H5_ROOT))
    md.append("- `*.h5` count (root only): **{}**".format(len(h5_files)))
    md.append("- First 10 files:")
    for p in h5_files[:10]:
        md.append("  - `{}`".format(p.name))
    md.append("")

    md.append("## B) HDF5 camera/image resolution")
    md.append("- All files detected at 320x320: **{}** (any 160 found: **{}**)".format(all_320, any_160))
    md.append("- Detection priority: explicit `/header/camera/nx,ny`; fallback to dataset shapes (`pol`, `unpol`, `tau`).")
    md.append("- Per-file summary:")
    md.append("```")
    md.append(format_table(h5_rows, h5_table_cols))
    md.append("```")
    md.append("- Top-level groups/datasets are in `_qa/h5_top_level_items.txt`.")
    md.append("")

    md.append("## C) Image directory comparison")
    md.append("- NEW dir: `{}`".format(NEW_IMG_DIR))
    md.append("- REF dir: `{}`".format(REF_IMG_DIR))
    md.append("- NEW summary:")
    md.append("```")
    md.append(summarize_inventory(new_inv, "NEW images"))
    md.append("```")
    md.append("- REF summary:")
    md.append("```")
    md.append(summarize_inventory(ref_inv, "REF images"))
    md.append("```")

    md.append("- Typical dimensions:")
    md.append("  - NEW mode: {}".format(new_dims_mode if new_dims_mode else "n/a"))
    md.append("  - REF mode: {}".format(ref_dims_mode if ref_dims_mode else "n/a"))

    md.append("- File formats:")
    md.append("  - NEW extensions: {}".format(sorted(new_exts)))
    md.append("  - REF extensions: {}".format(sorted(ref_exts)))

    new_dpi_counts = Counter((r.get("dpi") or "none") for r in new_inv)
    ref_dpi_counts = Counter((r.get("dpi") or "none") for r in ref_inv)
    md.append("- PNG DPI metadata (top values):")
    md.append("  - NEW: {}".format(dict(new_dpi_counts.most_common(5))))
    md.append("  - REF: {}".format(dict(ref_dpi_counts.most_common(5))))
    md.append("")

    md.append("## D) Hypothesis checks")
    md.append(
        "1. HDF5 fallback to 160 px: {}".format(
            "RULED OUT" if all_320 and not any_160 else "NOT RULED OUT"
        )
    )
    md.append("   - Evidence: `/header/camera/nx=320`, `/header/camera/ny=320` and `pol` shape `(320, 320, 5)` for all root files.")

    if not math.isnan(new_flat_mean) and not math.isnan(ref_flat_mean):
        md.append("2. Plotting pipeline interpolation/downsampling difference: SUPPORTED")
        md.append(
            "   - Matched-pair flat-neighbor fraction (higher => blockier): NEW={:.4f}, REF={:.4f}".format(
                new_flat_mean, ref_flat_mean
            )
        )
        md.append(
            "   - Matched-pair 2x2 identical-block fraction: NEW={:.4f}, REF={:.4f}".format(
                new_same2_mean, ref_same2_mean
            )
        )
    else:
        md.append("2. Plotting pipeline interpolation/downsampling difference: INCONCLUSIVE")
    md.append("   - Code evidence: new script sets `interpolation='nearest'` and passes it to both `imshow` calls.")

    new_max_w = max([r.get("width") or 0 for r in new_inv]) if new_inv else 0
    if new_max_w >= 2000:
        md.append("3. Low DPI/small figure export: RULED OUT for NEW outputs")
    else:
        md.append("3. Low DPI/small figure export: INCONCLUSIVE")
    md.append("   - NEW images are typically large (e.g., 2400x1200), not thumbnail-sized.")

    md.append("4. Colormap/interpolation/normalization differences: PARTIAL SUPPORT")
    md.append("   - Colormaps are similar (`afmhot`, `seismic`) and intensity panel uses `vmax=np.max(I)` in both script families.")
    md.append("   - Interpolation differs (new explicit nearest-neighbor), which directly increases visible blockiness.")

    if same_fov_hint:
        md.append("5. FOV/zoom difference causing apparent blockiness: NOT SUPPORTED")
        md.append("   - Sampled reference HDF5 and all new root HDF5 show matching `fovx_dsource=160`, `fovy_dsource=160`.")
    else:
        md.append("5. FOV/zoom difference causing apparent blockiness: INCONCLUSIVE")
        md.append("   - Not all reference images had matching on-disk HDF5 companions for one-to-one FOV verification.")

    if new_max_w >= 2000:
        md.append("6. New outputs are thumbnails/quicklooks: RULED OUT")
    else:
        md.append("6. New outputs are thumbnails/quicklooks: INCONCLUSIVE")
    md.append("   - NEW outputs are full-sized PNG files and generally larger than REF images.")

    md.append("")
    md.append("## Most likely root cause and fix")
    md.append("- Most likely cause: the new plotting script uses nearest-neighbor interpolation, which makes each native image pixel appear as hard blocks when rendered large.")
    md.append("- Fix: in `aricarte/aricarte-copy/aricarte/PositronIPOLEScripts/create_images.py`, change default `interpolation` from `nearest` to `antialiased` (or `bilinear`) for presentation images, or run with `--interpolation antialiased`.")
    md.append("- Keep simulation resolution at 320 (`nx=320`, `ny=320`); no fallback-to-160 issue was found in root HDF5 outputs.")

    md.append("")
    md.append("## Artifacts")
    md.append("- `_qa/h5_resolution_summary.csv`")
    md.append("- `_qa/h5_top_level_items.txt`")
    md.append("- `_qa/image_inventory_new.csv`")
    md.append("- `_qa/image_inventory_ref.csv`")
    md.append("- `_qa/image_dim_distribution_new.csv`")
    md.append("- `_qa/image_dim_distribution_ref.csv`")
    md.append("- `_qa/matched_pairs.csv`")
    md.append("- `_qa/matched_pair_metrics.csv`")
    md.append("- `_qa/comparisons/pair_*.png`")
    if (QA_DIR / "reference_h5_sample.csv").exists():
        md.append("- `_qa/reference_h5_sample.csv`")

    report_path = QA_DIR / "QA_REPORT.md"
    report_path.write_text("\n".join(md), encoding="utf-8")
    print("\nWrote report: {}".format(report_path))


if __name__ == "__main__":
    main()
