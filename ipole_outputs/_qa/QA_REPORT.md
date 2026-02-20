# QA REPORT: IPOLE Resolution and Image Quality

## A) Root-level HDF5 enumeration
- Root path: `/work/vmo703/ipole_outputs`
- `*.h5` count (root only): **20**
- First 10 files:
  - `Ma+0.94_4000_CRITBETAWJET_f0.5_bc0.01_0.00.h5`
  - `Ma+0.94_4000_CRITBETAWJET_f0.5_bc0.1_0.00.h5`
  - `Ma+0.94_4000_CRITBETA_f0.5_bc0.01_0.00.h5`
  - `Ma+0.94_4000_CRITBETA_f0.5_bc0.1_0.00.h5`
  - `Ma-0.5_4000_CRITBETAWJET_f0.5_bc0.01_0.00.h5`
  - `Ma-0.5_4000_CRITBETAWJET_f0.5_bc0.1_0.00.h5`
  - `Ma-0.5_4000_CRITBETA_f0.5_bc0.01_0.00.h5`
  - `Ma-0.5_4000_CRITBETA_f0.5_bc0.1_0.00.h5`
  - `Sa+0.94_4000_CRITBETAWJET_f0.5_bc0.01_0.00.h5`
  - `Sa+0.94_4000_CRITBETAWJET_f0.5_bc0.1_0.00.h5`

## B) HDF5 camera/image resolution
- All files detected at 320x320: **True** (any 160 found: **False**)
- Detection priority: explicit `/header/camera/nx,ny`; fallback to dataset shapes (`pol`, `unpol`, `tau`).
- Per-file summary:
```
filename                                      | camera_nx | camera_ny | detected_nx | detected_ny | pol_shape     | unpol_shape | tau_shape  | resolution_source                   
----------------------------------------------+-----------+-----------+-------------+-------------+---------------+-------------+------------+-------------------------------------
Ma+0.94_4000_CRITBETAWJET_f0.5_bc0.01_0.00.h5 | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Ma+0.94_4000_CRITBETAWJET_f0.5_bc0.1_0.00.h5  | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Ma+0.94_4000_CRITBETA_f0.5_bc0.01_0.00.h5     | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Ma+0.94_4000_CRITBETA_f0.5_bc0.1_0.00.h5      | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Ma-0.5_4000_CRITBETAWJET_f0.5_bc0.01_0.00.h5  | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Ma-0.5_4000_CRITBETAWJET_f0.5_bc0.1_0.00.h5   | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Ma-0.5_4000_CRITBETA_f0.5_bc0.01_0.00.h5      | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Ma-0.5_4000_CRITBETA_f0.5_bc0.1_0.00.h5       | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa+0.94_4000_CRITBETAWJET_f0.5_bc0.01_0.00.h5 | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa+0.94_4000_CRITBETAWJET_f0.5_bc0.1_0.00.h5  | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa+0.94_4000_CRITBETA_f0.5_bc0.01_0.00.h5     | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa+0.94_4000_CRITBETA_f0.5_bc0.1_0.00.h5      | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa+0.94_5000_CRITBETA_f0.5_bc0.01_0.00.h5     | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa+0.94_5000_CRITBETA_f0.5_bc0.1_0.00.h5      | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa-0.5_4000_CRITBETAWJET_f0.5_bc0.01_0.00.h5  | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa-0.5_4000_CRITBETAWJET_f0.5_bc0.1_0.00.h5   | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa-0.5_4000_CRITBETA_f0.5_bc0.01_0.00.h5      | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa-0.5_4000_CRITBETA_f0.5_bc0.1_0.00.h5       | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa-0.5_5000_CRITBETA_f0.5_bc0.01_0.00.h5      | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
Sa-0.5_5000_CRITBETA_f0.5_bc0.1_0.00.h5       | 320       | 320       | 320         | 320         | (320, 320, 5) | (320, 320)  | (320, 320) | /header/camera/nx, /header/camera/ny
```
- Top-level groups/datasets are in `_qa/h5_top_level_items.txt`.

## C) Image directory comparison
- NEW dir: `/work/vmo703/ipole_outputs/images`
- REF dir: `/work/vmo703/ipole_outputs/M87/images`
- NEW summary:
```
NEW images: count=20
NEW images: mean_size=219.8KB, median_size=229.0KB
NEW images: file types
ext | count
----+------
png | 20   
NEW images: dimensions (top 10)
width | height | count
------+--------+------
2400  | 1200   | 20   
NEW images: size examples
filename                                      | dims      | size    | dpi          
----------------------------------------------+-----------+---------+--------------
Sa-0.5_4000_CRITBETAWJET_f0.5_bc0.1_0.00.png  | 2400x1200 | 168.6KB | 300.00,300.00
Sa-0.5_4000_CRITBETAWJET_f0.5_bc0.01_0.00.png | 2400x1200 | 178.7KB | 300.00,300.00
Sa+0.94_4000_CRITBETA_f0.5_bc0.1_0.00.png     | 2400x1200 | 263.5KB | 300.00,300.00
Sa+0.94_5000_CRITBETA_f0.5_bc0.1_0.00.png     | 2400x1200 | 264.2KB | 300.00,300.00
```
- REF summary:
```
REF images: count=100
REF images: mean_size=125.2KB, median_size=114.2KB
REF images: file types
ext | count
----+------
png | 100  
REF images: dimensions (top 10)
width | height | count
------+--------+------
800   | 400    | 100  
REF images: size examples
filename                                          | dims    | size    | dpi          
--------------------------------------------------+---------+---------+--------------
output_Sa-0.5_5000_model_RBETA_Rhigh_40_0.000.png | 800x400 | 76.2KB  | 100.00,100.00
output_Ma-0.5_5000_model_RBETA_Rhigh_1_1.000.png  | 800x400 | 82.0KB  | 100.00,100.00
output_Sa-0.5_5000_model_RBETA_Rhigh_20_0.000.png | 800x400 | 186.7KB | 100.00,100.00
output_Sa-0.5_5000_model_RBETA_Rhigh_10_0.000.png | 800x400 | 203.8KB | 100.00,100.00
```
- Typical dimensions:
  - NEW mode: {'width': 2400, 'height': 1200, 'count': 20}
  - REF mode: {'width': 800, 'height': 400, 'count': 100}
- File formats:
  - NEW extensions: ['png']
  - REF extensions: ['png']
- PNG DPI metadata (top values):
  - NEW: {'300.00,300.00': 20}
  - REF: {'100.00,100.00': 100}

## D) Hypothesis checks
1. HDF5 fallback to 160 px: RULED OUT
   - Evidence: `/header/camera/nx=320`, `/header/camera/ny=320` and `pol` shape `(320, 320, 5)` for all root files.
2. Plotting pipeline interpolation/downsampling difference: SUPPORTED
   - Matched-pair flat-neighbor fraction (higher => blockier): NEW=0.9263, REF=0.7849
   - Matched-pair 2x2 identical-block fraction: NEW=0.8702, REF=0.6962
   - Code evidence: new script sets `interpolation='nearest'` and passes it to both `imshow` calls.
3. Low DPI/small figure export: RULED OUT for NEW outputs
   - NEW images are typically large (e.g., 2400x1200), not thumbnail-sized.
4. Colormap/interpolation/normalization differences: PARTIAL SUPPORT
   - Colormaps are similar (`afmhot`, `seismic`) and intensity panel uses `vmax=np.max(I)` in both script families.
   - Interpolation differs (new explicit nearest-neighbor), which directly increases visible blockiness.
5. FOV/zoom difference causing apparent blockiness: NOT SUPPORTED
   - Sampled reference HDF5 and all new root HDF5 show matching `fovx_dsource=160`, `fovy_dsource=160`.
6. New outputs are thumbnails/quicklooks: RULED OUT
   - NEW outputs are full-sized PNG files and generally larger than REF images.

## Most likely root cause and fix
- Most likely cause: the new plotting script uses nearest-neighbor interpolation, which makes each native image pixel appear as hard blocks when rendered large.
- Fix: in `aricarte/aricarte-copy/aricarte/PositronIPOLEScripts/create_images.py`, change default `interpolation` from `nearest` to `antialiased` (or `bilinear`) for presentation images, or run with `--interpolation antialiased`.
- Keep simulation resolution at 320 (`nx=320`, `ny=320`); no fallback-to-160 issue was found in root HDF5 outputs.

## Artifacts
- `_qa/h5_resolution_summary.csv`
- `_qa/h5_top_level_items.txt`
- `_qa/image_inventory_new.csv`
- `_qa/image_inventory_ref.csv`
- `_qa/image_dim_distribution_new.csv`
- `_qa/image_dim_distribution_ref.csv`
- `_qa/matched_pairs.csv`
- `_qa/matched_pair_metrics.csv`
- `_qa/comparisons/pair_*.png`
- `_qa/reference_h5_sample.csv`