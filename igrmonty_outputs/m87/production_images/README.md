# Production images (Phase-4 era, run_2026-09-15)

Home for production-quality IPOLE images of the Phase-4 branches.

- Which branches are production-grade: see `../_qa/p4_passfail.md` (live;
  auto-updated as spectra land). Use the tuned M_unit from each branch's
  spectrum h5 (`/params/M_unit`) or data/final_grmonty_paper.csv.
- IPOLE pars: /work/vmo703/ipole_pars/ ; workflow: notebooks/creating_images.ipynb
  (pattern also in sgrA/create_images.py).

** WARNING before generating anything labeled production: the ipole working
copy still needs BOTH one-line gates to match these spectra —
(1) constant_beta_paper_literal (P_B = B^2/8pi, model/iharm/model.c:547-class
fix) and (2) crit_floor = 3e-2. Until those land in ipole, images of wJET or
Crit-beta branches will NOT be consistent with the Phase-4 spectra. **
