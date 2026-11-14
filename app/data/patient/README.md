# Patient geometry served by the app

`AHMU1218142/` holds the six CFD-ready lumens shown in the app's **Patient CTA** mode:
* `<variant>.stl`: binary STL, mm, faces grouped by patch (wall, inlet, outlet_1, outlet_2)
* `<variant>.json`: patch counts, the M1 centerline station of each wall face, and reduced-order profiles
  along M1
* `manifest.json`: labels, predictions and attribution

**Source and licence.** These surfaces are derived from subject AHMU1218142 of the CMHA dataset:
Song M, Wang S, Qian Q, Zhou Y, Luo Y, Gong X. *Intracranial aneurysm CTA images and 3D models dataset
with clinical morphological and hemodynamic data.* Sci Data 2024. doi:10.1038/s41597-024-04056-8.
Data: Figshare doi:10.6084/m9.figshare.26965450, licensed **CC BY 4.0**.

**Changes made here.** Flow extensions and caps were added. Lesions (70% plaque, 95% thrombus) and
post-treatment lumens (thrombectomy, angioplasty, stent) were imposed computationally. The
post-treatment lumens are prescribed outcomes, not simulated device mechanics. Regenerate with
`python3 scripts/cmha_variants.py AHMU1218142`, which needs the dataset downloaded to `data/cmha/`.

**What is not committed.** The CTA head scans (NIfTI, about 100 MB each, 4.5 GB in total), the
authors' Fluent cases and the other 43 controls stay in the git-ignored `/data/` folder. They exceed
GitHub's file limits, and the app does not need them.
