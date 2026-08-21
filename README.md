# msg_coreg Toolbox

**MSG Coregistration Toolbox** — tools for generating anatomically informed mesh 
models for spinal cord simulations and concurrent cortico–spinal interaction studies.

Developed by **Maike Schmidt** at the **Department of Imaging Neuroscience, 
University College London**.

> For questions, issues, or contributions, please open an issue or pull request on GitHub.  
> Contact: maike.schmidt.23@ucl.ac.uk

---

## Directory Structure


```
msg_coreg/
├── coreg_path.m                    — path function (locates repository root)
├── cr_add_functions.m              — dependency setup (HBF, FieldTrip wrappers)
├── cr_check_registration.m         — main entry point: register, plot, return meshes
├── cr_load_meshes.m                — load the STL set and apply a transform
├── cr_register_torso.m             — fit the torso mesh to a subject surface
├── cr_register_brain.m             — fit the SPM brain/skull/scalp meshes
├── cr_get_fids.m                   — canonical / anatomical fiducial definitions
├── cr_generate_sensor_array_v4.m   — generate an OPM or electrode array
├── cr_generate_spine_center.m      — spinal cord centreline source model
├── determine_body_scan_units.mlx   — helper notebook: work out a scan's units
│
├── warping/                        — synthesise alternative body shapes
│   ├── cr_generate_warps.m         — reproducible set of affine warps
│   ├── cr_plot_warps.m             — INSPECT the warps before using them
│   ├── cr_build_warp_geometries.m  — apply each warp, write geometry files
│   ├── cr_find_intersections.m     — name the meshes that self-intersect
│   └── cr_scan_intersections.m     — run that check over a whole folder
│
├── example/
│   ├── example_script_1.m          — canonical workflow, existing sensor array
│   └── example_script_2            — anatomical workflow, generated array
│
├── meshes/                         — bundled STL meshes (see meshes/README.md)
└── README.md
```

---

## Overview

This toolbox supports both **canonical** and **anatomical** modelling approaches 
and is designed to integrate with MEG/OPM, EEG, and surface electrode simulations.

It allows you to:

- Generate torso, spinal cord, bone, and (optionally) brain meshes
- Register meshes into experimental sensor space
- Create or import sensor arrays (OPMs or surface electrodes)
- Export meshes and source models for forward modelling (BEM/FEM)

The core motivation is to **investigate how different bone geometries affect 
spinal cord forward modelling**, while enabling **simultaneous cortical and 
spinal simulations**.

---

## Requirements

1. **MATLAB** (R2020a or later recommended)

2. **SPM** — the developmental version is recommended  
   https://www.fil.ion.ucl.ac.uk/spm/

3. **FieldTrip** — required for sensor formatting and headshape reading  
   https://www.fieldtriptoolbox.org/

4. **Helsinki BEM Framework (HBF)** by Matti Stenroos  
   Add as a subfolder named `hbf_lc_p` inside this repository:  
   https://github.com/MattiStenroos/hbf_lc_p/tree/master/hbf_calc

5. **Optical / 3D surface scan of the participant**  
   Acquired in the experimental setup or scanner cast (depending on model choice)

---

## Getting Started

```matlab
% 1. Add the toolbox and all dependencies to your MATLAB path
cr_add_functions()

% 2. Set up your input struct and run the registration check
S.subject    = your_subject_mesh;   % struct with .vertices and .faces
S.torso_mode = 'canonical';         % or 'anatomical'
S.spine_mode = 'full';
S.bone_mode  = 'homo';

output_meshes = cr_check_registration(S);
```

See the `example/` folder for full worked workflows.

---

## Modelling Approaches

### 1. Canonical Model

Uses **canonical simulation meshes** with an optical/3D scan of the participant 
in the experimental setup. The user manually selects three fiducial points on 
the scan (left shoulder, right shoulder, chin) to transform the canonical meshes 
into experimental sensor space.

> **Note:** Canonical meshes are based on a seated subject, so spinal cord 
> localisation is approximate. This approach is suitable when subject-specific 
> MRI is unavailable.

### 2. Anatomical Model

Uses **subject-specific anatomical information** based on a custom-built MSG 
scanner cast designed from an anatomical MRI. The transform from MRI space to 
experimental sensor space is known.

The anatomical workflow expects a scanner-cast optical surface (`surface.stl`).
This scan is **not bundled** in `meshes/` — supply your own for your setup.

> For accurate spinal cord positioning, use the anatomical meshes together 
> with your scanner-cast `surface.stl`. If you have your own sensor array, 
> provide your own optical scan and use the canonical meshes instead.

---

## Bone Model Variants

A key feature of this toolbox is support for multiple bone geometries:

| Variant | Canonical | Anatomical |
|---|---|---|
| Continuous | ✓ | ✓ |
| Homogeneous toroidal | ✓ | ✓ |
| Inhomogeneous toroidal | ✓ | ✓ |
| Realistic MRI-segmented | ✗ | ✓ |

---

## Sensor Arrays

You can either **import** an existing experimental sensor array or **generate** 
one using the toolbox.

### Importing Experimental Sensor Arrays

Experimental sensor layouts can be imported directly. An example using SPM 
sensor structures is provided in `example/example_script_1.m`.

### Generating Sensor Arrays

Supported sensor types:

- **Magnetometers (OPMs)** — triaxial sensors aligned to the Cartesian 
  coordinate system (Z-axis labelled as radial due to mesh orientation)
- **Surface electrodes** — dual-axis electrodes with common-average reference

Supported array configurations:

- Front-only, back-only, or full 360° torso array
  (full torso uses surface normals as the radial direction)

Customisable parameters:

| Parameter | OPM default | Electrode default |
|---|---|---|
| Sensor spacing | 30 mm | 30 mm |
| Offset from body | 10 mm | 0 mm |
| Coverage (top/bottom/left/right) | 0.6 | 0.6 |

---

## Optional Brain Registration

To investigate **concurrent cortico–spinal interactions**, a brain model can 
be included using the SPM template brain. This requires selection of three 
fiducials: left preauricular, right preauricular, and nasion.

> To export the transformation matrix applied to the SPM brain template,
> uncomment the `output_meshes.brains_transform_mat` line near the end of
> `cr_check_registration.m`.

---

## Spinal Cord Source Model

`cr_generate_spine_center()` identifies the centreline of the spinal cord and 
places candidate source points along it. This step is optional and only required 
for simulating distributed spinal sources.

---

## Forward Modelling

For BEM forward modelling, export the following outputs to your pipeline:

- All registered meshes
- Spinal cord source locations
- The transformation matrix

Compatible forward modelling pipeline:
https://github.com/maikeschmidt/msg_fwd

---

## Optional: Forward Model Sensitivity Analysis

Both example scripts end with **optional, self-contained sections** that
generate shifted geometry files, so you can measure how sensitive a forward
solution is to registration uncertainty. Run or skip them independently of
the main coregistration workflow.

| Section | Models | Produces |
|---|---|---|
| **Source position sensitivity** | uncertainty in spinal cord localisation | 19 geometry files: the original plus fixed shifts of ±2, ±4, ±6 mm applied independently along X, Y and Z. Meshes and sensor array identical throughout. |
| **Sensor array sensitivity** | registration error in sensor placement | 25 geometry files: the original plus 24 random `[dx, dy, dz]` displacements of the whole array, in three bundles (~2 mm, ~5 mm, ~10 mm) with 8 realisations each. Meshes and source model identical throughout. |

For the sensor shifts, only `coilpos` and `chanpos` move; orientations
(`coilori`, `chanori`) and the transfer matrix (`tra`) are left alone, so the
triaxial orthogonal structure of the array is preserved. Shifts use `rng(42)`
and the exact vectors are printed at runtime, so a run can be reproduced
exactly on another machine.

Both sections require an experimental sensor array saved as
`experimental_sensors` in the geometry struct, which the example scripts set
up earlier.

Full workflow:

```matlab
% 1. Run the main coregistration workflow (example_script_1 or _2)
% 2. Run the sensitivity section(s) at the end of the same script
%    — these save geometry .mat files to the same output folder
% 3. In msg_fwd: run BEM leadfields for the shifted geometry files
% 4. In msg_pert: run pt_load_leadfields, pt_compute_rsq, then plot/table scripts
```

**For anything more than a quick check, use msg_pert instead.** It has more
systematic generators of its own (`pt_generate_source_shifts`,
`pt_generate_sensor_shifts` — 3 bundles x 8 random shifts each) plus
conductivity perturbation, and it is where the analysis lives:
https://github.com/maikeschmidt/msg_pert

---


## Example Scripts

### example_script_1.m — Register meshes with an existing sensor array

Demonstrates how to register canonical or anatomical simulation meshes into 
experimental sensor space and import an existing experimental OPM sensor layout. 
Recommended when you already have an experimentally defined sensor layout and 
want to run simulations in the same coordinate system as recorded data.

It also carries both optional sensitivity sections at the end — see
[Optional: Forward Model Sensitivity Analysis](#optional-forward-model-sensitivity-analysis)
above. Source-shift files are named `geometries_shift_<axis>_<sign><n>mm.mat`
alongside `geometries_original.mat`.

### example_script_2 — Build anatomical meshes and generate a sensor array

Demonstrates the full anatomical modelling pipeline using subject-specific 
geometry, realistic MRI-segmented bone, and a scanner-cast optical surface
(user-supplied `surface.stl`). Recommended when accurate spinal cord
positioning or realistic bone geometry is required.

Also includes the same optional sensitivity analysis sections as 
`example_script_1.m`, allowing sensitivity analyses to be run from either 
the canonical or anatomical modelling workflow.

---

## Synthesising Alternative Body Shapes

The `warping/` folder generates affine warps of the anatomical model —
taller/thinner, shorter/wider, and everything in between — and writes one
geometry file per warp. Use it to test whether a forward-modelling result
depends on one particular set of body proportions.

```matlab
W = cr_generate_warps(struct('n_warps', 30));

cr_plot_warps('geometries_original.mat');   % ALWAYS look at this first

S = struct('geom_file', 'geometries_anatom_full_realistic.mat', ...
           'warp_file', 'anatomical_warps.mat', ...
           'outdir',    'path/to/warp/geometries');
files = cr_build_warp_geometries(S);
```

Warps are affine, so vertebra count, spinal curvature and relative organ
placement are inherited unchanged from the source anatomy. They bound
robustness to body-shape variation; they are **not** a substitute for
scanning more participants, and should never be described as N subjects.

`cr_scan_intersections` checks a whole folder of geometries for the
self-intersections that make TetGen fail, and names the meshes responsible —
worth running before committing hours of FEM compute.

---

## Companion Repositories

| Repository | Role |
|---|---|
| [msg_fwd](https://github.com/maikeschmidt/msg_fwd) | Forward modelling — BEM and FEM lead fields from these geometries |
| [msg_pert](https://github.com/maikeschmidt/msg_pert) | Perturbation analysis of those forward models |
| [msg_analysis](https://github.com/maikeschmidt/msg_analysis) | Concurrent cortico-spinal OPM data analysis |
| [msg_error](https://github.com/maikeschmidt/msg_error) | MRI-derived head/neck geometries and error analysis |

---

## Copyright

Copyright (c) 2026 University College London  
Department of Imaging Neuroscience  
Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk  
Repository: https://github.com/maikeschmidt/msg_coreg
