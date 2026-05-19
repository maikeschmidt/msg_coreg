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
├── coreg_path.m
├── cr_add_functions.m
├── cr_check_registration.m
├── cr_generate_sensor_array_v4.m
├── cr_generate_spine_center.m
├── cr_get_fids.m
├── cr_load_meshes.m
├── cr_register_brain.m
├── cr_register_torso.m
├── example/
│   ├── example_script_1.m
│   └── example_script_2.m
├── meshes/
│   ├── back_muscle_temp.stl
│   ├── canonical_cervical_cont.stl
│   ├── canonical_cervical_homo.stl
│   ├── canonical_cervical_inhomo.stl
│   ├── canonical_full_cont.stl
│   ├── canonical_full_homo.stl
│   ├── canonical_full_inhomo.stl
│   ├── canonical_heart.stl
│   ├── canonical_lungs.stl
│   ├── canonical_torso.stl
│   ├── cervical_spine.stl
│   ├── heart.stl
│   ├── mri_cervical_cont.stl
│   ├── mri_cervical_homo.stl
│   ├── mri_cervical_inhomo.stl
│   ├── mri_cervical_spine.stl
│   ├── mri_full_cont.stl
│   ├── mri_full_homo.stl
│   ├── mri_full_inhomo.stl
│   ├── mri_full_spine.stl
│   ├── mri_lungs.stl
│   ├── mri_torso.stl
│   ├── realistic_cervical_bone.stl
│   ├── realistic_full_bone.stl
│   ├── spine.stl
│   └── vagus_nerve_temp.stl
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

An example optical scan is provided at `meshes/surface.stl`.

> For accurate spinal cord positioning, use the anatomical meshes together 
> with the provided `surface.stl`. If you have your own sensor array, provide 
> your own optical scan and use the canonical meshes instead.

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
> uncomment **line 345** in `cr_check_registration.m`.

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

> *[Insert forward modelling paper citation here]*

---

## Optional: Forward Model Sensitivity Analysis

Both example scripts include optional sections for generating shifted geometry 
files that can be used to assess how sensitive forward solutions are to 
registration uncertainty. These sections are **self-contained and clearly 
labelled** within each script — they can be run or skipped independently of 
the main coregistration workflow.

Two types of sensitivity analysis are supported, corresponding to two 
different sources of registration error:

### Source position sensitivity

Evaluates uncertainty in spinal cord localisation by shifting the source 
model by small fixed amounts independently along each anatomical axis 
(±2, ±4, ±6 mm in X, Y, and Z). This produces 18 shifted geometry files 
plus the original (19 total), each with an identical mesh and sensor array 
but a translated source model.

This is useful for quantifying how much the predicted sensor pattern 
changes if the spinal cord centre line is misregistered by a few millimetres.

**When to use:** when you want to assess the impact of anatomical 
uncertainty on forward model accuracy, for example when using canonical 
meshes where spinal cord positioning is approximate.

### Sensor array sensitivity

Evaluates uncertainty in sensor array registration by shifting the entire 
sensor array by random 3D displacements [dx, dy, dz]. Shifts are grouped 
into three bundles representing different registration error scales 
(~2 mm, ~5 mm, ~10 mm), with 8 random realisations per bundle. This 
produces 24 shifted geometry files plus the original (25 total), each with 
an identical mesh and source model but a translated sensor array.

Sensor orientations (`coilori`, `chanori`) and the transfer matrix (`tra`) 
are **not modified** — only `coilpos` and `chanpos` are shifted, so the 
triaxial orthogonal structure of the sensor array is fully preserved.

Shifts are generated with `rng(42)` for reproducibility. The exact 
[dx, dy, dz] vectors are printed at runtime and can be hardcoded in the 
script for exact reproduction across machines.

**When to use:** when you want to assess how sensitive forward solutions 
are to errors in sensor array placement or body scan registration, for 
example when the sensor-to-body transform has limited accuracy.

### Running the sensitivity analyses

Both sensitivity sections are at the end of each example script and can 
be run after the main coregistration workflow completes. The geometry 
`.mat` files they produce are passed directly to `run_bem_leadfields.m` 
in `msg_fwd` for leadfield computation, and then analysed using the 
sensitivity pipeline in `msg_fwd`. No additional configuration of 
`msg_coreg` is required.

Full workflow:

```matlab
% 1. Run the main coregistration workflow (example_script_1 or _2)
% 2. Run the sensitivity section(s) at the end of the same script
%    — these save geometry .mat files to the same output folder
% 3. In msg_fwd: run BEM leadfields for the shifted geometry files
% 4. In msg_fwd: run compute_sensitivity_rsq, then plot/table scripts
```

See the `msg_fwd` README for the full sensitivity analysis pipeline:  
https://github.com/maikeschmidt/msg_fwd

---

## Example Scripts

### example_script_1.m — Register meshes with an existing sensor array

Demonstrates how to register canonical or anatomical simulation meshes into 
experimental sensor space and import an existing experimental OPM sensor layout. 
Recommended when you already have an experimentally defined sensor layout and 
want to run simulations in the same coordinate system as recorded data.

**Optional sensitivity analysis sections** (at the end of the script):

**Source position sensitivity** — generates 19 geometry files (1 original 
+ 18 shifted) with source positions translated by ±2, ±4, and ±6 mm 
independently along X, Y, and Z. The meshes and sensor array are identical 
across all configurations.

The 19 geometry files produced are:
geometries_original.mat
geometries_shift_x_pos2mm.mat
geometries_shift_x_pos4mm.mat
geometries_shift_x_pos6mm.mat
geometries_shift_x_neg2mm.mat
geometries_shift_x_neg4mm.mat
geometries_shift_x_neg6mm.mat
geometries_shift_y_pos2mm.mat
geometries_shift_y_pos4mm.mat
geometries_shift_y_pos6mm.mat
geometries_shift_y_neg2mm.mat
geometries_shift_y_neg4mm.mat
geometries_shift_y_neg6mm.mat
geometries_shift_z_pos2mm.mat
geometries_shift_z_pos4mm.mat
geometries_shift_z_pos6mm.mat
geometries_shift_z_neg2mm.mat
geometries_shift_z_neg4mm.mat
geometries_shift_z_neg6mm.mat

**Sensor array sensitivity** — generates 25 geometry files (1 original + 
24 shifted) with the entire sensor array translated by random [dx, dy, dz] 
displacements in three bundles by error scale. The meshes and source model 
are identical across all configurations.

> **Note:** Both sensitivity sections require an experimental sensor array 
> saved as `experimental_sensors` in the geometry struct, which is set up 
> earlier in this script.

### example_script_2.m — Build anatomical meshes and generate a sensor array

Demonstrates the full anatomical modelling pipeline using subject-specific 
geometry, realistic MRI-segmented bone, and scanner-cast optical surface 
(`surface.stl`). Reproduces the simulation setup used in the publication. 
Recommended when accurate spinal cord positioning or realistic bone geometry 
is required.

Also includes the same optional sensitivity analysis sections as 
`example_script_1.m`, allowing sensitivity analyses to be run from either 
the canonical or anatomical modelling workflow.

---

---

## Citation

If you use this toolbox in your work, please cite:

> *[Add paper citation here]*

---


## Copyright

Copyright (c) 2026 University College London  
Department of Imaging Neuroscience  
Author: Maike Schmidt — maike.schmidt.23@ucl.ac.uk  
Repository: https://github.com/maikeschmidt/msg_coreg
