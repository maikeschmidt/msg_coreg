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

## Example Scripts

### example_script_1.m — Register meshes with an existing sensor array

Demonstrates how to register canonical or anatomical simulation meshes into 
experimental sensor space and import an existing sensor layout (MEG/OPM or EEG). 
Recommended when you already have an experimentally defined sensor layout and 
want to run simulations in the same coordinate system as recorded data.

### example_script_2.m — Build anatomical meshes and generate a sensor array

Demonstrates the full anatomical modelling pipeline using subject-specific 
geometry, realistic MRI-segmented bone, and scanner-cast optical surface 
(`surface.stl`). Reproduces the simulation setup used in the publication. 
Recommended when accurate spinal cord positioning or realistic bone geometry 
is required.

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
