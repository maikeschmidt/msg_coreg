# msg_coreg Toolbox

The **msg_coreg** toolbox provides tools to generate anatomically informed mesh models for **spinal cord simulations** and **concurrent cortico–spinal interaction studies**.  
It supports both **canonical** and **anatomical** modelling approaches and is designed to integrate with MEG/OPM, EEG, and surface electrode simulations.

---

## Overview

This toolbox allows you to:

- Generate torso, spinal cord, bone, and (optionally) brain meshes
- Register meshes into experimental sensor space
- Create or import sensor arrays (OPMs or surface electrodes)
- Export meshes and source models for forward modelling (BEM/FEM)

The core motivation behind this toolbox is to **investigate how different bone geometries affect spinal cord forward modelling**, while also enabling **simultaneous cortical and spinal simulations**.

---

## Requirements

To use this toolbox, you will need:

1. **SPM**  
   - The developmental version of SPM is recommended.

2. **Helsinki BEM Framework (HBF)** by Matti Stenroos  
   - Add the repository as a subfolder named `hbf_lc_p` inside this repository  
   - Source: https://github.com/MattiStenroos/hbf_lc_p/tree/master/hbf_calc

3. **Optical / 3D surface scan of the participant**  
   - Acquired in the experimental setup or scanner cast (depending on model choice)

---

## Modelling Approaches

Two modelling approaches are supported:

### 1. Canonical Model

- Uses **canonical simulation meshes**
- Requires an optical/3D scan of the participant in the experimental setup (MSG or ESG)
- The user manually selects **three fiducial points** on the scan:
  - Left shoulder
  - Right shoulder
  - Chin
- These fiducials are used to transform the canonical meshes into experimental sensor space

Note:  
Canonical meshes are based on a **seated subject**, so spinal cord localisation is approximate.  
However, this approach provides a good estimate and is suitable when subject-specific MRI is unavailable.

---

### 2. Anatomical Model

- Uses **subject-specific anatomical information**
- Based on a **custom-built MSG scanner cast**, designed from an anatomical MRI
- The transform from MRI space to experimental sensor space is known
- An example optical scan is provided:
  - `meshes/surface.stl`

If you want accurate spinal cord positioning within the torso, **use the anatomical meshes together with the provided `surface.stl`**.

If you have your own sensor array or setup, provide your own optical scan in sensor space and use the **canonical meshes** instead.

---

## Sensor Arrays

Once you have chosen your modelling approach, you must decide whether to:

1. **Generate a sensor array**, or  
2. **Import an experimental sensor array**

### Importing Experimental Sensor Arrays

- Experimental sensor layouts can be imported directly
- An example using SPM sensor structures is provided

### Generating Sensor Arrays

Supported sensor types:

- **Magnetometers (OPMs)**  
  - Tri-axial sensors aligned to the Cartesian coordinate system
  - The **Z-axis is labelled as radial** due to mesh orientation (can be modified)

- **Surface electrodes**  
  - Dual-axis electrodes

Supported array configurations:

- Front-only array
- Back-only array
- Full-torso array
  - Uses torso surface normals as the radial direction

#### Customisation Options

You can control:

- Sensor spacing  
  - Default OPM: **30 mm**
- Sensor offset from the body  
  - Default OPM: **10 mm**
  - Surface electrodes: **0 mm**
- Torso coverage exclusion percentages:
  - Top
  - Bottom
  - Left
  - Right

These parameters are demonstrated in the provided example scripts.

---

## Bone Model Variants

A key feature of this toolbox is support for multiple bone geometries:

- Continuous bone
- Homogeneous toroidal
- Inhomogeneous toroidal
- Realistic MRI-segmented bone

Availability:

- **Anatomical model**: all four bone types
- **Canonical model**: all except the realistic MRI-segmented bone

---

## Optional Brain Registration

To investigate **concurrent cortico–spinal interactions**, you may also include a brain model:

- Uses the **SPM template brain**
- Requires selection of three fiducials:
  - Left preauricular
  - Right preauricular
  - Nasion

If you wish to **export the transformation matrix** applied to the SPM brain template:
- Uncomment **line 345** in `cr_check_registration`

---

## Spinal Cord Source Model

The toolbox includes a function to:

- Identify the **centreline of the spinal cord**
- Place candidate source points along this centreline

This step is **optional** and only required if you wish to simulate distributed spinal sources.

---

## Forward Modelling

For boundary element forward modelling (as described in:

> *[Insert forward modelling paper citation here]*  
> *[Insert related GitHub repository here]*

you will need to export:

- All meshes
- Spinal cord source locations
- The transformation matrix

These outputs can then be passed to your BEM/FEM forward modelling pipeline.

---

## Citation

If you use this toolbox in your work, please cite:

> *[Add paper citation here]*

---

## Contact

For questions, issues, or contributions, please open an issue or pull request on GitHub.
