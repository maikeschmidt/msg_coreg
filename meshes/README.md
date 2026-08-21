# meshes/

The STL surface meshes shipped with msg_coreg. All are loaded by name
through `cr_load_meshes`; you should not normally need to open them
directly.

Naming follows `<space>_<extent>_<variant>.stl`:

| Part | Values | Meaning |
|---|---|---|
| space | `canonical`, `mri` | canonical (seated template) or MRI-derived anatomy |
| extent | `full`, `cervical` | whole torso, or cervical region only |
| variant | `cont`, `homo`, `inhomo` | continuous, homogeneous toroidal, inhomogeneous toroidal bone |

| File group | Contents |
|---|---|
| `canonical_torso.stl`, `mri_torso.stl` | outer body boundary (includes head and neck) |
| `canonical_heart.stl`, `heart.stl`, `mri_lungs.stl`, `canonical_lungs.stl` | thoracic organs |
| `*_full_*.stl`, `*_cervical_*.stl` | bone variants, per space and extent |
| `realistic_full_bone.stl`, `realistic_cervical_bone.stl` | MRI-segmented realistic vertebrae |
| `spine.stl`, `cervical_spine.stl`, `mri_full_spine.stl`, `mri_cervical_spine.stl` | spinal cord |
| `back_muscle_temp.stl`, `vagus_nerve_temp.stl` | optional extra structures |

The scanner-cast optical surface (`surface.stl`) used by the anatomical
workflow is **not** bundled here — supply your own for your setup.
