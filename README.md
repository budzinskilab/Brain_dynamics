# Brain Connectome Harmonics

Replication of Vohryzek et al. (2025, PNAS), *"Human brain dynamics are
shaped by rare long-range connections over and above cortical geometry"*
— the EDR+LR (Exponential Distance Rule + Long-Range exceptions)
harmonic-mode framework for reconstructing brain activity from
structural-connectome eigenmodes.

## What's in this repo

- **Root `.m` scripts**: the core replication pipeline — connectome
  construction (`connectome_harmonics.m`) and full-resolution task/FC
  reconstruction (`step4_5_*`, `step5_*`, `step6_*`, `step6b_*`).
- **Figure scripts**: reproductions of the paper's Fig. 2 and Fig. 3
  (`figure_2_fMRI_reconstruction.m`, `figure_2_our_replication.m`,
  `figure_3_tfMRI_reconstruction.m`, `figure_3D_surface_reconstruction.m`,
  `generate_our_figure3A.m`, `generate_our_figure3B.m`), and the
  structural and functional connectivity matrices
  (`generate_connectivity_matrices.m`, `figure_EDRLR_components.m`,
  `generate_functional_connectivity_matrices.m`,
  `figure_fc_longrange_mask.m`).
- **`functions/`**: supporting utilities — eigendecomposition and
  parcellation.
- **`results/`**: our figures and the smaller `.mat` outputs (the
  paper's own plots are included alongside ours for comparison).
- **`Project Update.pptx`**: project update deck (NEUR 4995).

## What's NOT in this repo, and why

- **HCP-derived data** (structural/functional connectomes, task-fMRI
  activation maps, resting-state scans): excluded because the Human
  Connectome Project's Data Use Agreement prohibits redistribution.
  Obtain your own approved access via the WU-Minn HCP consortium, and
  the precomputed connectome/eigenmode files from Pang et al. (2023)'s
  own OSF release.
- **Vendored third-party code**: `pang2023_BrainEigenmodes` (Pang et
  al. 2023's toolbox) and `Vohryzek2024` (this paper's own code
  release) — available at their original sources; not ours to
  redistribute.
- **Copyrighted PDFs** (`Articles/`): the literature referenced during
  this project.
- **Result files over 50 MB** (too large for GitHub), regenerable from
  the scripts here, given the data above:

  | File | Size | Regenerate with |
  |---|---|---|
  | `results/EDR_LR_connectome_full_resolution.mat` | 3.1 GB | Stage 2 of `connectome_harmonics.m` |
  | `results/step4_5_full_resolution_reconstruction.mat` | 85 MB | `step4_5_full_resolution_reconstruction.m` |

## How to run the code

All scripts are MATLAB and use paths relative to the repository root,
so **always run them from the repository root** (`cd` into it first).

### Start here: no data needed

These four scripts only read `.mat` files already included in
`results/`, so they run straight after cloning. They're a quick way to
check your setup and see the main results:

| Script | Makes | Shows |
|---|---|---|
| `figure_2_our_replication.m` | `results/our_figure2_replication.png` | Our Fig. 2: long-range FC reconstruction error vs. number of modes |
| `figure_EDRLR_components.m` | `results/connectivity_EDRLR_components.png` | EDR+LR connectome, EDR baseline (blue) vs. long-range splice (red) |
| `figure_fc_longrange_mask.m` | `results/fc_longrange_mask.png` | The functional long-range connections Fig. 2 is scored on |
| `generate_functional_connectivity_matrices.m` | `results/functional_connectivity_matrices.png` | Empirical vs. reconstructed functional connectivity (redraws from saved results) |

### Full pipeline: needs the data

To recompute everything, first place these folders in the repository
root (they are not in this repo; see above):

- `pang2023_BrainEigenmodes/BrainEigenmodes-main/`: Pang et al.'s
  toolbox (GitHub) plus its data (OSF), including the HCP group
  connectome and task maps. Its `functions_matlab/` also provides
  functions the scripts need (`calc_eigendecomposition`,
  `calc_parcellate`, `ft_read_cifti`).
- `Vohryzek2024/vohryzek2024_EDRLR/`: the paper's code release (surfaces,
  parcellation, `redblue` colormap, and the authors' precomputed results).
- `neurodynamics_data/empirical/`: the 255 subjects' full-resolution HCP
  resting-state scans (`rfMRI_raw/*.dtseries.nii`, ~105 GB) and
  `subject_list_HCP.txt`. Only steps 4 and 6 below need these.

Then run in this order (each step's outputs are the next step's inputs):

| # | Script | Needs | Produces |
|---|---|---|---|
| 1 | `connectome_harmonics.m` | Pang + Vohryzek data | `results/EDR_LR_connectome_full_resolution.mat` (the EDR+LR connectome). Stage 1 is a quick parcel-level check; Stage 2 is the full-resolution build and needs a lot of RAM (`run_stage2.m` runs it unattended with a log). |
| 2 | `step4_5_full_resolution_reconstruction.m` | step 1 | `results/step4_5_full_resolution_reconstruction.mat`: the four eigenmode bases + task reconstruction for 255 subjects |
| 3 | `step5_static_reconstruction.m` | Pang + Vohryzek data | `results/step5_static_reconstruction.mat` (parcel-level version; also provides parcel distances used later) |
| 4 | `step6b_fc_longrange_reconstruction_fullres.m` | steps 2, 3 + resting scans | `results/step6b_fc_longrange_reconstruction_fullres.mat`: the Fig. 2 long-range FC result and its statistics (~1 h) |
| 5 | `generate_connectivity_matrices.m` | step 1 | `results/connectivity_matrices_parc.mat` + the structural connectivity figure |
| 6 | `generate_functional_connectivity_matrices.m` | steps 2, 3, 5 + resting scans | Recomputes `results/functional_connectivity_matrices.mat` if you delete it (~30 min) |

Figures, once the steps above exist:

- `figure_2_our_replication.m` (after step 4), `generate_our_figure3A.m`
  and `generate_our_figure3B.m` (after step 2): our Fig. 2 / Fig. 3.
- `figure_EDRLR_components.m`, `figure_fc_longrange_mask.m`: see above.
- `figure_2_fMRI_reconstruction.m`, `figure_3_tfMRI_reconstruction.m`:
  the paper's Fig. 2 / Fig. 3 from the **authors'** precomputed results
  (only need `Vohryzek2024/`), for comparison with ours.
- `figure_3D_surface_reconstruction.m`: needs
  `results/synthetic_EDRLR_eigenmodes_fsLR_32k-lh_200.mat`, which Stage 2
  of `connectome_harmonics.m` is meant to save but has not produced yet
  (the run stopped at the eigendecomposition).

`step6_fc_longrange_reconstruction.m` is an earlier, parcel-resolution
approximation of step 4, kept for reference.

## Source paper

Vohryzek, J., Sanz-Perl, Y., Kringelbach, M.L., & Deco, G. (2025).
Human brain dynamics are shaped by rare long-range connections over and
above cortical geometry. *PNAS*, 122(1), e2415102122.
