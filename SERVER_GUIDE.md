# Running the project on parana

A short guide to running this project from your own account on the lab
server (`parana`). The code comes from GitHub; the data (~120 GB) stays in
Gop's project folder and is shared read-only, so nothing large is copied.

## 0. One-time setup (Gop does this)

Gop gives your account read access to the project, replacing `User` with
your parana username:

```bash
setfacl -m u:User:x ~
setfacl -R -m u:User:rX ~/brain_connectome_harmonics
setfacl -R -d -m u:User:rX ~/brain_connectome_harmonics
```

## 1. Get the code

```bash
ssh User@parana
git clone https://github.com/budzinskilab/Brain_dynamics.git ~/brain_dynamics
cd ~/brain_dynamics
```

## 2. Link the data

These point to Gop's copies instead of duplicating them:

```bash
G=/home/gop/brain_connectome_harmonics
ln -s $G/pang2023_BrainEigenmodes pang2023_BrainEigenmodes
ln -s $G/Vohryzek2024 Vohryzek2024
ln -s $G/neurodynamics_data neurodynamics_data
```

Optional: to skip the two slowest steps, also copy their outputs
(about 3.3 GB):

```bash
cp $G/results/EDR_LR_connectome_full_resolution.mat \
   $G/results/step4_5_full_resolution_reconstruction.mat \
   $G/results/synthetic_EDRLR_eigenmodes_fsLR_32k-lh_200.mat results/
```

## 3. Quick check (about a minute)

```bash
matlab -batch "figure_2_our_replication"
```

This should finish without errors and write
`results/our_figure2_replication.png`, the replication of the paper's
Fig. 2.

## 4. Full run (several hours)

Runs every step of the pipeline and then all the figures, in the
background so it keeps going after you log out:

```bash
cd ~/brain_dynamics
nohup matlab -batch "run_pipeline" < /dev/null > run_pipeline_console.txt 2>&1 &
disown
```

Check on it:

```bash
grep "==========" run_pipeline_log.txt     # which step it's on (14 in total)
tail -n 5 run_pipeline_log.txt            # latest output
sed -n '/Summary/,$p' run_pipeline_log.txt   # summary, once finished
```

If a pipeline step fails, the run stops there and the summary says which
one. Everything is written to your own `~/brain_dynamics/results/`;
Gop's files are never changed.

To run individual scripts or figures instead, see "How to run the code"
in the README for the order and what each one needs.

## Notes

- MATLAB needs the Statistics and Machine Learning, Signal Processing and
  Image Processing toolboxes.
- The resting-state scans are HCP data, covered by the Human Connectome
  Project's data use terms.
- To get the latest code later: `cd ~/brain_dynamics && git pull`.
