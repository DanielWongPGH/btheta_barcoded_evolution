## Bar-seq

### Inputs

- Paired R1/R2 FASTQs.
- `fastq_stems_*.txt`: sets of libraries (corresponding to different experiments/sequencing runs). First line is the FASTQ directory and its remaining lines are sample stems. Downloading the fastq files from the SRA may change these filenames, but they should be recoverablwe from the metadata.
- `library_mislabelings.tsv`: a small number of libraries were manually identified as *likely* mislabelings, which can be corrected for during processing.
- Bartender 1.1 and the executables listed in `project_config.yaml`.

The checked-in manifests preserve the sample stems. Generate their first-line
FASTQ directories from `project_config.yaml` before submission:

```bash
bash step0_generate_fastq_stems.sh ../../project_config.yaml
```

The script maps each `fastq_stems_<name>.txt` file to
`raw_processing.barseq.manifest_fastq_dirs.<name>`, falling back to the Bar-seq
`fastq_dir` when no specific key exists. It replaces only the first line and is safe to rerun. Specific manifests
can be supplied after the config path if only a subset should be updated.

### Ordered workflow

1. Configure the pipeline.

   Set `raw_processing.barseq.execution.scheduler` to `slurm` on a cluster or
   `local` to run sequentially. Set `load_modules: false` where the configured
   Python and Bartender executables are already available. Both plain and
   gzip-compressed FASTQs are accepted; compressed raw files are not modified.

   ```bash
   cd raw_processing/process_barseq
   # Edit ../../project_config.yaml, then generate portable manifest headers.
   bash step0_generate_fastq_stems.sh ../../project_config.yaml
   ```

2. Cluster barcodes for every sequencing batch.

   ```bash
   bash step1_sbatch_bartender.sh ../../project_config.yaml fastq_stems_E1_barseq.txt -UMI false
   bash step1_sbatch_bartender.sh ../../project_config.yaml fastq_stems_E2_barseq.txt -UMI false
   bash step1_sbatch_bartender.sh ../../project_config.yaml fastq_stems_EV_barseq.txt -UMI false
   bash step1_sbatch_bartender.sh ../../project_config.yaml fastq_stems_rebarseq.txt -UMI true
   ```

   `-UMI false` runs `step1b_run_bartender_no_umi.sh` directly on R1 and writes
   `bartender_noUMI/*_cluster.csv`; this is the required mode for E1, E2, and
   the in-vitro libraries. `-UMI true` runs
   `step1a_run_bartender_with_umi.sh`, which merges R1 with the R2 UMI before
   clustering and writes under `bartender/`; this is the required rebarseq
   mode. Omitting `-UMI` defaults to `false`, while `-UMI` without a value is
   shorthand for `-UMI true`.

3. Wait for all Bartender cluster files to finish, then estimate noise for the
   E1, E2, and in-vitro libraries.

   ```bash
   bash step3_sbatch_noise_estimation.sh ../../project_config.yaml fastq_stems_E1_barseq.txt
   bash step3_sbatch_noise_estimation.sh ../../project_config.yaml fastq_stems_E2_barseq.txt
   bash step3_sbatch_noise_estimation.sh ../../project_config.yaml fastq_stems_EV_barseq.txt
   ```

   This stage reads the paired FASTQs and no-UMI clusters and writes
   `pseudo_derep/*_noise.out` plus diagnostic PDFs.

4. After every noise job completes, aggregate the library counts and effective
   depths. The aggregator expects the historical batch layout beneath one data
   root: `E1_barseq`, `E2_barseq`, `vitro_barseq`, and `rebarseq`.

   ```bash
   python make_barcode_counts_csv.py DATA_ROOT \
       --noise-estimate-dir NOISE_ROOT \
       --output-dir OUTPUT_DIR
   ```

   Its principal output is `all_barcode_reads.csv` (with trimmed/demo variants
   also written). E1, E2, and in-vitro effective depths come from the noise
   results; rebarseq counts come from the UMI-aware Bartender output.

5. Optionally assemble the retained noise-inference pickle.

   ```bash
   python process_noise_estimates.py NOISE_DIR \
       --output pickled_noise_inference.pkl
   ```

   The aggregator checks that each filename agrees with the library ID stored
   inside the file and rejects duplicate library IDs. It does not require the
   legacy `mouse_col_map.pkl` or `vitro_col_map.pkl` intermediates.

6. Transfer the count table into the local repository at
   `data/rebarseq/all_barcode_reads.csv`, or the corresponding path configured
   in the root `project_config.yaml`.

7. From the repository root, build the current notebook inputs.

   ```bash
   python raw_processing/process_barseq/make_barcode_arrays.py
   ```

   This is the authoritative array builder. It reads the configured
   `local.rebarseq` directory and writes the barcode arrays,
   depth arrays, sample metadata, column maps, and overlap maps under
   `local.pickles`. The historical `notebooks/make_barcode_arrays.py` command is
   retained as a compatibility entry point to the same implementation.

8. Run `notebooks/0a_barcode_sorting.ipynb` for barcode pool assignments and
   index-hopping diagnostics. Run `notebooks/estimate_noise_floors.ipynb` where
   its noise-floor products are required, then proceed to the numbered analysis
   notebooks.

### Bar-seq completion checks

- Every requested sample has its expected `*_cluster.csv`.
- Every non-rebarseq library has a corresponding `*_noise.out`.
- Rebarseq clusters are in the UMI-aware `bartender/` tree.
- `all_barcode_reads.csv` has been copied to the configured local input path.
- `raw_processing/process_barseq/make_barcode_arrays.py` completes and refreshes
  the configured `local.pickles` directory.

### Bar-seq scripts

| Script | Role |
|---|---|
| `step1_sbatch_bartender.sh` | Select UMI/no-UMI mode and submit barcode-clustering jobs |
| `step1a_run_bartender_with_umi.sh` | Merge paired reads with R2 UMIs and cluster them |
| `step1b_run_bartender_no_umi.sh` | Cluster R1 barcodes without UMI handling |
| `merge_barseq_fastq.py` | Construct the UMI-aware merged FASTQ |
| `step3_sbatch_noise_estimation.sh`, `step3_srun_noise_estimate.sh` | Submit and run per-library noise fits |
| `estimate_noise.py` | Fit the technical-noise model |
| `process_noise_estimates.py` | Assemble retained noise estimates into a pickle |
| `make_barcode_counts_csv.py` | Build the cross-library barcode-count table |
| `step0_generate_fastq_stems.sh` | Populate manifest FASTQ directories from `project_config.yaml` |
| `../parse_config.sh` | Shared loader for the Bar-seq YAML subsection |
