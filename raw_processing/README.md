# Raw processing workflows

The workflows below process the raw fastq libraries ([PRJNA1111532](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1111532)) to generate cleaned data used in notebooks. SLURM-scheduling scripts used on the Stanford Sherlock cluster are included (which should be adjusted for other systems).

## Overview

| Workflow | Raw input | Main output | 
|---|---|---|
| [Barcode-seq](#barcode-seq) | paired Bar-seq FASTQs | barcode clusters, noise estimates, `all_barcode_reads.csv` |
| [Isolate sequencing](#isolate-sequencing) | paired isolate FASTQs | breseq GD files and BAMs |
| [Metagenomics](#metagenomics) | paired metagenomic FASTQs | merged mutation timecourses and coverage timecourses | 


## Configuration

All local and raw-processing settings live in the root `project_config.yaml`. The shared `raw_processing/parse_config.sh` loader reads only the requested `raw_processing.<pipeline>` subsection. 

| Pipeline | Software (environment)| 
|---|---|
| `process_barseq` | Python 3.9+, NumPy, SciPy, Matplotlib, Bartender 1.1, `yq` | 
| `process_isolates` | Python 3, Java, Trimmomatic, breseq, `yq` | 
| `process_metagenomics` | Python 2 environment below, Trimmomatic, breseq, SAMtools, gdtools, `yq` | 

Before running a workflow:

1. Replace every placeholder path, executable, account, partition, and module
   under `raw_processing` in `project_config.yaml`.
2. Generate the FASTQ-manifest headers as described by each workflow. The
   metagenomics `population_samples.csv` is generated from its E1/E2 manifests
   by `make_timecourse_fastq_lists.py`.
3. Invoke the workflow scripts from either the repository root or their own
   pipeline directory; bundled scripts and metadata are resolved relative to
   the script location.
4. Pass `../../project_config.yaml` to every orchestration command.
5. Keep generated reads, logs, and results outside this source directory.

For example:

```bash
cd raw_processing/process_barseq
# Edit ../../project_config.yaml for this system before submitting jobs.
```

The shell wrappers require the mikefarah implementation of `yq`. Install it
with Conda, Homebrew, or the package method appropriate to the cluster, and
confirm that the configured module names exist before submission. If a job
fails immediately with a module or executable error, fix `modules.*` or the
corresponding executable path in `project_config.yaml`. Bartender configuration must
point to the directory containing `bartender_extractor_com` and
`bartender_single_com`.

The Bar-seq and isolate workflows use Python 3. The historical metagenomic
parsers require Python 2 and their frozen NumPy/SciPy versions:

```bash
conda env create -f raw_processing/environment-python2.yaml
conda activate btheta-remote-python2
```

Set `raw_processing.metagenomics.python2.executable: python` in
`project_config.yaml` after activating that environment. Porting those parsers to Python 3 should only be
done with output-level regression fixtures because they define the published
parsing and fuzzy junction-merging behavior.

## Isolate sequencing

### Inputs

- Paired isolate FASTQs.
- `clone_barcode.tsv`, which maps sequencing samples to clones/barcodes.
- The reference genome and adapter/executable paths in `project_config.yaml`.

### Ordered workflow

1. Configure the pipeline and update `clone_barcode.tsv` if necessary.

   ```bash
   cd raw_processing/process_isolates
   # Edit ../../project_config.yaml.
   bash generate_fastq_stems.sh ../../project_config.yaml
   ```

   The generator updates the directory header in `E1_isolates.txt` from
   `raw_processing.isolates.manifest_fastq_dirs.E1_isolates`. Step 0 also
   discovers current samples directly from its FASTQ-directory argument and
   writes `samples_in_dir.txt`; the retained manifest is useful when starting
   at step 1 or checking the historical sample set.

2. Discover samples and trim adapters.

   ```bash
   bash step0_trim_adapters.sh ../../project_config.yaml /path/to/untrimmed_fastq
   ```

   This creates `samples_in_dir.txt` and submits one trimming job per sample.
   The trimmed FASTQs are written beneath the configured `trimmed` output.
   Run from this directory because the sample parser reads `clone_barcode.tsv`
   relative to the working directory.

3. Wait for every trimmed read pair, then run breseq.

   ```bash
   bash step1_breseq_genomes.sh ../../project_config.yaml samples_in_dir.txt
   ```

   Each sample produces a breseq result directory containing, among other
   files, `output/evidence/evidence.gd` and `data/reference.bam`.

4. Optionally run the auxiliary Delly structural-variant workflow after the
   BAMs exist.

   ```bash
   bash call_structural_variants_with_delly.sh ../../project_config.yaml samples_in_dir.txt
   ```

   Delly outputs are not an identified prerequisite of the current manuscript
   notebooks; those notebooks use breseq evidence, including JC records. Treat
   this as an optional parallel analysis unless its downstream role is restored.

5. Export the clone map, breseq GD/evidence files, and required BAM/coverage
   products into the local paths configured under `data/isolate_wgs/` and
   `data/BAM_coverages/`.

6. Run `notebooks/4_isolate_sequencing.ipynb` to parse isolate mutations and
   generate the mutation/driver pickles and supplementary tables. Notebooks
   `4a_isolate_seq_plotting.ipynb` and `4c_BAM_coverages.ipynb` consume those
   processed calls and coverage products for the manuscript figures.

### Unresolved isolate handoff

The curated raw scripts end at ordinary breseq outputs. The exact historical
step that annotated or renamed GD files and copied them into the local
`data/isolate_wgs/` layout was not found. Until that export is scripted, compare
the expected filenames in notebook `4` with the generated per-sample breseq
directories before replacing the checked-in data.

### Isolate scripts

| Script | Role |
|---|---|
| `step0_trim_adapters.sh`, `step0_trim_adapters.sbatch` | Discover samples and trim paired reads with Trimmomatic |
| `step1_breseq_genomes.sh`, `step1_breseq_genomes.sbatch` | Submit and run per-isolate breseq calls |
| `call_structural_variants_with_delly.sh`, `call_structural_variants_with_delly.sbatch` | Optional Delly structural-variant calls |
| `parse_fastq_directory.py` | Build the sample list from FASTQs and `clone_barcode.tsv` |
| `generate_fastq_stems.sh` | Populate isolate manifest directories from `project_config.yaml` |
| `../parse_config.sh` | Shared loader for the isolate YAML subsection |

## Metagenomics

Activate the Python 2 environment described above before running this
workflow.

### Inputs

- Paired metagenomic FASTQs.
- `population_samples.csv`, which assigns samples and time points to
  populations.
- `all_mgx_timecourse_lists.tsv` and the relevant FASTQ-stem lists.
- Reference genome, BED regions, and executable paths in `project_config.yaml`.

`population_samples.csv` is generated by `make_timecourse_fastq_lists.py` from
the E1/E2 manifests and includes their configured directory headers. See
`process_metagenomics/README.md` for generation and portability details.

### Ordered workflow

For steps that take `POPULATION`, repeat the command for each desired population
or supported comma-separated group. `SEQUENCE_TYPE` controls the metadata
subset; retained examples include `non_stragglers`, `stragglers`, and
`non_straggler_non_clones`.

1. Configure the pipeline and trim the reads.

   ```bash
   cd raw_processing/process_metagenomics
   # Edit ../../project_config.yaml.
   bash generate_fastq_stems.sh ../../project_config.yaml
   bash step0_trim_adapters.sh ../../project_config.yaml /path/to/untrimmed_fastq
   ```

   The generator updates the E1 and E2 FASTQ-manifest headers from
   `raw_processing.metagenomics.manifest_fastq_dirs` without changing their
   sample stems.

   This generates `samples_in_dir.txt` and submits the trimming jobs.

2. Wait for trimming, then call each sample with breseq in polymorphism mode.

   ```bash
   bash step1_breseq_genomes.sh ../../project_config.yaml POPULATION SEQUENCE_TYPE
   ```

3. Build the union of candidate junction evidence for each population.

   ```bash
   bash step2_merge_candidate_junctions.sh ../../project_config.yaml POPULATION SEQUENCE_TYPE
   ```

   This writes `${POPULATION}_merged_breseq_output.gd`. However, the current
   stage-3 wrapper rebuilds the same junction union internally. Stage 2 is thus
   useful as an inspectable checkpoint but is not a strict prerequisite when
   stage 3 is run unchanged.

4. Re-run breseq for all time points using the population-wide junction union.

   ```bash
   bash step3_rebreseq_genomes.sh ../../project_config.yaml POPULATION SEQUENCE_TYPE
   ```

5. Wait for every rebreseq GD and BAM. Then launch two independent branches;
   they may run in parallel.

   Mutation-frequency branch:

   ```bash
   bash step4_create_distributed_timecourse.sh ../../project_config.yaml POPULATION SEQUENCE_TYPE
   ```

   This submits distributed `samtools mpileup`/parsing jobs and writes numbered
   `${POPULATION}_*_timecourse.txt` chunks.

   Coverage branch:

   ```bash
   bash step4b_make_coverage_timecourse.sh ../../project_config.yaml POPULATION SEQUENCE_TYPE
   ```

   This writes `${POPULATION}_coverage_timecourse.txt`, including read counts,
   genome-wide coverage, and depth in configured BED regions.

6. After every mutation-timecourse chunk has completed, merge and compress it.

   ```bash
   bash step5_create_merged_timecourse.sh ../../project_config.yaml POPULATION SEQUENCE_TYPE
   ```

   This combines the chunks with rebreseq evidence, performs the historical
   fuzzy junction merge, and writes `${POPULATION}_merged_timecourse.bz2`.

7. Export/decompress the products into the local configured layout:

   - merged mutation timecourses feed `data/mgx_seq/*_merged_timecourse` and are
     consumed by `4b` and `invivo_metagenomics`;
   - coverage timecourses feed `data/breseq_timecourse_files/`;
   - BAM-derived coverage summaries used by `4c` live under
     `data/BAM_coverages/`.

### Unresolved metagenomic handoff

The exact historical copying, decompression, and renaming commands between the
remote outputs and those three local directories were not retained. In
particular, `4c_BAM_coverages.ipynb` reads `all_bam_coverage.txt`, which is not
directly produced under that name by `step4b_make_coverage_timecourse.sh`.
Validate sample ordering and expected filenames against the notebooks before
replacing the checked-in local files.

## Troubleshooting

- If `yq` is not found, install the mikefarah implementation and ensure it is
  available inside submitted jobs as well as the login shell.
- If cluster modules fail to load, update `modules.*` in the relevant pipeline
  configuration.
- If isolate FASTQs are not discovered, check the expected
  `{sample}_R1_001.fastq.gz`/`{sample}_R2_001.fastq.gz` naming and
  `clone_barcode.tsv` entries.
- If Delly cannot find its reference, confirm that the preceding breseq job
  completed and produced `data/reference.fasta` and `data/reference.bam`.
- If a downstream stage sees missing files, first check whether the preceding
  SLURM jobs have actually completed; the orchestration commands are
  asynchronous.

## Software citations

The processing workflows use Bartender (Zhao et al., 2018), breseq (Deatherage
and Barrick, 2014), Delly (Rausch et al., 2012), and Trimmomatic (Bolger et al.,
2014). Cite the tools used in a particular rerun and record their exact
versions with the output.

## Snapshot classification

| Snapshot material | Classification | Disposition |
|---|---|---|
| `process_barseq_clean/` | newest Bar-seq orchestration; canonical source for the port | essential scripts and manifests ported |
| `process_barseq/` | same core algorithms; older hard-coded wrappers plus archived scripts | redundant except as provenance |
| top-level `bartender-1.1-master/` and the two pipeline copies | three identical vendor trees | not copied; install Bartender 1.1 externally and record its source/version |
| `process_isolates/` | newest isolate workflow | trimming and breseq stages ported |
| `isolates_scripts/` | older duplicate of isolate workflow | redundant |
| isolate Delly scripts | auxiliary SV workflow with no identified downstream manuscript handoff | ported with the isolate workflow but marked auxiliary |
| `process_metagenomics/` | newest metagenomic workflow | essential processing and sample metadata ported |
| `metagenomics_scripts/` | mostly byte-identical predecessors | redundant |
| `Kim_todo_submissions/cluster_scripts/` | still earlier mixed isolate/metagenomic generation | obsolete/redundant |
| `prep_fastq_for_SRA/` | archive-deposition utilities | useful for submission provenance, not analysis processing |
| `invivo_metagenomics.ipynb` | downstream analysis/plotting copy | excluded; maintained under local `notebooks/` |
| `barseq_noise/`, `slurm_out/`, `logs/`, pipeline `data/`, CSV/pickle outputs | generated artifacts | excluded from source port |
| one-off relabel/refactor/comparison scripts | migration or exploratory utilities | excluded; essential relabel map retained |

The immutable snapshot remains the source of truth for historical files that
were intentionally excluded.
