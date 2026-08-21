# B. thetaiotaomicron Barcoded Evolution

Code, processed data, and analysis notebooks associated for the paper: "High-resolution lineage tracking of within-host evolution and strain transmission in a human gut symbiont across ecological scales" ([bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.02.17.580834v1.abstract))

Associated raw barcode, isolate whole genome, and population/metagenomic sequencing data is deposited on the [SRA: BioProject accession PRJNA1111532](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1111532).

## Repository structure

- `notebooks/` — analysis and figure-generation notebooks
- `data/` — processed data and intermediate analysis files
- `figures/` — generated figure files
- `papers/` — manuscript, supplementary information, and figure/code crosswalk
- `remote_processing/` — preprocessing workflows for barcode sequencing, isolate sequencing, and metagenomics. This processing takes a ~day with access to ~dozens of cores.

## Environment

The local analysis environment uses Python 3.9 and is described in `environment.yaml`.

Create the environment with Conda:

```bash
conda env create -f environment.yaml
conda activate btheta-barcoded-evolution
```

If the environment already exists, update it with:

```bash
conda env update -f environment.yaml
```

Restart the Jupyter kernel after changing package versions.

## Configuration

Project paths are defined in `project_config.yaml`. Update the system-specific entries before running analyses.

Python code can load the configuration through `project_config.py`. Bash workflows can load the same configuration with:

```bash
source ./load_project_config.sh
```

## Running analyses

Start Jupyter from the repository root:

```bash
jupyter lab
```

The notebooks are organized approximately by manuscript figure. More detailed execution order, input dependencies, and expected outputs will be documented here as the workflow is consolidated.

## Data processing

The repository separates upstream processing from downstream analysis:

- `remote_processing/process_barseq/` — barcode sequencing
- `remote_processing/process_isolates/` — isolate sequencing
- `remote_processing/process_metagenomics/` — metagenomic sequencing
- `notebooks/` — downstream calculations, analyses, and plotting

Some workflows may require data generated on a computing cluster or copied from archival storage. These dependencies will be documented as they are verified.

## Reproducibility status

The manuscript-to-code audit is maintained in `papers/FIGURE_CALCULATION_CODE_CROSSWALK.md`. It records the notebook or script associated with each figure, relevant inputs, and current execution status.

## Citation

Citation information will be added when the manuscript record is finalized.

## License

License information has not yet been specified.
