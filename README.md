# Destination Earth - Use Case Energy System: Pilot Implementation

This is the pilot implementation of the REMix model instance in the Destination
Earth - Use Case Energy System contract. The pilot consists of two parts:

1. data collection and preprocessing
2. model building and solving

Both parts are implemented in a respective automatic workflow using
[Snakemake](https://snakemake.readthedocs.io/en/stable/).

## Requirements

- Python 3.10
- GAMS 41.5
- [Python third party packages](requirements.txt)

Set up a python environment and GAMS, then clone or download this repository.
From the root of the repository create a new environment and install the
required python packages:

    conda create -n destine-env python=3.10 -y
    conda activate destine-env
    pip install -r requirements.txt

## Execute the workflows

To run the workflows, first run the ENTSOE workflow, then run the remix
workflow. The workflow for the pilot is configured in a way, that 10 scenarios
are executed, i.e. target years 2027 and 2030 and climatic years 1997, 2000,
2007, 2010, 2012.

    cd ENTSOE
    snakemake --cores 10
    cd ../remix
    snakemake --cores 10

Change the respective `config.yaml` file inside the `ENTSOE` folder to change
the scenarios to be optimized.

## Adequacy assessment

Exemplary results of the adequacy assessment for the scenarios are available in
[this](remix/results.ipynb) jupyter notebook.
