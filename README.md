# Destination Earth - Use Case Energy System: Pilot Implementation

This is the pilot implementation of the REMix model instance in the Destination
Earth - Use Case Energy System contract. The pilot consists of two parts:

1. data collection and preprocessing
2. model building and solving

Both parts are implemented in a respective automatic workflow using
[Snakemake](https://snakemake.readthedocs.io/en/stable/). Exemplary evaluation
has been carried out based on 10 scenarios.

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
2007, 2010, 2012. To change the target years and climatic years, adjust the
respective `config.yaml` file inside the `ENTSOE`.

### Data collection and Preparation

The data collection automatically downloads the ERAA 2022 datasets and extracts
the respective relevant scenario data from them. To run the data collection and
preparation workflow use

    cd ENTSOE
    snakemake --cores 1 download
    snakemake --cores 1 extract_ntc extract_demand extract_climate_data
    snakemake --cores 10

![Preprocessing workflow](ENTSOE/workflow.png)

### Model building and solving

The model building reads the prepared data and reformats them in a way that is
compatible with the REMix energy system optimization software. After that, the
datasets for the different energy system components (i.e. storage, renewable
generators, AC and DC lines, etc.) are collected in a singe REMix instance and
a dispatch optimization is carried out. The workflow prepares the output data
for the adequacy assessment.

    cd ../remix
    snakemake --cores 10

![REMix workflow](ENTSOE/workflow.png)

## Adequacy assessment

The adequacy assessment is not part of the automatic workflow. Exemplary results
of the adequacy assessment for the scenarios are available in
[this](remix/results.ipynb) jupyter notebook.
