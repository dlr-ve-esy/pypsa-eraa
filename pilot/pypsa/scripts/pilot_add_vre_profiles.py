#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pypsa
import yaml

def read_vre_p_max_pu(network, config):

    TY = snakemake.wildcards[0]
    simulation_year = int(snakemake.wildcards[1])
    short_names = config['carriers_short_names']
    basedir = config['io']['databasedir']+f'scenario_{TY}_{simulation_year}/'

    vre_p_max_pu = pd.read_csv(basedir+f'generation_vre_{simulation_year}.csv', header=0, index_col=0)

    vre_gens = pd.DataFrame(index=vre_p_max_pu.columns)
    vre_gens['carrier'] = [c.split(' - ')[1] for c in vre_p_max_pu.columns]
    vre_gens['carrier_short'] = vre_gens['carrier'].map(short_names)
    vre_gens['node'] = [c.split(' - ')[0] for c in vre_p_max_pu.columns]
    vre_gens['name_short'] = vre_gens[['node', 'carrier_short']].agg(' '.join, axis=1)
    
    vre_p_max_pu.columns = vre_p_max_pu.columns.map(vre_gens['name_short'])
    vre_p_max_pu.index = network.snapshots

    return(vre_p_max_pu)

### add inflow to run of river plants
def read_inflow(network, config):

    TY = snakemake.wildcards[0]
    simulation_year = int(snakemake.wildcards[1])
    short_names = config['carriers_short_names']
    basedir = config['io']['databasedir']+f'scenario_{TY}_{simulation_year}/'
    
    inflow = pd.read_csv(basedir+f'hydro_inflows_{simulation_year}.csv', header=0, index_col=0)

    hydro_gens = pd.DataFrame(index=inflow.columns)
    hydro_gens['carrier'] = [c.split(' - ')[1] for c in inflow.columns]
    hydro_gens['carrier_short'] = hydro_gens['carrier'].map(short_names)
    hydro_gens['node'] = [c.split(' - ')[0] for c in inflow.columns]
    hydro_gens['name_short'] = hydro_gens[['node', 'carrier_short']].agg(' '.join, axis=1)
    
    inflow.columns = inflow.columns.map(hydro_gens['name_short'])
    inflow.index = network.snapshots
    
    inflow = inflow.loc[:, inflow.columns[inflow.columns.isin(network.generators.index)]]
    inflow_pu = (inflow / network.generators.p_nom[inflow.columns]).clip(upper=1.)

    return(inflow_pu)

def add_generators_p_max_pu(network, config):

    vre_p_max_pu = read_vre_p_max_pu(network, config)

    inflow_pu = read_inflow(network, config)    

    p_max_pu = pd.concat([vre_p_max_pu, inflow_pu], axis=1)
    network.generators_t.p_max_pu = p_max_pu[p_max_pu.columns[p_max_pu.columns.isin(network.generators.index)]]

if __name__ == "__main__":

    network = pypsa.Network(snakemake.input[0])
    network.name = 'DestinE Pilot with VRE profiles'

    add_generators_p_max_pu(network, snakemake.config)

    network.export_to_netcdf(snakemake.output[0])




