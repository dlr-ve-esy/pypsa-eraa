#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pypsa

def map_column_names(p_max_pu, short_names):
    vre_gens = pd.DataFrame(index=p_max_pu.columns)
    vre_gens['carrier'] = [c.split(' - ')[1] for c in p_max_pu.columns]
    vre_gens['carrier_short'] = vre_gens['carrier'].map(short_names)
    vre_gens['node'] = [c.split(' - ')[0] for c in p_max_pu.columns]
    vre_gens.loc[vre_gens.node == 'GB', 'node'] = 'UK'
    vre_gens['name_short'] = vre_gens[['node', 'carrier_short']].agg(' '.join, axis=1)
    
    p_max_pu.columns = p_max_pu.columns.map(vre_gens['name_short'])

    return(p_max_pu)

if __name__ == "__main__":

    network = pypsa.Network(snakemake.input['network'])
    replace_with = snakemake.wildcards['CDS'] #snakemake.output[0].split('/')[-1].split('_')[0]
    CY = snakemake.wildcards['simulation_year'] #snakemake.output[0].split('_')[-1].split('.')[0]

    network.name = f'DestinE pilot with {replace_with}'
    
    data_basedir = snakemake.input['cf_basedir']
    dat_f = data_basedir+f'/{replace_with}/CY{CY}/generation_vre_timeseries.csv'

    short_names = snakemake.config['carriers_short_names']

    p_max_pu = pd.read_csv(dat_f, header=0, index_col=0, parse_dates=True)
    p_max_pu = map_column_names(p_max_pu, short_names)

    if len(p_max_pu.index) > 8760:
        p_max_pu = p_max_pu.iloc[:8760,:]
        p_max_pu.index = network.snapshots

    p_max_pu.drop(p_max_pu.columns[~p_max_pu.columns.isin(network.generators_t.p_max_pu.columns)], axis=1, inplace=True)
    
    network.generators_t.p_max_pu.loc[:,p_max_pu.columns] = p_max_pu#.loc[network.snapshots,:]
    network.export_to_netcdf(snakemake.output[0])