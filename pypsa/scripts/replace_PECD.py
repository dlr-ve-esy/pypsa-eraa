#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import pypsa

from pypsa.clustering.spatial import get_clustering_from_busmap
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

def simplify_network(network):

    #### cluster network
    stores_strategies = {
        'e_max_pu': "capacity_weighted_average",
        'e_min_pu': "capacity_weighted_average",
        'e_initial': "sum",
        'e_nom': 'sum'    
    }

    network.loads.carrier = 'total_demand'
    
    busmap = {bus: bus[:2] for bus in network.buses.index}
    busmap['DEKF'] = 'KF' # Kriegers Flak
    busmap['DKKF'] = 'KF'
    busmap['LUV1'] = 'DE' # pumped storage Vianden (located in Luxenbourg but connected to the German grid)

    stores = network.stores
    for store in stores.index:
        busmap[stores.bus[store]] = f'{stores.bus[store][:2]} {stores.carrier[store]}'

    Clustering = get_clustering_from_busmap(network, busmap, with_time=True, aggregate_one_ports={'Generator', 'Load', 'Store'}, one_port_strategies={'Store': stores_strategies})

    n_clustered = Clustering.network
    
    n_clustered.links_t.p_max_pu.drop(n_clustered.links_t.p_max_pu.columns[~n_clustered.links_t.p_max_pu.columns.isin(n_clustered.links.index)], axis=1, inplace=True)
    n_clustered.links_t.p_min_pu.drop(n_clustered.links_t.p_min_pu.columns[~n_clustered.links_t.p_min_pu.columns.isin(n_clustered.links.index)], axis=1, inplace=True)

    return(n_clustered)

if __name__ == "__main__":

    network = pypsa.Network(snakemake.input['network'])
    replace_with = snakemake.wildcards['ALT'] #snakemake.output[0].split('/')[-1].split('_')[0]
    CY = snakemake.wildcards['simulation_year'] #snakemake.output[0].split('_')[-1].split('.')[0]

    ### replace onwind p_max_pu
    data_basedir = snakemake.input['cf_basedir']
    dat_f = data_basedir+f'/{replace_with}.csv'

    p_max_pu = pd.read_csv(dat_f, header=0, index_col=0, parse_dates=True)
    p_max_pu = p_max_pu.loc[network.snapshots, :]
    p_max_pu.columns = [c+' onwind' for c in p_max_pu.columns]

    p_max_pu.drop(p_max_pu.columns[~p_max_pu.columns.isin(network.generators_t.p_max_pu.columns)], axis=1, inplace=True)
    
    network.generators_t.p_max_pu.loc[:,p_max_pu.columns] = p_max_pu#.loc[network.snapshots,:]

    ### aggregate to country scale
    n_clustered = simplify_network(network)
    n_clustered.name = f'DestinE pilot with {replace_with}'
    
    n_clustered.export_to_netcdf(snakemake.output[0])